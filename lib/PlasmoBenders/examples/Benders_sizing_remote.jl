using Revise
using Random, Distributions, Distributed
t = 0:1:30

Random.seed!(10)
cost1 = sin.(t .* 2 .* pi ./ 10 ) .* 2 .+ 10 .+ rand(Uniform(-1,1),31) .* .5
demand1 = sin.(t .* 2 .* pi ./ 15 ) .* 30 .+ 100 .+ rand(Uniform(-1,1),31) .* 5

using Plasmo, HiGHS, PlasmoBenders

if nprocs() == 1
    addprocs(1)
end

@everywhere begin
    using Revise
    using Random, Distributions, Plasmo, HiGHS, PlasmoBenders
end

# Define graph and add a node
g_root = RemoteOptiGraph()
@optinode(g_root, n_root)

# Set the possible storage and reactor sizes to choose from
storage_sizes = [10, 50, 100, 500, 1000]
reactor_sizes = [10, 20, 50, 100, 200]

# Define variables for the sizes and define binary variables for choosing sizes
@variable(n_root, storage_size >= 0)
@variable(n_root, reactor_size >= 0)
@variable(n_root, storage_bin[1:5], Bin)
@variable(n_root, reactor_bin[1:5], Bin)

# Ensure that we choose one and only one size
@constraint(n_root, sum(storage_bin) == 1)
@constraint(n_root, sum(reactor_bin) == 1)

# Set the size variables equal to the binary variables times their respective options
@constraint(n_root, storage_size == sum(storage_bin[i] * storage_sizes[i] for i in 1:5))
@constraint(n_root, reactor_size == sum(reactor_bin[i] * reactor_sizes[i] for i in 1:5))

# Set the objective on the node
@objective(n_root, Min, 5 * storage_size + 20 * reactor_size)

# Set the objective on the graph
set_to_node_objectives(g_root)

function build_scenario_graph(cost, demand, num_scenarios = 1)
    # Define time points
    T = length(cost)
    
    # Define OptiGraph and nodes (one node for each time point)
    g = RemoteOptiGraph(worker=2)
    @optinode(g, n[1:T])
    # Loop through the nodes  and add variables and constraints
    for (t, node) in enumerate(n)

        # Define variables on the nodes
        @variable(node, 0 <= x_buy )
        @variable(node, -50 <= x_save <= 50)
        @variable(node, 0 <= x_store)
        @variable(node, 0 <= y_product)
        @variable(node, 0 <= unmet_demand)

        # Add constraints to the nodes
        @constraint(node, unmet_demand >= demand[t] - y_product)
        @constraint(node, y_product == 5 * (x_buy - x_save))

        # Define objective on the node
        @objective(node, Min, (unmet_demand * 1000 + x_buy * cost[t]) / num_scenarios)
    end

    # Set the initial storage amount
    @constraint(n[1], n[1][:x_store] == 10)

    # Link the storage variables acorss time points
    @linkconstraint(g, [t = 1:(T-1)], n[t + 1][:x_store] - n[t][:x_store] == n[t][:x_save])
    
    # Set graph objective to summation of node objectives
    set_to_node_objectives(g)

    # Return the graph
    return g
end

# Define overall graph
g = RemoteOptiGraph()

# Define operations level subgraph
g1 = build_scenario_graph(cost1, demand1)

# Add subgraphs to graph `g`
add_subgraph!(g, g_root)
add_subgraph!(g, g1)

# Define time points
T = length(all_nodes(g1))
g1_nodes = all_nodes(g1)

# Define constraints on storage and reactor sizes
l1 = @linkconstraint(g, [t = 1:T], g1_nodes[t][:x_store] <= g_root[:n_root][:storage_size])
l2 = @linkconstraint(g, [t = 1:T], g1_nodes[t][:x_buy] - g1_nodes[t][:x_save] <= g_root[:n_root][:reactor_size])

# Define a subproblem object to use for the subgraphs
solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)

# Define the BendersAlgorithm object and set the subproblem solver
benders_alg = BendersAlgorithm(g, g_root, solver = solver)

run_algorithm!(benders_alg)
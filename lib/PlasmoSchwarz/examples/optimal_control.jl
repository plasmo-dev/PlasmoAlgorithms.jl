#Example demonstrating the use of overlap to solve a long horizon control problem
using Plasmo, Ipopt
using PlasmoSchwarz

T = 200             # number of time points
d = sin.(1:T)       # a disturbance vector
imbalance = 0.1     # partition imbalance
distance = 2        # expansion distance
n_parts = 10        # number of partitions

# create the optigraph
graph = Plasmo.OptiGraph()
@optinode(graph, state[1:T])
@optinode(graph, control[1:(T - 1)])
for (i, node) in enumerate(state)
    @variable(node, x)
    @constraint(node, x >= 0)
    @objective(node, Min, x^2) #- 2*x*d[i])
end
for node in control
    @variable(node, u)
    @constraint(node, u >= -1000)
    @objective(node, Min, u^2)
end

# initial condition
n1 = state[1]
@constraint(n1, n1[:x] == 0)

# dynamics
for i in 1:(T - 1)
    @linkconstraint(graph, state[i + 1][:x] == state[i][:x] + control[i][:u] + d[i])
end

# subproblem optimizer
sub_optimizer = Plasmo.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

# optimize using overlapping schwarz decomposition
optimizer = SchwarzAlgorithm(
    graph;
    n_partitions=n_parts,
    overlap_distance=1,
    subproblem_optimizer=sub_optimizer,
    max_iterations=100,
    mu=10.0,
)

# run the algorithm
run_algorithm!(optimizer)

# check termination status
@show Plasmo.termination_status(optimizer)

# check objective value
@show Plasmo.objective_value(optimizer)

# check first state and control values
@show Plasmo.value(optimizer, state[1][:x])
@show Plasmo.value(optimizer, control[1][:u])

# you can also access the primal and dual feasibility vectors using
# algorithm methods.
prf = PlasmoSchwarz.calculate_primal_feasibility(optimizer)
duf = PlasmoSchwarz.calculate_dual_feasibility(optimizer)

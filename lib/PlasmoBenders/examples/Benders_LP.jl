using Revise
using Plasmo, JuMP, HiGHS, Plots, PlasmoBenders

g = OptiGraph()
set_optimizer(g, HiGHS.Optimizer)

@optinode(g, nodes[1:5])

for (i, node) in enumerate(nodes)
    @variable(node, x[1:2] >= 0)
    @constraint(node, x[1] + 2 * x[2] >= 1)
    @objective(node, Min, x[1] + x[2])
end

for i in 1:4
    @linkconstraint(g, nodes[i][:x][1] + nodes[i + 1][:x][1] == i + 1)
    @linkconstraint(g, nodes[i][:x][2] + nodes[i + 1][:x][2] >= i * .3)
end

part_vector = [1, 2, 3, 4, 5]
partition = Partition(g, part_vector)

apply_partition!(g, partition)
set_optimizer(g, HiGHS.Optimizer)
set_to_node_objectives(g)

subgraphs = local_subgraphs(g)
for i in 1:length(subgraphs)
    set_optimizer(subgraphs[i], optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false))
end

BendersOpt = BendersAlgorithm(g, subgraphs[1]; add_slacks = true, slack_penalty = 100000)

t2 = @elapsed run_algorithm!(BendersOpt)

xsDDP = zeros(5, 2)

lbs = DDPOpt.lower_bounds
ubs = DDPOpt.upper_bounds
best_upper_bound = [minimum(ubs[1:i]) for i in 1:BendersOpt.current_iter]

plot(1:BendersOpt.current_iter, lbs)
plot!(1:BendersOpt.current_iter, best_upper_bound)

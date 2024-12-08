using Revise
using Plasmo, JuMP, HiGHS, Plots, PlasmoBenders

g = OptiGraph()
set_optimizer(g, HiGHS.Optimizer)

@optinode(g, nodes[1:5])

for (i, node) in enumerate(nodes)
    @variable(node, x[1:2] >= 0, Bin)
    @variable(node, y[1:2] >= 0)
    @constraint(node, x[1] + 2 * x[2] >= 1)
    @constraint(node, y[1] + y[2] * 2 >= 2)
    @objective(node, Min, x[1] + x[2] + y[1] + y[2])
end

for i in 1:4
    @linkconstraint(g, nodes[i][:x][1] + nodes[i + 1][:x][1] >= 1)
    @linkconstraint(g, nodes[i][:x][2] + nodes[i + 1][:x][2] >= .5)
    @linkconstraint(g, nodes[i][:y][1] + nodes[i + 1][:y][2] >= i * .5)
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

BendersOpt = BendersAlgorithm(g, subgraphs[1]; max_iters = 30, add_slacks = true, fix_slacks = true, strengthened = true)

t2 = @elapsed begin
    run_algorithm!(BendersOpt)
end

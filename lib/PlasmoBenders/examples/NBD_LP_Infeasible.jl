using Revise
using Plasmo, JuMP, HiGHS, Plots, PlasmoBenders
a=1

g = OptiGraph()
set_optimizer(g, HiGHS.Optimizer)

@optinode(g, nodes[1:2])

@variable(nodes[1], x >= 0)
@variable(nodes[2], 1 >= y >= 0)
@objective(nodes[1], Min, x)
@objective(nodes[2], Min, 2 * y)
@linkconstraint(g, x + y >= 2)

set_optimizer(nodes[1], optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false))
set_optimizer(nodes[2], optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false))

solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)

BendersOpt = BendersAlgorithm(g, nodes[1]; max_iters = 30, feasibility_cuts = true, solver = solver)

t2 = @elapsed begin
    run_algorithm!(BendersOpt)
end

using Revise
using Plasmo, JuMP, HiGHS, Plots, PlasmoBenders
using Plasmo, HiGHS, JuMP

g0 = OptiGraph()
g1 = OptiGraph()
g21 = OptiGraph()
g22 = OptiGraph()

@optinode(g1, n)
@variable(n, x[1:2] >= 0)
@variable(n, y, Bin)
@constraint(n, 2 * x[1] + x[2] + y == 3)
@objective(n, Min, x[1] + 2 * x[2] + y)

@optinode(g21, n)
@variable(n, x >= 0)
@objective(n, Min, 4 * x)

@optinode(g22, n)
@variable(n, x >= 1)
@objective(n, Min, x)

add_subgraph!(g0, g1)
add_subgraph!(g0, g21)
add_subgraph!(g0, g22)

set_to_node_objectives(g1)
set_to_node_objectives(g21)
set_to_node_objectives(g22)

@linkconstraint(g0, g1[:n][:x][2] + g21[:n][:x] + y >= 3)
@linkconstraint(g0, g1[:n][:x][1] + g22[:n][:x] >= 1)

gs = [g0, g1, g21, g22]
solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
for i in gs
    set_optimizer(i, solver)
end

BendersOpt = BendersAlgorithm(g0, g1; max_iters = 20, parallelize_benders = true, regularize = false)

t2 = @elapsed begin
    run_algorithm!(BendersOpt)
end

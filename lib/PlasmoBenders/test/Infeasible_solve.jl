module Test_Infeasible
using Plasmo, HiGHS, JuMP, PlasmoBenders, Test

function build_graph(MILP = false)
    g = OptiGraph()
    set_optimizer(g, HiGHS.Optimizer)

    @optinode(g, nodes[1:2])

    @variable(nodes[1], x >= 0)

    if MILP
        @variable(nodes[2], y, Bin)
    else
        @variable(nodes[2], 1 >= y >= 0)
    end

    @objective(nodes[1], Min, x)
    @objective(nodes[2], Min, 2 * y)
    @linkconstraint(g, x + y >= 2)

    return g
end

function get_gap(DDPOpt::BendersAlgorithm)
    ub = DDPOpt.best_upper_bound
    lb = DDPOpt.lower_bounds[end]

    gap = (ub - lb) / (abs(lb))

    return gap
end

function test_graphs()
    solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)

    gtest = build_graph()
    DDPOpt = BendersAlgorithm(gtest, gtest[:nodes][1], max_iters = 20, feasibility_cuts = true, solver = solver);
    run_algorithm!(DDPOpt)
    @test isapprox(DDPOpt.best_upper_bound, 2.0, rtol = 1e-6)
    @test DDPOpt.status == MOI.OPTIMAL

    gtest = build_graph(true)
    DDPOpt = BendersAlgorithm(gtest, gtest[:nodes][1], max_iters = 20, feasibility_cuts = true, solver = solver);
    run_algorithm!(DDPOpt)
    @test isapprox(DDPOpt.best_upper_bound, 2.0, rtol = 1e-6)
    @test DDPOpt.status == MOI.OPTIMAL
end

function run_tests()
    for name in names(@__MODULE__; all=true)
        if !startswith("$(name)", "test_")
            continue
        end
        @testset "$(name)" begin
            getfield(@__MODULE__, name)()
        end
    end
end

end


Test_Infeasible.run_tests()

using Plasmo, HiGHS, JuMP, PlasmoBenders, Test, Distributed

if nprocs() < 2
    addprocs(1)
end
@everywhere begin
    using Plasmo, HiGHS, JuMP, PlasmoBenders, Test, Distributed
end

module Test_Infeasible
using Plasmo, HiGHS, JuMP, PlasmoBenders, Test, Distributed

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

function build_remote_graph(MILP = false)
    rg = RemoteOptiGraph(worker = 2)
    rg1 = RemoteOptiGraph(worker = 2)
    rg2 = RemoteOptiGraph(worker = 2)

    @optinode(rg1, nodes1)
    @optinode(rg2, nodes2)

    @variable(nodes1, x >= 0)

    if MILP
        @variable(nodes2, y, Bin)
    else
        @variable(nodes2, 0 <= y <= 1)
    end

    @objective(rg1, Min, x)
    @objective(rg2, Min, 2 * y)
    add_subgraph(rg, rg1)
    add_subgraph(rg, rg2)
    @linkconstraint(rg, x + y >= 2)

    return rg
end

function get_gap(BendersAlg::BendersAlgorithm)
    ub = BendersAlg.best_upper_bound
    lb = BendersAlg.lower_bounds[end]

    gap = (ub - lb) / (abs(lb))

    return gap
end

function test_graphs()
    solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)

    gtest = build_graph()
    BendersAlg = BendersAlgorithm(gtest, gtest[:nodes][1], max_iters = 20, feasibility_cuts = true, solver = solver);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 2.0, rtol = 1e-6)
    @test BendersAlg.status == MOI.OPTIMAL

    gtest = build_graph(true)
    BendersAlg = BendersAlgorithm(gtest, gtest[:nodes][1], max_iters = 20, feasibility_cuts = true, solver = solver);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 2.0, rtol = 1e-6)
    @test BendersAlg.status == MOI.OPTIMAL
    
    gtest = build_remote_graph(false)
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, feasibility_cuts = true, solver = solver);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 2.0, rtol = 1e-6)
    @test BendersAlg.status == MOI.OPTIMAL

    gtest = build_remote_graph(true)
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, feasibility_cuts = true, solver = solver);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 2.0, rtol = 1e-6)
    @test BendersAlg.status == MOI.OPTIMAL
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

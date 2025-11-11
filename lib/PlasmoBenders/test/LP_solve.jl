using Plasmo, HiGHS, JuMP, PlasmoBenders, Test, Distributed

if nprocs() < 2
    addprocs(1)
end
@everywhere begin
    using Plasmo, HiGHS, JuMP, PlasmoBenders, Test, Distributed
end

module Test_LP_solves

using Plasmo, HiGHS, JuMP, PlasmoBenders, Test

function build_graph(single_objective = false)
    solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
    graph = OptiGraph()
    set_optimizer(graph, solver)
    @optinode(graph, nodes[1:3])

    for (i, node) in enumerate(nodes)
        @variable(node, 0 <= x <= 1)
        @variable(node, y >= 0)
        @constraint(node, x + y >= 1.3)

        if single_objective
            @objective(node, Min, x)
        else
            @objective(node, Min, x + 2 * y)
        end
    end

    for i in 1:2
        @linkconstraint(graph, nodes[i][:x] + nodes[i + 1][:y] >= i)
    end

    part_vector = [1, 2, 3]
    partition = Partition(graph, part_vector)

    apply_partition!(graph, partition)
    set_optimizer(graph, HiGHS.Optimizer)
    set_to_node_objectives(graph)

    subgraphs = local_subgraphs(graph)
    for i in 1:length(subgraphs)
        set_to_node_objectives(subgraphs[i])
        set_optimizer(subgraphs[i], solver)
    end

    return graph
end

function build_remote_graph(single_objective = false)
    solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
    rgraph = RemoteOptiGraph(worker=2)
    set_optimizer(rgraph, solver)
    rg1 = RemoteOptiGraph(worker=2)
    rg2 = RemoteOptiGraph(worker=2)
    rg3 = RemoteOptiGraph(worker=2)
    @optinode(rg1, nodes1)
    @optinode(rg2, nodes2)
    @optinode(rg3, nodes3)
    add_subgraph(rgraph, rg1)
    add_subgraph(rgraph, rg2)
    add_subgraph(rgraph, rg3)

    nodes = [nodes1, nodes2, nodes3]

    for (i, node) in enumerate(nodes)
        @variable(node, 0 <= x <= 1)
        @variable(node, y >= 0)
        @constraint(node, x + y >= 1.3)

        if single_objective
            @objective(node, Min, x)
        else
            @objective(node, Min, x + 2 * y)
        end
    end

    for i in 1:2
        @linkconstraint(rgraph, nodes[i][:x] + nodes[i + 1][:y] >= i)
    end

    subgraphs = local_subgraphs(rgraph)
    for i in 1:length(subgraphs)
        set_to_node_objectives(subgraphs[i])
        set_optimizer(subgraphs[i], solver)
    end

    return rgraph
end

function get_gap(BendersAlg::BendersAlgorithm)
    ub = BendersAlg.best_upper_bound
    lb = BendersAlg.lower_bounds[end]

    gap = (ub - lb) / (abs(lb))

    return gap
end

function test_graphs()
    # OptiGraphs
    gtest = build_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test BendersAlg.status == MOI.OPTIMAL

    gtest = build_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 2);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test BendersAlg.status == MOI.ITERATION_LIMIT

    gtest = build_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, time_limit = 0.00001);
    run_algorithm!(BendersAlg)
    @test BendersAlg.status == MOI.TIME_LIMIT

    gtest = build_graph(true)
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, add_slacks = true);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 0, rtol = 1e-6)
    @test BendersAlg.status == MOI.OPTIMAL

    gtest = build_graph(true)
    DDPOpt = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, add_slacks = false);
    run_algorithm!(DDPOpt)
    @test isapprox(DDPOpt.best_upper_bound, 0, rtol = 1e-6)

    gtest = build_graph(true)
    DDPOpt = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, add_slacks = false, sequential_backward_pass = true);
    run_algorithm!(DDPOpt)
    @test isapprox(DDPOpt.best_upper_bound, 0, rtol = 1e-6)

    gtest = build_graph()
    DDPOpt = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, multicut = true);
    run_algorithm!(DDPOpt)
    @test isapprox(DDPOpt.best_upper_bound, 5.5, rtol = 1e-6)
    @test isapprox(get_gap(DDPOpt), 0, rtol = 1e-6)

    gtest = build_graph()
    DDPOpt = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, multicut = true, sequential_backward_pass = true);
    run_algorithm!(DDPOpt)
    @test isapprox(DDPOpt.best_upper_bound, 5.5, rtol = 1e-6)
    @test isapprox(get_gap(DDPOpt), 0, rtol = 1e-6)

    gtest = build_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, multicut = false);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test isapprox(get_gap(BendersAlg), 0, rtol = 1e-6)

    gtest = build_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, multicut = true, regularize = false, parallelize_backward = true);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test isapprox(get_gap(BendersAlg), 0, rtol = 1e-6)

    gtest = build_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, multicut = true, regularize = false, parallelize_benders = false);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test isapprox(get_gap(BendersAlg), 0, rtol = 1e-6)

    gtest = build_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[2], max_iters = 20, multicut = true, parallelize_benders = true);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test isapprox(get_gap(BendersAlg), 0, rtol = 1e-6)

    # RemoteOptiGraphs
    gtest = build_remote_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test BendersAlg.status == MOI.OPTIMAL

    gtest = build_remote_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 2);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test BendersAlg.status == MOI.ITERATION_LIMIT

    gtest = build_remote_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, time_limit = 0.00001);
    run_algorithm!(BendersAlg)
    @test BendersAlg.status == MOI.TIME_LIMIT

    gtest = build_remote_graph(true)
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, add_slacks = true);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 0, rtol = 1e-6)
    @test BendersAlg.status == MOI.OPTIMAL

    gtest = build_remote_graph(true)
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, add_slacks = false);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 0, rtol = 1e-6)

    gtest = build_remote_graph(true)
    DDPOpt = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, add_slacks = false, sequential_backward_pass = true);
    run_algorithm!(DDPOpt)
    @test isapprox(DDPOpt.best_upper_bound, 0, rtol = 1e-6)

    gtest = build_remote_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, multicut = true);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test isapprox(get_gap(BendersAlg), 0, rtol = 1e-6)

    gtest = build_remote_graph()
    DDPOpt = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, multicut = true, sequential_backward_pass = true);
    run_algorithm!(DDPOpt)
    @test isapprox(DDPOpt.best_upper_bound, 5.5, rtol = 1e-6)
    @test isapprox(get_gap(DDPOpt), 0, rtol = 1e-6)

    gtest = build_remote_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, multicut = false);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test isapprox(get_gap(BendersAlg), 0, rtol = 1e-6)

    gtest = build_remote_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, multicut = true, regularize = false, parallelize_backward = true);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test isapprox(get_gap(BendersAlg), 0, rtol = 1e-6)

    gtest = build_remote_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[1], max_iters = 20, multicut = true, regularize = false, parallelize_benders = false);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test isapprox(get_gap(BendersAlg), 0, rtol = 1e-6)

    gtest = build_remote_graph()
    BendersAlg = BendersAlgorithm(gtest, local_subgraphs(gtest)[2], max_iters = 20, multicut = true, parallelize_benders = true);
    run_algorithm!(BendersAlg)
    @test isapprox(BendersAlg.best_upper_bound, 5.5, rtol = 1e-6)
    @test isapprox(get_gap(BendersAlg), 0, rtol = 1e-6)
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

Test_LP_solves.run_tests()

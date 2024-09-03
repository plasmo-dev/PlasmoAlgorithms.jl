module TestDDP_MIP_solves

    using Plasmo, HiGHS, JuMP, PlasmoBenders, Test

    function build_graph()
        solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
        graph = OptiGraph()
        set_optimizer(graph, solver)
        @optinode(graph, nodes[1:3])

        for (i, node) in enumerate(nodes)
            set_optimizer(graph, solver)
            @variable(node, x, Bin)
            @variable(node, y >= 0)
            @constraint(node, x + y >= 1.3)

            @objective(node, Min, x + 2 * y)
        end

        for i in 1:2
            @linkconstraint(graph, nodes[i][:x] + nodes[i + 1][:y] >= i)
        end

        part_vector = [1, 2, 3]
        partition = Partition(graph, part_vector)

        apply_partition!(graph, partition)
        set_optimizer(graph, HiGHS.Optimizer)
        set_to_node_objectives(graph)

        subgraphs = getsubgraphs(graph)
        for i in 1:length(subgraphs)
            set_to_node_objectives(subgraphs[i])
            set_optimizer(subgraphs[i], solver)
        end

        return graph
    end

    function get_gap(DDPOpt::BendersOptimizer)
        ub = DDPOpt.best_upper_bound
        lb = DDPOpt.lower_bounds[end]

        gap = (ub - lb) / (abs(lb))

        return gap
    end

    gtest = build_graph()
    DDPOpt = BendersOptimizer(gtest, getsubgraphs(gtest)[1], max_iters = 20);
    optimize!(DDPOpt)
    @test isapprox(DDPOpt.best_upper_bound, 5.8, rtol = 1e-6)

    gtest = build_graph()
    DDPOpt = BendersOptimizer(gtest, getsubgraphs(gtest)[1], max_iters = 20, strengthened = true);
    optimize!(DDPOpt)
    @test isapprox(DDPOpt.best_upper_bound, 5.8, rtol = 1e-6)
    @test isapprox(get_gap(DDPOpt), 0, rtol = 1e-6)

    gtest = build_graph()
    DDPOpt = BendersOptimizer(gtest, getsubgraphs(gtest)[1], max_iters = 20, strengthened = true, multicut = true);
    optimize!(DDPOpt)
    @test isapprox(DDPOpt.best_upper_bound, 5.8, rtol = 1e-6)
    @test isapprox(get_gap(DDPOpt), 0, rtol = 1e-6)

    gtest = build_graph()
    DDPOpt = BendersOptimizer(gtest, getsubgraphs(gtest)[1], max_iters = 20, strengthened = true, multicut = true, regularize = false, parallelize_backward = true
    );
    optimize!(DDPOpt)
    @test isapprox(DDPOpt.best_upper_bound, 5.8, rtol = 1e-6)
    @test isapprox(get_gap(DDPOpt), 0, rtol = 1e-6)

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

TestDDP_MIP_solves.run_tests()

module Test_Tree

using PlasmoBenders, Plasmo, JuMP, HiGHS, Test

function build_graph()
    solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
    graph = OptiGraph()
    set_optimizer(graph, solver)
    @optinode(graph, nodes[1:4])

    for (i, node) in enumerate(nodes)
        set_optimizer(graph, solver)
        @variable(node, x, Bin)
        @variable(node, y >= 0)
        @constraint(node, x + y >= 1.3)

        @objective(node, Min, x + 2 * y)
    end

    @linkconstraint(graph, nodes[1][:x] + nodes[2][:y] >= 1)
    @linkconstraint(graph, nodes[1][:x] + nodes[3][:y] >= 2)
    @linkconstraint(graph, nodes[3][:x] + nodes[4][:y] >= 3)

    part_vector = [1, 2, 3, 4]
    partition = Partition(graph, part_vector)

    apply_partition!(graph, partition)
    set_optimizer(graph, HiGHS.Optimizer)
    set_to_node_objectives(graph)

    subgraphs = local_subgraphs(graph)
    for i in 1:length(subgraphs)
        set_to_node_objectives(subgraphs[i])
        set_optimizer(subgraphs[i], optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false))
    end

    return graph, nodes, subgraphs
end

function test_DDP_Tree_initialization()

    graph, nodes, subgraphs = build_graph()

    max_iters = 50
    tol = 1e-5
    M = -100
    DDPOpt = BendersAlgorithm(graph, subgraphs[1]; max_iters = max_iters, tol = tol, M = M)

    @test length(DDPOpt.solve_order) == 4
    @test length(keys(DDPOpt.solve_order_dict)) == 4

    @test length(DDPOpt.solve_order_dict[subgraphs[1]]) == 2
    @test length(DDPOpt.solve_order_dict[subgraphs[3]]) == 1

    @test DDPOpt.parent_objects[subgraphs[2]] == subgraphs[1]
    @test DDPOpt.parent_objects[subgraphs[3]] == subgraphs[1]
    @test DDPOpt.parent_objects[subgraphs[4]] == subgraphs[3]

    @test get_parallelize_benders(DDPOpt) == false
    @test get_parallelize_forward(DDPOpt) == false
    @test get_parallelize_backward(DDPOpt) == false

    graph, nodes, subgraphs = build_graph()
    @test_throws ErrorException BendersAlgorithm(graph, subgraphs[1], parallelize_benders = true)

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

Test_Tree.run_tests()

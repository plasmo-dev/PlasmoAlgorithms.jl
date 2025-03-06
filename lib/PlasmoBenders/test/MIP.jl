module Test_MIP

using PlasmoBenders, Plasmo, JuMP, HiGHS, Test

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

subgraphs = local_subgraphs(graph)
for i in 1:length(subgraphs)
    set_to_node_objectives(subgraphs[i])
    set_optimizer(subgraphs[i], optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false))
end

graphtest = OptiGraph()
@optinode(graphtest, nodetest)


function test_DDP_MIP_initialization()
    @test_throws ErrorException BendersAlgorithm(graph, graphtest)

    max_iters = 50
    tol = 1e-5
    M = -100
    DDPOpt = BendersAlgorithm(graph, subgraphs[1]; max_iters = max_iters, tol = tol, M = M)

    @test DDPOpt.graph == graph
    @test DDPOpt.root_object == subgraphs[1]
    @test DDPOpt.is_MIP == true
    @test DDPOpt.max_iters == max_iters
    @test DDPOpt.tol == tol
    @test DDPOpt.M == M
    @test DDPOpt.solve_order == OptiGraph[i for i in subgraphs]

    @test "_theta_node[:_theta[1]]" in name.(all_variables(subgraphs[1]))
    @test "_theta_node[:_theta[1]]" in name.(all_variables(subgraphs[2]))

    @test JuMP.lower_bound(subgraphs[1][:_theta_node][:_theta][1]) == M
    @test JuMP.lower_bound(subgraphs[2][:_theta_node][:_theta][1]) == M

    obj_func1 = objective_function(subgraphs[1])
    obj_func2 = objective_function(subgraphs[2])
    @test subgraphs[1][:_theta_node][:_theta][1] in keys(obj_func1.terms)
    @test subgraphs[2][:_theta_node][:_theta][1] in keys(obj_func2.terms)

    @test length(DDPOpt.comp_vars) == 2
    @test DDPOpt.comp_vars[subgraphs[2]] == NodeVariableRef[nodes[1][:x]]
    @test DDPOpt.comp_vars[subgraphs[3]] == NodeVariableRef[nodes[2][:x]]

    var_copy_map1 = DDPOpt.var_copy_map[subgraphs[2]]
    var_copy_map2 = DDPOpt.var_copy_map[subgraphs[3]]
    @test var_copy_map1[nodes[1][:x]] == nodes[2][:_comp_vars_copy][1]
    @test var_copy_map2[nodes[2][:x]] == nodes[3][:_comp_vars_copy][1]

    for (i, node) in enumerate(nodes)
        @test DDPOpt.binary_map[subgraphs[i]][node[:x]] == 0
    end

    # Test API
    @test get_graph(DDPOpt) == graph
    @test get_root_object(DDPOpt) == subgraphs[1]
    @test get_max_iters(DDPOpt) == max_iters
    @test get_tol(DDPOpt) == tol
    @test get_time_forward_pass(DDPOpt) == 0
    @test get_time_backward_pass(DDPOpt) == 0
    @test get_time_init(DDPOpt) > 0
    @test typeof(get_time_iterations(DDPOpt)) <: Vector
    @test typeof(get_lower_bounds(DDPOpt)) <: Vector
    @test typeof(get_upper_bounds(DDPOpt, monotonic = true)) <: Vector
    @test typeof(get_upper_bounds(DDPOpt, monotonic = false)) <: Vector

    # Test options
    set_strengthened!(DDPOpt, true)
    @test get_strengthened(DDPOpt) == true
    set_multicut!(DDPOpt, false)
    @test get_multicut(DDPOpt) == false
    set_regularize!(DDPOpt, true)
    @test get_regularize(DDPOpt) == true
    set_parallelize_benders!(DDPOpt, true)
    @test get_parallelize_benders(DDPOpt) == true
    set_parallelize_forward!(DDPOpt, true)
    @test get_parallelize_forward(DDPOpt) == true
    set_parallelize_backward!(DDPOpt, true)
    @test get_parallelize_backward(DDPOpt) == true
    set_add_slacks!(DDPOpt, true)
    @test get_add_slacks(DDPOpt) == true
    set_fix_slacks!(DDPOpt, true)
    @test get_fix_slacks(DDPOpt) == true
    set_warm_start!(DDPOpt, false)
    @test get_warm_start(DDPOpt) == false
    set_relaxed_init_cuts!(DDPOpt, true)
    @test get_relaxed_init_cuts(DDPOpt) == true
    set_slack_penalty!(DDPOpt, 1e5)
    @test isapprox(get_slack_penalty(DDPOpt), 1e5, rtol = 1e-6)
    set_regularize_param!(DDPOpt, 0.1)
    @test isapprox(get_regularize_param(DDPOpt), 0.1, rtol = 1e-6)
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

Test_MIP.run_tests()

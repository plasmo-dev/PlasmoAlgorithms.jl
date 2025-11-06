using Plasmo, HiGHS, JuMP, PlasmoBenders, Test, Distributed

if nprocs() < 2
    addprocs(1)
end
@everywhere begin
    using Plasmo, HiGHS, JuMP, PlasmoBenders, Test, Distributed
end

module Test_Initialization

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


solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
rgraph = RemoteOptiGraph(worker=2)
set_optimizer(graph, solver)
rg1 = RemoteOptiGraph(worker=2)
rg2 = RemoteOptiGraph(worker=2)
rg3 = RemoteOptiGraph(worker=2)
@optinode(rg1, rnodes1)
@optinode(rg2, rnodes2)
@optinode(rg3, rnodes3)
add_subgraph(rgraph, rg1)
add_subgraph(rgraph, rg2)
add_subgraph(rgraph, rg3)

rnodes = [rnodes1, rnodes2, rnodes3]

for (i, node) in enumerate(rnodes)
    set_optimizer(graph, solver)
    @variable(node, x, Bin)
    @variable(node, y >= 0)
    @constraint(node, x + y >= 1.3)

    @objective(node, Min, x + 2 * y)
end

for i in 1:2
    @linkconstraint(rgraph, rnodes[i][:x] + rnodes[i + 1][:y] >= i)
end

rsubgraphs = local_subgraphs(rgraph)
for i in 1:length(rsubgraphs)
    set_to_node_objectives(rsubgraphs[i])
    set_optimizer(rsubgraphs[i], optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false))
end

rgraphtest = RemoteOptiGraph(worker = 2)
@optinode(rgraphtest, rnodetest)


function test_DDP_MIP_initialization()
    @test_throws ErrorException BendersAlgorithm(graph, graphtest)

    max_iters = 50
    tol = 1e-5
    M = -100
    BendersAlg = BendersAlgorithm(graph, subgraphs[1]; max_iters = max_iters, tol = tol, M = M)

    @test BendersAlg.graph == graph
    @test BendersAlg.root_object == subgraphs[1]
    @test BendersAlg.is_MIP == true
    @test BendersAlg.max_iters == max_iters
    @test BendersAlg.tol == tol
    @test BendersAlg.M == M
    @test BendersAlg.solve_order == OptiGraph[i for i in subgraphs]

    @test "_theta_node[:_theta[1]]" in name.(all_variables(subgraphs[1]))
    @test "_theta_node[:_theta[1]]" in name.(all_variables(subgraphs[2]))

    @test JuMP.lower_bound(subgraphs[1][:_theta_node][:_theta][1]) == M
    @test JuMP.lower_bound(subgraphs[2][:_theta_node][:_theta][1]) == M

    obj_func1 = objective_function(subgraphs[1])
    obj_func2 = objective_function(subgraphs[2])
    @test subgraphs[1][:_theta_node][:_theta][1] in keys(obj_func1.terms)
    @test subgraphs[2][:_theta_node][:_theta][1] in keys(obj_func2.terms)

    @test length(BendersAlg.comp_vars) == 2
    @test BendersAlg.comp_vars[subgraphs[2]] == NodeVariableRef[nodes[1][:x]]
    @test BendersAlg.comp_vars[subgraphs[3]] == NodeVariableRef[nodes[2][:x]]

    var_copy_map1 = BendersAlg.var_copy_map[subgraphs[2]]
    var_copy_map2 = BendersAlg.var_copy_map[subgraphs[3]]
    @test var_copy_map1[nodes[1][:x]] == nodes[2][:_comp_vars_copy][1]
    @test var_copy_map2[nodes[2][:x]] == nodes[3][:_comp_vars_copy][1]

    for (i, node) in enumerate(nodes)
        @test BendersAlg.binary_map[subgraphs[i]][node[:x]] == 0
    end

    # Test API
    @test PlasmoBenders.get_graph(BendersAlg) == graph
    @test get_root_object(BendersAlg) == subgraphs[1]
    @test get_max_iters(BendersAlg) == max_iters
    @test get_tol(BendersAlg) == tol
    @test get_time_forward_pass(BendersAlg) == 0
    @test get_time_backward_pass(BendersAlg) == 0
    @test get_time_init(BendersAlg) > 0
    @test get_time_root_problem_solve(BendersAlg) == 0
    @test get_time_subproblem_solves(BendersAlg) == 0
    @test typeof(get_time_iterations(BendersAlg)) <: Vector
    @test typeof(get_lower_bounds(BendersAlg)) <: Vector
    @test typeof(get_upper_bounds(BendersAlg, monotonic = true)) <: Vector
    @test typeof(get_upper_bounds(BendersAlg, monotonic = false)) <: Vector

    # Test options
    set_strengthened!(BendersAlg, true)
    @test get_strengthened(BendersAlg) == true
    set_multicut!(BendersAlg, false)
    @test get_multicut(BendersAlg) == false
    set_regularize!(BendersAlg, true)
    @test get_regularize(BendersAlg) == true
    set_parallelize_benders!(BendersAlg, true)
    @test get_parallelize_benders(BendersAlg) == true
    set_parallelize_forward!(BendersAlg, true)
    @test get_parallelize_forward(BendersAlg) == true
    set_parallelize_backward!(BendersAlg, true)
    @test get_parallelize_backward(BendersAlg) == true
    set_add_slacks!(BendersAlg, true)
    @test get_add_slacks(BendersAlg) == true
    set_warm_start!(BendersAlg, false)
    @test get_warm_start(BendersAlg) == false
    set_relaxed_init_cuts!(BendersAlg, true)
    @test get_relaxed_init_cuts(BendersAlg) == true
    set_slack_penalty!(BendersAlg, 1e5)
    @test isapprox(get_slack_penalty(BendersAlg), 1e5, rtol = 1e-6)
    set_regularize_param!(BendersAlg, 0.1)
    @test isapprox(get_regularize_param(BendersAlg), 0.1, rtol = 1e-6)

    set_objective_sense(subgraphs[1], MOI.MAX_SENSE)
    @test_throws ErrorException BendersAlgorithm(graph, subgraphs[1])
end

function test_DDP_MIP_initialization_remote()
    @test_throws ErrorException BendersAlgorithm(rgraph, rgraphtest)

    max_iters = 50
    tol = 1e-5
    M = -100
    BendersAlg = BendersAlgorithm(rgraph, rsubgraphs[1]; max_iters = max_iters, tol = tol, M = M)

    @test BendersAlg.graph == rgraph
    @test BendersAlg.root_object == rsubgraphs[1]
    @test BendersAlg.is_MIP == true
    @test BendersAlg.max_iters == max_iters
    @test BendersAlg.tol == tol
    @test BendersAlg.M == M
    @test BendersAlg.solve_order == rsubgraphs

    @test JuMP.lower_bound(rsubgraphs[1][:_theta_node][:_theta][1]) == M
    @test JuMP.lower_bound(rsubgraphs[2][:_theta_node][:_theta][1]) == M

    obj_func1 = objective_function(rsubgraphs[1])
    obj_func2 = objective_function(rsubgraphs[2])
    @test rsubgraphs[1][:_theta_node][:_theta][1] in keys(obj_func1.terms)
    @test rsubgraphs[2][:_theta_node][:_theta][1] in keys(obj_func2.terms)

    @test length(BendersAlg.comp_vars) == 2
    @test BendersAlg.comp_vars[rsubgraphs[2]] == RemoteVariableRef[rnodes[1][:x]]
    @test BendersAlg.comp_vars[rsubgraphs[3]] == RemoteVariableRef[rnodes[2][:x]]

    var_copy_map1 = BendersAlg.var_copy_map[rsubgraphs[2]]
    var_copy_map2 = BendersAlg.var_copy_map[rsubgraphs[3]]
    @test var_copy_map1[rnodes[1][:x]] == rnodes[2][:_comp_vars_copy][1]
    @test var_copy_map2[rnodes[2][:x]] == rnodes[3][:_comp_vars_copy][1]

    for (i, node) in enumerate(rnodes)
        @test BendersAlg.binary_map[rsubgraphs[i]][node[:x]] == 0
    end

    # Test API
    @test PlasmoBenders.get_graph(BendersAlg) == rgraph
    @test get_root_object(BendersAlg) == rsubgraphs[1]
    @test get_max_iters(BendersAlg) == max_iters
    @test get_tol(BendersAlg) == tol
    @test get_time_forward_pass(BendersAlg) == 0
    @test get_time_backward_pass(BendersAlg) == 0
    @test get_time_init(BendersAlg) > 0
    @test get_time_forward_pass(BendersAlg) == 0
    @test get_time_backward_pass(BendersAlg) == 0
    @test typeof(get_time_iterations(BendersAlg)) <: Vector
    @test typeof(get_lower_bounds(BendersAlg)) <: Vector
    @test typeof(get_upper_bounds(BendersAlg, monotonic = true)) <: Vector
    @test typeof(get_upper_bounds(BendersAlg, monotonic = false)) <: Vector

    # Test options
    set_strengthened!(BendersAlg, true)
    @test get_strengthened(BendersAlg) == true
    set_multicut!(BendersAlg, false)
    @test get_multicut(BendersAlg) == false
    set_regularize!(BendersAlg, true)
    @test get_regularize(BendersAlg) == true
    set_parallelize_benders!(BendersAlg, true)
    @test get_parallelize_benders(BendersAlg) == true
    set_parallelize_forward!(BendersAlg, true)
    @test get_parallelize_forward(BendersAlg) == true
    set_parallelize_backward!(BendersAlg, true)
    @test get_parallelize_backward(BendersAlg) == true
    set_add_slacks!(BendersAlg, true)
    @test get_add_slacks(BendersAlg) == true
    set_warm_start!(BendersAlg, false)
    @test get_warm_start(BendersAlg) == false
    set_relaxed_init_cuts!(BendersAlg, true)
    @test get_relaxed_init_cuts(BendersAlg) == true
    set_slack_penalty!(BendersAlg, 1e5)
    @test isapprox(get_slack_penalty(BendersAlg), 1e5, rtol = 1e-6)
    set_regularize_param!(BendersAlg, 0.1)
    @test isapprox(get_regularize_param(BendersAlg), 0.1, rtol = 1e-6)
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

Test_Initialization.run_tests()


"""

"""
function expand_solve_order_objects!(
    optimizer::BendersAlgorithm{Plasmo.OptiGraph},
    expansion_set::Vector{Vector{Plasmo.OptiNode}}
)
    if length(expansion_set) != length(optimizer.solve_order)
        error("expansion_set is not the right length")
    end

    graph = optimizer.graph
    Plasmo._init_graph_backend(graph)
    hypergraph, hyper_map = Plasmo.graph_backend_data(graph)

    new_solve_order = Plasmo.OptiGraph[]

    optimizer.ext["extra_objective"] = Dict{Plasmo.OptiGraph, AffExpr}()

    #TODO: test that new node set does not contain any of the original nodes in the subgraph
    for (i, sub) in enumerate(optimizer.solve_order)
        nodes = all_nodes(sub)
        new_node_set = expansion_set[i]

        if length(new_node_set) != 0
            hypernodes = getindex.(Ref(hyper_map), nodes)
            hypernodes_to_add = getindex.(Ref(hyper_map), new_node_set)
            new_nodes = union(hypernodes, hypernodes_to_add)
            new_edges = induced_edges(hypergraph, new_nodes)

            new_optinodes = getindex.(Ref(hyper_map), new_nodes)
            new_optiedges = getindex.(Ref(hyper_map), new_edges)

            expanded_subgraph = OptiGraph(new_optinodes, new_optiedges)
            push!(new_solve_order, expanded_subgraph)
            theta_var = sub[:_theta_node][:_theta]
            optimizer.ext["theta_vars"][expanded_subgraph] = theta_var
            new_sub = expanded_subgraph
        else
            push!(new_solve_order, sub)
            new_sub = sub
        end

        extra_objective = AffExpr()
        for n in new_node_set
            add_to_expression!(extra_objective, objective_function(n))
        end

        optimizer.ext["extra_objective"][new_sub] = extra_objective

        if i != length(optimizer.solve_order) # TODO: make this storage more efficient
            optimizer.dual_iters[new_sub] = optimizer.dual_iters[sub]
            optimizer.primal_iters[new_sub] = optimizer.primal_iters[sub]
            optimizer.phis[new_sub] = optimizer.phis[sub]
            optimizer.phis_LR[new_sub] = optimizer.phis_LR[sub]
            optimizer.comp_vars[new_sub] = optimizer.comp_vars[sub]
            optimizer.comp_var_map[new_sub] = optimizer.comp_var_map[sub]
            optimizer.var_copy_map[new_sub] = optimizer.var_copy_map[sub]

        end

        # update these:
        optimizer.binary_map[new_sub] = optimizer.binary_map[sub]
        optimizer.integer_map[new_sub] = optimizer.integer_map[sub]
        optimizer.last_solutions[new_sub] = optimizer.last_solutions[sub]
    end

    optimizer.ext["is_overlapped"] = true

    optimizer.solve_order = new_solve_order

    for (i, object) in enumerate(optimizer.solve_order)

        object_var_solution_map = Dict{VariableRef, Int}()
        all_vars = all_variables(object)

        bin_value_dict = Dict{VariableRef, Float64}()
        int_value_dict = Dict{VariableRef, Float64}()
        for (j, var) in enumerate(all_vars)
            if JuMP.is_binary(var)
                bin_value_dict[var] = 0
            elseif JuMP.is_integer(var)
                int_value_dict[var] = 0
            end
            object_var_solution_map[var] = j
        end

        optimizer.var_solution_map[object] = object_var_solution_map
        optimizer.binary_map[object] = bin_value_dict
        optimizer.integer_map[object] = int_value_dict

        optimizer.last_solutions[object] = zeros(Float64, length(all_vars))
    end
end

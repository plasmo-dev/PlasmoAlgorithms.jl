function _init_ext!(optimizer::BendersAlgorithm{Plasmo.OptiGraph})
    graph = optimizer.graph
    subgraphs = local_subgraphs(graph)

    # Define mappings for variables and nodes to their subgraph
    var_to_graph_map = Dict{NodeVariableRef, Plasmo.OptiGraph}()
    node_to_graph_map = Dict{Plasmo.OptiNode, Plasmo.OptiGraph}()

    # Build mappings
    for (i, g) in enumerate(subgraphs)
        vars = JuMP.all_variables(g)
        for (j, var) in enumerate(vars)
            var_to_graph_map[var] = g
        end

        nodes = all_nodes(g)
        for (j, node) in enumerate(nodes)
            node_to_graph_map[node] = g
        end
    end

    # Save mappings to the optimizer object
    optimizer.ext["var_to_graph"] = var_to_graph_map #TODO: Make a new structure
    optimizer.ext["node_to_graph"] = node_to_graph_map
    optimizer.ext["theta_vars"] = Dict{Plasmo.OptiGraph, Vector{NodeVariableRef}}()
    optimizer.ext["is_overlapped"] = false
    optimizer.ext["incident_edges"] = Dict{Plasmo.OptiGraph, Vector{Plasmo.OptiEdge}}()

    return nothing
end

function _add_second_object!(optimizer::BendersAlgorithm{Plasmo.OptiGraph}, relaxed)
    graph = optimizer.graph
    root_object = optimizer.root_object
    parent_objects = optimizer.parent_objects
    solve_order_dict = optimizer.solve_order_dict

    var_to_graph_map = optimizer.ext["var_to_graph"]
    node_to_graph_map = optimizer.ext["node_to_graph"]

    optimizer.ext["search_next"] = Vector{typeof(root_object)}()

    # Get all nodes in the start object
    start_nodes = all_nodes(root_object)

    # Get incident edges
    edges = Plasmo.incident_edges(optimizer, graph, root_object)

    # Loop through the nodes in the incident edges, and add their owning subgraphs to `graphs`
    graphs = []
    for (i, edge) in enumerate(edges)
        nodes = edge.nodes
        for (j, node) in enumerate(nodes)
            node_graph = node_to_graph_map[node]
            if node_graph != root_object
                if !(node_graph in graphs)
                    push!(graphs, node_graph)
                end
            end
        end
    end

    # Error if the subgraphs are not linear
    if length(graphs) == 0
        error("The starting subgraph has no connections to other subgraphs in the graph")
    end

    solve_order_dict[root_object] = graphs

    # Add cost-to-go variable and add to objective function
    _add_cost_to_go!(optimizer, root_object, relaxed)
    #Plasmo._init_graph_backend(optimizer.graph)

    for g in graphs
        parent_objects[g] = root_object
        push!(optimizer.solve_order, g)
        push!(optimizer.ext["search_next"], g)
    end
end

function _add_cost_to_go!(
    optimizer::BendersAlgorithm{Plasmo.OptiGraph},
    last_object::Plasmo.OptiGraph,
    relaxed::Bool
)
    num_thetas = length(optimizer.solve_order_dict[last_object])

    # Define a new node on the subgraph and add theta to that node
    @optinode(last_object, _theta_node)
    @variable(last_object[:_theta_node], _theta[1:num_thetas] >= 0)

    # Add theta to the node's objective
    theta_sum = sum(_theta[i] for i in 1:num_thetas)
    @objective(last_object[:_theta_node], Min, theta_sum)

    # If the initial relaxation is not being solved, set lower bound as M
    if !relaxed || !(optimizer.is_MIP)
        for i in 1:num_thetas
            JuMP.set_lower_bound(last_object[:_theta_node][:_theta][i], optimizer.M)
        end
    end

    # Reset graph backend
    #Plasmo.set_graph_backend(optimizer.graph)
    #Plasmo.set_graph_backend(last_object)

    optimizer.ext["theta_vars"][last_object] = last_object[:_theta_node][:_theta]
end

function _add_complicating_variables!(
    optimizer::BendersAlgorithm{Plasmo.OptiGraph},
    last_object::Plasmo.OptiGraph,
    next_object::Plasmo.OptiGraph,
    add_slacks::Bool,
    slack_penalty::Real
)
    graph = optimizer.graph
    last_nodes = all_nodes(last_object)
    next_nodes = all_nodes(next_object)

    last_object_vars = all_variables(last_object)
    next_object_vars = all_variables(next_object)

    # Get incident edges between last and next objects
    last_object_optiedges = Plasmo.incident_edges(optimizer, graph, last_object)
    next_object_optiedges = Plasmo.incident_edges(optimizer, graph, next_object)

    complicating_edges = [e for e in next_object_optiedges if e in last_object_optiedges]
    comp_vars = NodeVariableRef[]

    # Map complicating variables to their owning node and vice versa
    node_to_var = Dict{Plasmo.OptiNode, Vector{NodeVariableRef}}()
    var_to_node = Dict{NodeVariableRef, Plasmo.OptiNode}()

    # Loop through each edge
    for (i, edge) in enumerate(complicating_edges)
        # Loop through each linkref on each edge
        for (j, link) in enumerate(all_constraints(edge))

            con_obj = constraint_object(link)
            vars = con_obj.func.terms.keys

            next_object_link_vars = [var for var in vars if source_graph(JuMP.owner_model(var)) == next_object]
            last_object_link_vars = [var for var in vars if source_graph(JuMP.owner_model(var)) == last_object]

            # Get the optinodes containing the next set of variables
            #next_optinode = optinode(next_object_vars[1]) #NEXT: Fix this!

            # Get the set of nodes in the next_object that are included in the constraint
            next_optinode = JuMP.owner_model(next_object_link_vars[1])

            # Map each complicating variable to the node where the
            # complicating variable copy will live
            for (k, var) in enumerate(last_object_link_vars)
                if !(var in comp_vars)
                    push!(comp_vars, var)
                end
                if haskey(var_to_node, var)
                    continue
                else
                    var_to_node[var] = next_optinode
                    if haskey(node_to_var, next_optinode)
                        push!(node_to_var[next_optinode], var)
                    else
                        node_to_var[next_optinode] = NodeVariableRef[]
                        push!(node_to_var[next_optinode], var)
                    end
                end
            end
        end
    end


    # Get the set of all complicating variables on last object
    comp_var_map = Dict{NodeVariableRef, Int}()

    # Map complicating variables to their index
    for (i, var) in enumerate(comp_vars)
        comp_var_map[var] = i
    end

    # Save the complicating variable data
    optimizer.comp_vars[next_object] = comp_vars
    optimizer.comp_var_map[next_object] = comp_var_map

    # Save structures used for cut data
    dual_iters = Matrix{Float64}(undef, length(comp_vars), 0)
    primal_iters = Matrix{Float64}(undef, length(comp_vars), 0)
    phis = Vector{Float64}[]
    phis_LR = Vector{Float64}[]

    # Map dual/primal/objective info to nodes
    optimizer.dual_iters[next_object] = dual_iters
    optimizer.primal_iters[next_object] = primal_iters
    optimizer.phis[next_object] = phis
    optimizer.phis_LR[next_object] = phis_LR

    ############## Add constraints for complicating variables #####################

    # Define map from variable copies to their original variables
    var_copy_map = Dict{NodeVariableRef, NodeVariableRef}()

    # Loop through the nodes and add a variable copy on each node
    for node in keys(node_to_var)
        node_comp_vars = node_to_var[node]

        varref_copy = @variable(node, _comp_vars_copy[1:length(node_comp_vars)])

        for (i, var) in enumerate(node_comp_vars)
            if JuMP.has_lower_bound(var)
                JuMP.set_lower_bound(varref_copy[i], lower_bound(var))
            end

            if JuMP.has_upper_bound(var)
                JuMP.set_upper_bound(varref_copy[i], upper_bound(var))
            end
            var_copy_map[var] = varref_copy[i]
        end
    end

    # Save the var_copy_map
    optimizer.var_copy_map[next_object] = var_copy_map

    # Define dictionary structures to map from constraints to nodes
    # con = normal constraints (does not require linking)
    # link = linking constriants (requires linking)
    # Linking constraints are required for constraints in the subgraph that have
    # complicating variables on more than one node
    con_to_node = Dict{ConstraintRef, Plasmo.OptiNode}()
    node_to_con = Dict{Plasmo.OptiNode, Vector{ConstraintRef}}()
    linking_to_node = Dict{ConstraintRef, Plasmo.OptiNode}()
    node_to_linking = Dict{Plasmo.OptiNode, Vector{ConstraintRef}}()

    # Loop through the complicating edges; decide if a constraint on the edge
    # requires a normal or a linking constraint
    for (i, edge) in enumerate(complicating_edges)
        for (j, link) in enumerate(all_constraints(edge))
            con_obj = constraint_object(link)
            vars = con_obj.func.terms.keys

            # Get the variables from the constraint in the next_object
            next_object_link_vars = [var for var in vars if source_graph(JuMP.owner_model(var)) == next_object]
            last_object_link_vars = [var for var in vars if source_graph(JuMP.owner_model(var)) == last_object]
            next_object_copy_vars = [var_copy_map[var] for var in last_object_link_vars]

            # Get the set of nodes in the next_object that are included in the constraint
            var_copy_optinodes = JuMP.owner_model.(next_object_copy_vars)
            next_object_optinodes = JuMP.owner_model.(next_object_link_vars)
            next_optinodes = unique(union(var_copy_optinodes, next_object_optinodes))

            # Get the first optinode in the set
            next_optinode = next_optinodes[1]

            # if the length of next_optinodes is 1, then use a normal constraint
            if length(next_optinodes) == 1
                con_to_node[link] = next_optinode
                if haskey(node_to_con, next_optinode)
                    push!(node_to_con[next_optinode], link)
                else
                    node_to_con[next_optinode] = ConstraintRef[]
                    push!(node_to_con[next_optinode], link)
                end

            # if the length of next_optinodes is > 1, use a linking constraint
            elseif length(next_optinodes) > 1
                linking_to_node[link] = next_optinode
                if haskey(node_to_linking, next_optinode)
                    push!(node_to_linking[next_optinode], link)
                else
                    node_to_linking[next_optinode] = ConstraintRef[]
                    push!(node_to_linking[next_optinode], link)
                end
            else
                error("Link does not connect to next object")
            end
        end
    end

    # Add slacks to the constraints if add_slacks i true
    if add_slacks
        # Loop through nodes and add a normal constraint for constraints with complicating variables
        for node in keys(node_to_con)
            cons = node_to_con[node]
            _add_slack_to_node(optimizer, next_object, node, length(cons), slack_penalty)

            for (j, con) in enumerate(cons)
                con_obj = constraint_object(con)
                _add_constraint_to_subproblem!(con_obj, comp_vars, var_copy_map,
                                              next_object, node, false; slack = true,
                                              slack_up = node[:_slack_up][j],
                                              slack_down = node[:_slack_down][j]
                )
            end
        end
        # Loop through nodes and add a linking constraint for constraints with complicating variables
        for node in keys(node_to_linking)
            links = node_to_linking[node]
            _add_slack_to_node_for_links(optimizer, next_object, node, length(links), slack_penalty)

            for (j, link) in enumerate(links)
                con_obj = constraint_object(link)
                _add_constraint_to_subproblem!(con_obj, comp_vars, var_copy_map,
                                              next_object, node, true; slack = true,
                                              slack_up = node[:_slack_up_link][j],
                                              slack_down = node[:_slack_down_link][j]
                )
            end
        end
    else
        # Loop through nodes and add a normal constraint for constraints with complicating variables
        for node in keys(node_to_con)
            cons = node_to_con[node]
            for (j, con) in enumerate(cons)
                con_obj = constraint_object(con)
                _add_constraint_to_subproblem!(con_obj, comp_vars, var_copy_map,
                                               next_object, node, false
                )
            end
        end
        # Loop through nodes and add a linking constraint for constraints with complicating variables
        for node in keys(node_to_linking)
            links = node_to_linking[node]
            for (j, link) in enumerate(links)
                con_obj = constraint_object(link)
                _add_constraint_to_subproblem!(con_obj, comp_vars, var_copy_map,
                                               next_object, node, true
                )
            end
        end
    end
end

function _add_next_object!(
    optimizer::BendersAlgorithm{Plasmo.OptiGraph},
    last_object::Plasmo.OptiGraph,
    next_object::Plasmo.OptiGraph,
    relaxed::Bool
)
    # Get the nodes from the last and next objects
    next_object_nodes = all_nodes(next_object)
    last_object_nodes = all_nodes(last_object)

    # Get the map from node to graph
    node_to_graph = optimizer.ext["node_to_graph"]

    # Get the incident edges for next object nodes and last object nodes
    last_object_incident_edges = incident_edges(optimizer, optimizer.graph, last_object)
    next_object_incident_edges = incident_edges(optimizer, optimizer.graph, next_object)

    # Get the edges that are incident to next object but not connected to last object
    next_object_edge_diff = setdiff(next_object_incident_edges, last_object_incident_edges)

    next_graphs = Plasmo.OptiGraph[]

    # Loop through the edges and get the subgraphs they map to
    for (i, edge) in enumerate(next_object_edge_diff)
        for (j, node) in enumerate(edge.nodes)
            connected_graph = node_to_graph[node]
            if (connected_graph != next_object) && !(connected_graph in next_graphs)
                push!(next_graphs, connected_graph)
            end
        end
    end

    optimizer.solve_order_dict[next_object] = next_graphs
    popfirst!(optimizer.ext["search_next"])

    if length(next_graphs) == 0
        return false
    end

    ############# Add Cost-to-Go Variable #################
    # Add cost-to-go variable and add to objective function
    _add_cost_to_go!(optimizer, next_object, relaxed)
    #Plasmo._init_graph_backend(optimizer.graph)

    for g in next_graphs
        if g in optimizer.solve_order
            error("The subproblems do not form a tree")
        end
        push!(optimizer.solve_order, g)
        optimizer.parent_objects[g] = next_object
        push!(optimizer.ext["search_next"], g)
    end
end

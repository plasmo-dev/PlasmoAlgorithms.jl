# Enable printing the optimizer
function Base.string(optimizer::BendersAlgorithm)
    return return @sprintf(
        """
        A BendersAlgorithm
        -------------------------------------------
        %32s %9s
        %32s %9s
        %32s %9s
        %32s %9s
        """,
        "Num subproblem subgraphs:",
        length(optimizer.solve_order),
        "MIP subproblems (nonroot):",
        optimizer.is_MIP,
        "Absolute Tolerance:",
        optimizer.tol,
        "Maximum Iterations:",
        optimizer.max_iters,
    )
end

Base.print(io::IO, optimizer::BendersAlgorithm) = Base.print(io, Base.string(optimizer))
Base.show(io::IO, optimizer::BendersAlgorithm) = Base.print(io, optimizer)

function _get_hyper_projection(optimizer::BendersAlgorithm, graph::Plasmo.OptiGraph)
    if haskey(optimizer.ext, "_projection")
        return optimizer.ext["_projection"]
    else
        proj = Plasmo.hyper_projection(graph)
        optimizer.ext["_projection"] = proj
        return proj
    end
end

function Plasmo.incident_edges(projection::T, subgraph::Plasmo.OptiGraph) where T <: Plasmo.GraphProjection
    return incident_edges(projection, all_nodes(subgraph))
end

function _get_incident_edges(optimizer::BendersAlgorithm, graph::Plasmo.OptiGraph, subgraph::Plasmo.OptiGraph)
    projection = _get_hyper_projection(optimizer, graph)
    edges = Plasmo.incident_edges(projection, subgraph)
    optimizer.ext["incident_edges"][subgraph] = edges
    return edges
end

function _get_incident_edges(optimizer::BendersAlgorithm, graph::Plasmo.RemoteOptiGraph, subgraph::Plasmo.RemoteOptiGraph)
    edges = Plasmo.incident_edges(subgraph)
    return edges
end

function Plasmo.incident_edges(optimizer::BendersAlgorithm, graph::T, subgraph::T) where {T <: Plasmo.AbstractOptiGraph}
    if haskey(optimizer.ext["incident_edges"], subgraph)
        return optimizer.ext["incident_edges"][subgraph]
    else
        return _get_incident_edges(optimizer, graph, subgraph)
    end
end

function _theta_value(optimizer::BendersAlgorithm, object::T) where {T <: Plasmo.AbstractOptiGraph}
    theta_var = optimizer.ext["theta_vars"][object]
    return sum(JuMP.value(object, theta) for theta in theta_var)
end

function _get_theta(optimizer::BendersAlgorithm, object::T) where {T <: Plasmo.AbstractOptiGraph}
    return optimizer.ext["theta_vars"][object]
end

function _get_theta(optimizer::BendersAlgorithm, object::T, idx::Int) where {T <: Plasmo.AbstractOptiGraph}
    return optimizer.ext["theta_vars"][object][idx]
end

# Define a function for getting the dual value of the linking constraints
function _get_next_duals(optimizer::BendersAlgorithm, next_object::T) where {T <: Plasmo.AbstractOptiGraph}
    var_copy_map = optimizer.var_copy_map[next_object]
    comp_vars = optimizer.comp_vars[next_object]

    copy_vars = [var_copy_map[comp_vars[i]] for i in 1:length(comp_vars)]
    dual_vector = _get_dual_vector(next_object, copy_vars)
    return dual_vector
end

# TODO: Need to improve the interface for getting values
# Define function for getting the best solution from the optimizer
function JuMP.value(
    optimizer::BendersAlgorithm{T},
    var::V,
    object::T
) where {T <: Plasmo.AbstractOptiGraph, V <: JuMP.AbstractVariableRef}
    var_solution_map = optimizer.var_solution_map[object]
    # Get the variable index for the variable on the node
    var_index = var_solution_map[var]
    # Get the vector of solutions for the node
    graph_values = optimizer.last_solutions[object]
    # Get the variable value from the vector of solutions

    var_value = graph_values[var_index]

    return var_value
end

"""
    JuMP.value(opt::BendersAlgorithm, var::VariableRef)
    
Returns the value of `var` from the BendersAlgorithm object. The value corresponds to the
best upper bound of the optimizer.
"""
function JuMP.value(optimizer::BendersAlgorithm, var::V) where {V <: JuMP.AbstractVariableRef}
    best_solutions = optimizer.best_solutions
    var_solution_map = optimizer.var_solution_map
    var_to_graph_map = optimizer.var_to_graph_map

    if !(var in keys(var_to_graph_map))
        error("$var is not defined in any of the subgraphs")
    end
    owning_object = var_to_graph_map[var]
    object_sols = best_solutions[owning_object]
    var_idx = var_solution_map[owning_object][var]

    var_value = object_sols[var_idx]

    return var_value
end

"""
    JuMP.value(opt::BendersAlgorithm, vars::Vector{VariableRef})
    
Returns a vector of variables contained in the `vars` vector from the BendersAlgorithm
object. The values correspond to the best upper bound of the optimizer.
"""
function JuMP.value(optimizer::BendersAlgorithm, vars::Vector{V}) where {V <: JuMP.AbstractVariableRef}
    best_solutions = optimizer.best_solutions
    var_solution_map = optimizer.var_solution_map
    var_to_graph_map = optimizer.var_to_graph_map

    value_vector = Vector{Float64}(undef, length(vars))

    for (i, var) in enumerate(vars)
        if !(var in keys(var_to_graph_map))
            error("$var is not defined in any of the subgraphs")
        end
        owning_object = var_to_graph_map[var]
        object_sols = best_solutions[owning_object]
        var_idx = var_solution_map[owning_object][var]

        var_value = object_sols[var_idx]

        value_vector[i] = var_value
    end

    return value_vector
end

function JuMP.objective_value(optimizer::BendersAlgorithm)
    if optimizer.current_iter == 0
        error("Optimize has not been called")
    end
    return optimizer.best_upper_bound
end

function JuMP.dual_objective_value(optimizer::BendersAlgorithm)
    if optimizer.current_iter == 0
        error("Optimize has not been called")
    end
    return optimizer.lower_bounds[end]
end

function relative_gap(optimizer::BendersAlgorithm)
    if optimizer.current_iter == 0
        error("Optimize has not been called")
    end
    lb = optimizer.lower_bounds[end]
    ub = optimizer.best_upper_bound

    gap = abs((ub - lb)/ub)

    return gap
end

# Define function for warm starting
function _warm_start(optimizer::BendersAlgorithm{Plasmo.OptiGraph}, object::Plasmo.OptiGraph)
    all_vars = all_variables(object)
    all_vals = optimizer.best_solutions[object]

    # Set the start value for each variable
    for j in 1:length(all_vars)
        JuMP.set_start_value(all_vars[j], all_vals[j])
    end
end

function _warm_start(optimizer::BendersAlgorithm{Plasmo.RemoteOptiGraph}, object::Plasmo.RemoteOptiGraph)
    all_vals = optimizer.best_solutions[object]
    darray = rgraph.graph

    f = @spawnat object.worker begin
        lg = Plasmo.local_graph(darray)
        all_vars = all_variables(lg)
        JuMP.set_start_value.(all_vars, all_vals)
    end
end

# Define functions for adding slacks, and adding them to the objectives
function _add_slack_to_node(optimizer::BendersAlgorithm, next_object, node::N, num_links, slack_penalty) where {N <: Union{Plasmo.OptiNode, Plasmo.RemoteNodeRef}}
    # Define slack variables
    @variable(node, _slack_up[1:num_links] >= 0)
    @variable(node, _slack_down[1:num_links] >= 0)
    V = Plasmo.variable_type(optimizer.graph)
    slack_vars_dict = optimizer.slack_vars

    if !(haskey(slack_vars_dict, next_object))
        slack_vars_dict[next_object] = V[]
    end

    for var in _slack_up
        push!(slack_vars_dict[next_object], var)
    end
    for var in _slack_down
        push!(slack_vars_dict[next_object], var)
    end

    # Get the objective function
    obj_func = objective_function(node)

    # Ensure the objective is an Affine Expression
    if isa(obj_func, V)
        obj_func = GenericAffExpr{Float64, N}(0, obj_func => 1)
    elseif typeof(obj_func) == nothing
        obj_func = GenericAffExpr{Float64, N}()
    end

    # Add the slacks to the objective function
    for i in 1:num_links
        JuMP.add_to_expression!(obj_func, slack_penalty * _slack_up[i])
        JuMP.add_to_expression!(obj_func, slack_penalty * _slack_down[i])
    end

    # Reset the objective so it has the slacks
    #JuMP.set_objective_function(node, obj_func)
    @objective(node, Min, obj_func)
end

function _add_slack_to_node_for_links(optimizer::BendersAlgorithm, next_object::T, node::N, num_links, slack_penalty) where {T <: Plasmo.AbstractOptiGraph, N <: Union{Plasmo.OptiNode, Plasmo.RemoteNodeRef}}
    V = Plasmo.variable_type(optimizer.graph)
    # Define slack variables
    @variable(node, _slack_up_link[1:num_links] >= 0)
    @variable(node, _slack_down_link[1:num_links] >= 0)

    slack_vars_dict = optimizer.slack_vars

    if !(haskey(slack_vars_dict, next_object))
        slack_vars_dict[next_object] = V[]
    end

    for var in _slack_up_link
        push!(slack_vars_dict[next_object], var)
    end
    for var in _slack_down_link
        push!(slack_vars_dict[next_object], var)
    end

    # Get the objective function
    obj_func = objective_function(node)

    # Ensure the objective is an Affine Expression
    if isa(obj_func, V)
        obj_func = GenericAffExpr{Float64, V}(0, obj_func => 1)
    elseif isnothing(obj_func)
        obj_func = GenericAffExpr{Float64, V}()
    end

    # Add the slacks to the objective function
    for i in 1:num_links
        JuMP.add_to_expression!(obj_func, slack_penalty * _slack_up_link[i])
        JuMP.add_to_expression!(obj_func, slack_penalty * _slack_down_link[i])
    end

    # Reset the objective so it has the slacks
    #JuMP.set_objective_function(node, obj_func)
    @objective(node, Min, obj_func)

end

# Test if the problem is a MIP
function _set_is_MIP(optimizer::BendersAlgorithm{Plasmo.OptiGraph})
    graph = optimizer.graph

    graph_vars = setdiff(JuMP.all_variables(graph), JuMP.all_variables(optimizer.root_object))

    for var in graph_vars
        if JuMP.is_binary(var) || JuMP.is_integer(var)
            optimizer.is_MIP = true
            return nothing
        end
    end
    optimizer.is_MIP = false
    return nothing
end

# Test if the problem is a MIP
function _set_is_MIP(optimizer::BendersAlgorithm{Plasmo.RemoteOptiGraph})
    graph = optimizer.graph

    subgraphs = graph.subgraphs
    for g in subgraphs
        if PlasmoBenders._check_is_MIP(g)
            optimizer.is_MIP = true
            return nothing
        end
    end
    optimizer.is_MIP = false
    return nothing
end

function _check_is_MIP(rgraph::Plasmo.RemoteOptiGraph)
    darray = rgraph.graph
    f = @spawnat rgraph.worker begin
        lg = Plasmo.local_graph(darray)
        vars = all_variables(lg)
        bool_val = any(JuMP.is_binary.(vars)) || any(JuMP.is_integer.(vars))
        bool_val
    end
    return fetch(f)
end

function _get_objects(optimizer::BendersAlgorithm{T}) where {T <: Plasmo.AbstractOptiGraph}
    return local_subgraphs(optimizer.graph)
end

function _get_objects(graph::T) where {T <: Plasmo.AbstractOptiGraph}
    return local_subgraphs(graph)
end

function _check_termination_status(optimizer::BendersAlgorithm, object, count)
    if termination_status(object) == MOI.INFEASIBLE
        if get_add_slacks(optimizer)
            if haskey(optimizer.slack_vars, object)
                slack_vars = optimizer.slack_vars[object]
                if !(get_fix_slacks(optimizer))
                    error("Model on node/graph $count is infeasible; `add_slacks` is true, but the model is still infeasible!")
                else
                    for var in slack_vars
                        JuMP.unfix(var)
                        JuMP.set_lower_bound(var, 0)
                    end
                    @warn("Model on node/graph $count is infeasible; unfixing slack variables")

                    optimize!(object)

                    if termination_status(object) != MOI.OPTIMAL
                        error("Model on node/graph $count still did not reach optimal solution; status = $(termination_status(object))")
                    end
                end
            else
                error("Model on node/graph $count is infeasible; `add_slacks` is true, but there are not slacks to unfix!")
            end
        elseif get_feasibility_cuts(optimizer) && count != 1
            return false
        else
            error("Model on node/graph $count is infeasible")
        end
    elseif termination_status(object) == MOI.INFEASIBLE_OR_UNBOUNDED
        error("Model on node/graph $count is unbounded or infeasible")
    elseif (termination_status(object) != MOI.OPTIMAL) && (termination_status(object) != MOI.LOCALLY_SOLVED) && (termination_status(object) != MOI.TIME_LIMIT)
        error("Model on node/graph $count terminated with status $(termination_status(object))")
    end
    return true
end

function JuMP.objective_value(optimizer::BendersAlgorithm, object, last_obj::Bool)
    if optimizer.parallelize_Benders && last_obj
        subs = get_subgraphs(object)
        if length(subs) > 0
            obj_val = [0.]
            for sub in subs
                obj_val[1] += JuMP.objective_value(sub)
            end
            return obj_val[1]
        else
            return JuMP.objective_value(object)
        end
    else
        return JuMP.objective_value(object)
    end
end

function _check_fixed_slacks!(optimizer::BendersAlgorithm, object)
    if get_add_slacks(optimizer)
        if get_fix_slacks(optimizer)
            slack_vars = optimizer.slack_vars
            if haskey(slack_vars, object)
                if !(JuMP.is_fixed(slack_vars[object][1]))
                    for var in slack_vars[object]
                        JuMP.fix(var, 0, force = true)
                    end
                end
            end
        end
    end
end

function _check_two_stage_tree(optimizer::BendersAlgorithm)
    if get_parallelize_benders(optimizer) || get_regularize(optimizer) || get_feasibility_cuts(optimizer)
        graph = optimizer.graph
        root_object = optimizer.root_object
        solve_order = optimizer.solve_order

        next_objects = optimizer.solve_order_dict[root_object]
        if length(next_objects) != (length(solve_order) - 1)
            unconnected = length(solve_order) - length(next_objects) - 1
            if get_parallelize_benders(optimizer)
                error("There are $(length(solve_order) - 1) subproblem objects, but $unconnected are not connected to start object")
            elseif get_regularize(optimizer)
                error(
                    "Regularization only works for Benders type problems; " *
                    "There are $(length(solve_order) - 1) subproblem objects, but $unconnected are not connected to start object"
                )
            end
        end
    end
end

_options_bool_fields = [
    :strengthened,
    :multicut,
    :feasibility_cuts,
    :regularize,
    :parallelize_benders,
    :parallelize_forward,
    :parallelize_backward,
    :add_slacks,
    :fix_slacks,
    :warm_start,
    :relaxed_init_cuts
]

_options_real_fields = [
    :slack_penalty,
    :regularize_param
]

_algorithm_fields = [
    :graph,
    :root_object,
    :is_MIP,
    :max_iters,
    :time_limit,
    :tol,
    :time_forward_pass,
    :time_backward_pass,
    :time_init,
    :time_iterations,
    :lower_bounds,
]

for field in union(_options_bool_fields, _options_real_fields)
    method = Symbol("get_", field)
    @eval begin
        @doc """
            $($method)(optimizer::BendersAlgorithm)
        Return the value of $($(QuoteNode(field))) from the `options` field of the `BendersAlgorithm`
        """
        $method(optimizer::BendersAlgorithm) = getproperty(optimizer.options, $(QuoteNode(field)))
    end
    @eval export $method
end

for field in _algorithm_fields
    method = Symbol("get_", field)
    @eval begin
        @doc """
            $($method)(optimizer::BendersAlgorithm)
        Return the value of the $($(QuoteNode(field))) attribute from the `BendersAlgorithm`
        """
        $method(optimizer::BendersAlgorithm) = getproperty(optimizer, $(QuoteNode(field)))
    end
    @eval export $method
end

@doc """
    get_upper_bounds(optimizer::BendersAlgorithm; monotonic = true)
Return the value of `upper_bounds` from the `BendersAlgorithm` object. This is a vector of the upper bounds at each iteration. The upper bound returned by each iteration is not guaranteed to decrease monotonically. If `monotonic=true`, this function returns the best upper bound seen up to each iteration (rather than the value of the feasible solution seen at that iteration)
"""
function get_upper_bounds(optimizer::BendersAlgorithm; monotonic = true)
    ubs = optimizer.upper_bounds
    if monotonic
        if optimizer.current_iter == 0
            return ubs
        else
            return [minimum(ubs[1:i] for i in 1:length(ubs))]
        end
    else
        return ubs
    end
end

for field in _options_bool_fields
    method = Symbol("set_", field, "!")
    @eval begin
        @doc """
            $($method)(optimizer::BendersAlgorithm, val::Bool)
        Set the value of $($(QuoteNode(field))) from the `options` field of the `BendersAlgorithm` to `val`
        """
        $method(optimizer::BendersAlgorithm, val::Bool) = optimizer.options.$field = val
    end
    @eval export $method
end

for field in _options_real_fields
    method = Symbol("set_", field, "!")
    @eval begin
        @doc """
            $($method)(optimizer::BendersAlgorithm, val::Real)
        Set the value of $($(QuoteNode(field))) from the `options` field of the `BendersAlgorithm` to `val`
        """
        $method(optimizer::BendersAlgorithm, val::Real) = optimizer.options.$field = val
    end
    @eval export $method
end

function _get_binary_bool_vector(object::RemoteOptiGraph)
    darray = object.graph
    f = @spawnat object.worker begin
        lg = Plasmo.local_graph(darray)
        PlasmoBenders._get_binary_bool_vector(lg)
    end
    return fetch(f)
end

function _get_integer_bool_vector(object::RemoteOptiGraph)
    darray = object.graph
    f = @spawnat object.worker begin
        lg = Plasmo.local_graph(darray)
        PlasmoBenders._get_integer_bool_vector(lg)
    end
    return fetch(f)
end

function _get_binary_bool_vector(object::OptiGraph)
    all_vars = all_variables(object)
    return JuMP.is_binary.(all_vars)
end

function _get_integer_bool_vector(object::OptiGraph)
    all_vars = all_variables(object)
    return JuMP.is_integer.(all_vars)
end

function _get_object_last_solutions(object::OptiGraph)
    vars = JuMP.all_variables(object)
    return [JuMP.value(object, var) for var in vars]
end

function _get_object_last_solutions(object::RemoteOptiGraph)
    darray = object.graph
    f = @spawnat object.worker begin
        lg = Plasmo.local_graph(darray)
        PlasmoBenders._get_object_last_solutions(lg)
    end
    return fetch(f)
end

function _get_variable_values(object::OptiGraph, variables::Vector{Plasmo.NodeVariableRef})
    return [value(object, var) for var in variables]
end

function _get_variable_values(object::RemoteOptiGraph, variables::Vector{Plasmo.RemoteVariableRef})
    darray = object.graph
    pvars = [Plasmo._convert_remote_to_proxy(object, var) for var in variables]
    f = @spawnat object.worker begin
        lg = Plasmo.local_graph(darray)
        lvars = [Plasmo._convert_proxy_to_local(lg, var) for var in pvars]
        PlasmoBenders._get_variable_values(lg, lvars)
    end
    return fetch(f)
end

function _get_dual_vector(object::OptiGraph, variables::Vector{Plasmo.NodeVariableRef})
    return Float64[JuMP.dual(object, FixRef(var)) for var in variables]
end

function _get_dual_vector(object::RemoteOptiGraph, variables::Vector{Plasmo.RemoteVariableRef})
    darray = object.graph
    pvars = [Plasmo._convert_remote_to_proxy(object, var) for var in variables]
    f = @spawnat object.worker begin
        lg = Plasmo.local_graph(darray)
        lvars = [Plasmo._convert_proxy_to_local(lg, var) for var in pvars]
        PlasmoBenders._get_dual_vector(lg, lvars)
    end
    return fetch(f)
end

function _fix_variables(object::OptiGraph, variables::Vector{Plasmo.NodeVariableRef}, values::Vector{Float64})
    JuMP.fix.(variables, values, force = true)
end

function _fix_variables(object::RemoteOptiGraph, variables::Vector{Plasmo.RemoteVariableRef}, values::Vector{Float64})
    darray = object.graph
    pvars = [Plasmo._convert_remote_to_proxy(object, var) for var in variables]
    f = @spawnat object.worker begin
        lg = Plasmo.local_graph(darray)
        lvars = [Plasmo._convert_proxy_to_local(lg, var) for var in pvars]

        PlasmoBenders._fix_variables(lg, lvars, values)
    end
    return nothing
end

function _unfix_variables(object::OptiGraph, variables::Vector{Plasmo.NodeVariableRef})
    JuMP.unfix.(variables)
end

function _unfix_variables(object::RemoteOptiGraph, variables::Vector{Plasmo.RemoteVariableRef})
    darray = object.graph
    pvars = [Plasmo._convert_remote_to_proxy(object, var) for var in variables]
    f = @spawnat object.worker begin
        lg = Plasmo.local_graph(darray)
        lvars = [Plasmo._convert_proxy_to_local(lg, var) for var in pvars]

        PlasmoBenders._unfix_variables(lg, lvars)
    end
    return nothing
end
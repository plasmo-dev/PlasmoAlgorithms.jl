function Plasmo.incident_edges(graph::Plasmo.OptiGraph, subgraph::Plasmo.OptiGraph)
    projection = hyper_projection(graph)
    return incident_edges(projection, all_nodes(subgraph))
end

function _theta_value(optimizer::BendersOptimizer, object::Plasmo.OptiGraph)
    theta_var = optimizer.ext["theta_vars"][object]
    return sum(JuMP.value(object, theta) for theta in theta_var)
end

function _get_theta(optimizer::BendersOptimizer, object::Plasmo.OptiGraph)
    return optimizer.ext["theta_vars"][object]
end

function _get_theta(optimizer::BendersOptimizer, object::Plasmo.OptiGraph, idx::Int)
    return optimizer.ext["theta_vars"][object][idx]
end

# Define a function for getting the dual value of the linking constraints
function _get_next_duals(optimizer::BendersOptimizer, next_object::Plasmo.OptiGraph)
    var_copy_map = optimizer.var_copy_map[next_object]
    comp_vars = optimizer.comp_vars[next_object]

    dual_vector = Float64[JuMP.dual(next_object, FixRef(var_copy_map[comp_vars[i]])) for i in 1:length(comp_vars)]
    return dual_vector
end

# TODO: Need to improve the interface for getting values
# Define function for getting the best solution from the optimizer
function JuMP.value(
    optimizer::BendersOptimizer{Plasmo.OptiGraph},
    var::NodeVariableRef,
    object::Plasmo.OptiGraph
)
    var_solution_map = optimizer.var_solution_map[object]
    # Get the variable index for the variable on the node
    var_index = var_solution_map[var]
    # Get the vector of solutions for the node
    graph_values = optimizer.last_solutions[object]
    # Get the variable value from the vector of solutions

    var_value = graph_values[var_index]

    return var_value
end

# Define function for warm starting
function _warm_start(optimizer::BendersOptimizer, object)
    all_vars = all_variables(object)
    all_vals = optimizer.best_solutions[object]

    # Set the start value for each variable
    for j in 1:length(all_vars)
        JuMP.set_start_value(all_vars[j], all_vals[j])
    end
end

# Define functions for adding slacks, and adding them to the objectives
function _add_slack_to_node(optimizer::BendersOptimizer, next_object, node::Plasmo.OptiNode, num_links, slack_penalty)
    # Define slack variables
    @variable(node, _slack_up[1:num_links] >= 0)
    @variable(node, _slack_down[1:num_links] >= 0)

    slack_vars_dict = optimizer.slack_vars

    if !(haskey(slack_vars_dict, next_object))
        slack_vars_dict[next_object] = NodeVariableRef[]
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
    if typeof(obj_func) == NodeVariableRef
        obj_func = AffExpr(0, obj_func => 1)
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

function _add_slack_to_node_for_links(optimizer::BendersOptimizer, next_object::Plasmo.OptiGraph, node::Plasmo.OptiNode, num_links, slack_penalty)
    # Define slack variables
    @variable(node, _slack_up_link[1:num_links] >= 0)
    @variable(node, _slack_down_link[1:num_links] >= 0)

    slack_vars_dict = optimizer.slack_vars

    if !(haskey(slack_vars_dict, next_object))
        slack_vars_dict[next_object] = NodeVariableRef[]
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
    if typeof(obj_func) == NodeVariableRef
        obj_func = AffExpr(0, obj_func => 1)
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
function _set_is_MIP(optimizer::BendersOptimizer)
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

function _get_objects(optimizer::BendersOptimizer{Plasmo.OptiGraph})
    return getsubgraphs(optimizer.graph)
end

function _get_objects(graph::Plasmo.OptiGraph)
    return getsubgraphs(graph)
end

function _check_termination_status(optimizer::BendersOptimizer, object, count)
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
        else
            error("Model on node/graph $count is infeasible")
        end
    elseif termination_status(object) == MOI.INFEASIBLE_OR_UNBOUNDED
        error("Model on node/graph $count is unbounded or infeasible")
    elseif (termination_status(object) != MOI.OPTIMAL) && (termination_status(object) != MOI.LOCALLY_SOLVED)
        error("Model on node/graph $count terminated with status $(termination_status(object))")
    end
end

function JuMP.objective_value(optimizer::BendersOptimizer, object, last_obj::Bool)
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

function _check_fixed_slacks!(optimizer::BendersOptimizer, object)
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

function _check_parallelize_Benders(optimizer::BendersOptimizer)
    if get_parallelize_benders(optimizer) || get_regularize(optimizer)
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

for field in union(_options_bool_fields, _options_real_fields)
    method = Symbol("get_", field)
    @eval begin
        @doc """
            $($method)(optimizer::BendersOptimizer)
        Return the value of $($(QuoteNode(field))) from the `options` field of the `BendersOptimizer`
        """
        $method(optimizer::BendersOptimizer) = getproperty(optimizer.options, $(QuoteNode(field)))
    end
    @eval export $method
end

for field in _options_bool_fields
    method = Symbol("set_", field, "!")
    @eval begin
        @doc """
            $($method)(optimizer::BendersOptimizer, val::Bool)
        Set the value of $($(QuoteNode(field))) from the `options` field of the `BendersOptimizer` to `val`
        """
        $method(optimizer::BendersOptimizer, val::Bool) = optimizer.options.$field = val
    end
    @eval export $method
end

for field in _options_real_fields
    method = Symbol("set_", field, "!")
    @eval begin
        @doc """
            $($method)(optimizer::BendersOptimizer, val::Real)
        Set the value of $($(QuoteNode(field))) from the `options` field of the `BendersOptimizer` to `val`
        """
        $method(optimizer::BendersOptimizer, val::Real) = optimizer.options.$field = val
    end
    @eval export $method
end

# TODO: Extend JuMP functions like termination status, objective value, etc.

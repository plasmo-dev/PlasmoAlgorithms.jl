function _add_constraint_to_subproblem!(
    con_obj::ScalarConstraint{GenericAffExpr{Float64, Plasmo.NodeVariableRef}, MOI.LessThan{Float64}},#ConstraintRef{Plasmo.OptiEdge{Plasmo.OptiGraph}, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}}},
    comp_vars,
    var_copy_map,
    next_object::Plasmo.OptiGraph,
    next_node::Plasmo.OptiNode,
    link::Bool;
    slack::Bool = false,
    slack_up = nothing,
    slack_down = nothing
)
    # Create empty expression (for a constraint)
    new_expr = GenericAffExpr{Float64, Plasmo.NodeVariableRef}()

    # Add variables or variable copy to the expression
    for (i, var) in enumerate(con_obj.func.terms.keys)
        if var in comp_vars
            add_to_expression!(new_expr, var_copy_map[var] * con_obj.func.terms[var])
        else
            add_to_expression!(new_expr, var * con_obj.func.terms[var])
        end
    end

    # If slacks are required, add them to the expression
    if slack
        add_to_expression!(new_expr, slack_up)
        add_to_expression!(new_expr, - slack_down)
    end

    # If link, make it a link constraint on next object
    if link
        @linkconstraint(next_object, new_expr <= con_obj.set.upper)
    # Otherwise, make it a normal constraint on a node
    else
        @constraint(next_node, new_expr <= con_obj.set.upper)
    end
end

function _add_constraint_to_subproblem!(
    con_obj::ScalarConstraint{GenericAffExpr{Float64, Plasmo.NodeVariableRef}, MOI.EqualTo{Float64}},#ConstraintRef{Plasmo.OptiEdge{Plasmo.OptiGraph}, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}}},
    comp_vars,
    var_copy_map,
    next_object::Plasmo.OptiGraph,
    next_node::Plasmo.OptiNode,
    link::Bool;
    slack::Bool = false,
    slack_up = nothing,
    slack_down = nothing
)
    # Create empty expression (for a constraint)
    new_expr = GenericAffExpr{Float64, Plasmo.NodeVariableRef}()

    # Add variables or variable copy to the expression
    for (i, var) in enumerate(con_obj.func.terms.keys)
        if var in comp_vars
            add_to_expression!(new_expr, var_copy_map[var] * con_obj.func.terms[var])
        else
            add_to_expression!(new_expr, var * con_obj.func.terms[var])
        end
    end

    # If slacks are required, add them to the expression
    if slack
        add_to_expression!(new_expr, slack_up)
        add_to_expression!(new_expr, - slack_down)
    end

    # If link, make it a link constraint on next object
    if link
        @linkconstraint(next_object, new_expr == con_obj.set.value)
    # Otherwise, make it a normal constraint on a node
    else
        @constraint(next_node, new_expr == con_obj.set.value)
    end
end

function _add_constraint_to_subproblem!(
    con_obj::ScalarConstraint{GenericAffExpr{Float64, Plasmo.NodeVariableRef}, MOI.GreaterThan{Float64}},#ConstraintRef{Plasmo.OptiEdge{Plasmo.OptiGraph}, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}}},
    comp_vars,
    var_copy_map,
    next_object::Plasmo.OptiGraph,
    next_node::Plasmo.OptiNode,
    link::Bool;
    slack::Bool = false,
    slack_up = nothing,
    slack_down = nothing
)
    # Create empty expression (for a constraint)
    new_expr = GenericAffExpr{Float64, Plasmo.NodeVariableRef}()

    # Add variables or variable copy to the expression
    for (i, var) in enumerate(con_obj.func.terms.keys)
        if var in comp_vars
            add_to_expression!(new_expr, var_copy_map[var] * con_obj.func.terms[var])
        else
            add_to_expression!(new_expr, var * con_obj.func.terms[var])
        end
    end

    # If slacks are required, add them to the expression
    if slack
        add_to_expression!(new_expr, slack_up)
        add_to_expression!(new_expr, - slack_down)
    end

    # If link, make it a link constraint on next object
    if link
        @linkconstraint(next_object, new_expr >= con_obj.set.lower)
    # Otherwise, make it a normal constraint on a node
    else
        @constraint(next_node, new_expr >= con_obj.set.lower)
    end
end

# Define function to add theta to the objective
function _add_theta_to_objective!(
    obj_func::AffExpr,
    node::Plasmo.OptiNode,
    _theta::Vector{NodeVariableRef}
)
    for i in 1:length(_theta)
        add_to_expression!(obj_func, _theta[i] * 1)
    end
    JuMP.set_objective_function(node, obj_func)

end

function _add_theta_to_objective!(
    obj_func::NodeVariableRef,
    node::Plasmo.OptiNode,
    _theta::Vector{NodeVariableRef}
)
    obj_func = obj_func + sum(_theta[i] for i in 1:length(_theta))
    JuMP.set_objective_function(node, obj_func)
end

function _set_theta_lower_bound!(object::Plasmo.OptiGraph, M)
    theta = object[:_theta_node][:_theta]
    for i in 1:length(theta)
        JuMP.set_lower_bound(object[:_theta_node][:_theta][i], M)
    end
end

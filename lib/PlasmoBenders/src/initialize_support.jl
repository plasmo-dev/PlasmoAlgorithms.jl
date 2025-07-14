function _add_constraint_set_to_subproblem!(
    optimizer::BendersAlgorithm{T}, 
    cons::Vector{ConstraintRef}, 
    next_object::T, 
    node::N, 
    comp_vars::Vector{V}, 
    var_copy_map::Dict, 
    link::Bool, 
    slack::Bool
) where {V <: JuMP.AbstractVariableRef, T <: Plasmo.AbstractOptiGraph, N <: Union{Plasmo.OptiNode, Plasmo.RemoteNodeRef}}
    if get_add_slacks(optimizer)
        _add_slack_to_node(optimizer, next_object, node, length(cons), slack_penalty)
        for (j, con) in enumerate(cons)
            con_obj = JuMP.constraint_object(con)
            _add_constraint_to_subproblem!(con_obj, comp_vars, var_copy_map,
                                          next_object, node, link; slack = slack,
                                          slack_up = node[:_slack_up][j],
                                          slack_down = node[:_slack_down][j]
            )
        end
    else
        for (j, con) in enumerate(cons)
            con_obj = JuMP.constraint_object(con)
                _add_constraint_to_subproblem!(con_obj, comp_vars, var_copy_map,
                                              next_object, node, link; slack = slack,
                                              slack_up = nothing,
                                              slack_down = nothing
                )
        end
    end
    return nothing
end

# function _add_constraint_set_to_subproblem!(
#     optimizer::BendersAlgorithm{Plasmo.RemoteOptiGraph}, 
#     cons::Vector{ConstraintRef}, 
#     next_object::Plasmo.RemoteOptiGraph, 
#     node::Plasmo.RemoteNodeRef, 
#     comp_vars::Vector{Plasmo.RemoteVariableRef}, 
#     var_copy_map::Dict, 
#     link::Bool, 
#     slack::Bool
# )
#     if get_add_slacks(optimizer)
#         _add_slack_to_node(optimizer, next_object, node, length(cons), slack_penalty)
#         for (j, con) in enumerate(cons)
#             con_obj = JuMP.constraint_object(con)
#             _add_constraint_to_subproblem!(con_obj, comp_vars, var_copy_map,
#                                           next_object, node, link; slack = slack,
#                                           slack_up = node[:_slack_up][j],
#                                           slack_down = node[:_slack_down][j]
#             )
#         end
#     else
#         con_objs = JuMP.constraint_object.(cons)
#         exprs = GenericAffExpr[]
#         constants = Number[]
#         funcs = []
#         for con in con_objs
#             new_expr = GenericAffExpr{Float64, Plasmo.RemoteVariableRef}()

#             for (i, var) in enumerate(con.func.terms.keys)
#                 if var in comp_vars
#                     add_to_expression!(new_expr, var_copy_map[var] * con.func.terms[var])
#                 else
#                     add_to_expression!(new_expr, var * con.func.terms[var])
#                 end
#             end #TODO: Make this its own function?
#             push!(exprs, new_expr)
#             push!(constants, _get_function_constant(con))
#             push!(funcs, con.set)
#         end

#         darray = next_object.graph
#         pnode = Plasmo._convert_remote_to_proxy(next_object, node)
#         #pexpr = Plasmo._convert_remote_to_proxy(next_object, expr)
#         f = @spawnat next_object.worker begin
#             lgraph = Plasmo.local_graph(darray)
#             if !link 
#                 lnode = Plasmo._convert_proxy_to_local(lgraph, pnode)
#             end
#             for (i, expr) in enumerate(exprs)
#                 lexpr = Plasmo._convert_remote_to_local(lgraph, pexpr)
#                 if link
#                     if isa(funcs[i], MOI.GreaterThan)
#                         @linkconstraint(lgraph, lexpr >= constants[i])
#                     elseif isa(funcs[i], MOI.LessThan)
#                         @linkconstraint(lgraph, lexpr <= constants[i])
#                     elseif isa(funcs[i], MOI.EqualTo)
#                         @linkconstraint(lgraph, lexpr == constants[i])
#                     end
#                 else
#                     if isa(funcs[i], MOI.GreaterThan)
#                         @constraint(lnode, lexpr >= constants[i])
#                     elseif isa(funcs[i], MOI.LessThan)
#                         @constraint(lnode, lexpr <= constants[i])
#                     elseif isa(funcs[i], MOI.EqualTo)
#                         @constraint(lnode, lexpr == constants[i])
#                     end
#                 end
#             end
#         end
#     end
#     return nothing
# end

function _get_function_constant(con::ScalarConstraint)
    if isa(con.set, MOI.GreaterThan)
        return con.set.lower
    elseif isa(con.set, MOI.LessThan)
        return con.set.upper
    elseif isa(con.set, MOI.EqualTo)
        return con.set.value
    else 
        error("Constraint is of type $(typeof(con.set)) which is not implemented")
    end
end

function _add_constraint_to_subproblem!(
    con_obj::ScalarConstraint{GenericAffExpr{Float64, V}, MOI.LessThan{Float64}},#ConstraintRef{Plasmo.OptiEdge{Plasmo.OptiGraph}, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}}},
    comp_vars,
    var_copy_map,
    next_object::T,
    next_node::N,
    link::Bool;
    slack::Bool = false,
    slack_up = nothing,
    slack_down = nothing
) where {V <: JuMP.AbstractVariableRef, T <: Plasmo.AbstractOptiGraph, N <: Union{Plasmo.OptiNode, Plasmo.RemoteNodeRef}}
    # Create empty expression (for a constraint)
    new_expr = GenericAffExpr{Float64, V}()

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
    con_obj::ScalarConstraint{GenericAffExpr{Float64, V}, MOI.EqualTo{Float64}},#ConstraintRef{Plasmo.OptiEdge{Plasmo.OptiGraph}, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}}},
    comp_vars,
    var_copy_map,
    next_object::T,
    next_node::N,
    link::Bool;
    slack::Bool = false,
    slack_up = nothing,
    slack_down = nothing
) where {V <: JuMP.AbstractVariableRef, T <: Plasmo.AbstractOptiGraph, N <: Union{Plasmo.OptiNode, Plasmo.RemoteNodeRef}}
    # Create empty expression (for a constraint)
    new_expr = GenericAffExpr{Float64, V}()

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
    con_obj::ScalarConstraint{GenericAffExpr{Float64, V}, MOI.GreaterThan{Float64}},#ConstraintRef{Plasmo.OptiEdge{Plasmo.OptiGraph}, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}}},
    comp_vars,
    var_copy_map,
    next_object::T,
    next_node::N,
    link::Bool;
    slack::Bool = false,
    slack_up = nothing,
    slack_down = nothing
) where {V <: JuMP.AbstractVariableRef, T <: Plasmo.AbstractOptiGraph, N <: Union{Plasmo.OptiNode, Plasmo.RemoteNodeRef}}
    # Create empty expression (for a constraint)
    new_expr = GenericAffExpr{Float64, V}()

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

################ QUADRATICS: 
function _add_constraint_to_subproblem!(
    con_obj::ScalarConstraint{GenericQuadExpr{Float64, V}, MOI.LessThan{Float64}},#ConstraintRef{Plasmo.OptiEdge{Plasmo.OptiGraph}, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}}},
    comp_vars,
    var_copy_map,
    next_object::T,
    next_node::N,
    link::Bool;
    slack::Bool = false,
    slack_up = nothing,
    slack_down = nothing
) where {V <: JuMP.AbstractVariableRef, T <: Plasmo.AbstractOptiGraph, N <: Union{Plasmo.OptiNode, Plasmo.RemoteNodeRef}}
    # Create empty expression (for a constraint)
    new_expr = GenericQuadExpr{Float64, V}()

    # Add variables or variable copy to the expression
    for (i, var) in enumerate(con_obj.func.aff.terms.keys)
        if var in comp_vars
            add_to_expression!(new_expr, var_copy_map[var] * con_obj.func.aff.terms[var])
        else
            add_to_expression!(new_expr, var * con_obj.func.aff.terms[var])
        end
    end

    for (i, pair) in enumerate(con_obj.func.terms.keys)
        var1 = pair.a
        var2 = pair.b

        expr_var1 = var1 in comp_vars ? var_copy_map[var1] : var1
        expr_var2 = var2 in comp_vars ? var_copy_map[var2] : var2
        
        add_to_expression!(new_expr, expr_var1 * expr_var2 * con_obj.func.terms[pair])
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
    con_obj::ScalarConstraint{GenericQuadExpr{Float64, V}, MOI.EqualTo{Float64}},#ConstraintRef{Plasmo.OptiEdge{Plasmo.OptiGraph}, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}}},
    comp_vars,
    var_copy_map,
    next_object::T,
    next_node::N,
    link::Bool;
    slack::Bool = false,
    slack_up = nothing,
    slack_down = nothing
) where {V <: JuMP.AbstractVariableRef, T <: Plasmo.AbstractOptiGraph, N <: Union{Plasmo.OptiNode, Plasmo.RemoteNodeRef}}
    # Create empty expression (for a constraint)
    new_expr = GenericQuadExpr{Float64, V}()

    # Add variables or variable copy to the expression
    for (i, var) in enumerate(con_obj.func.aff.terms.keys)
        if var in comp_vars
            add_to_expression!(new_expr, var_copy_map[var] * con_obj.func.aff.terms[var])
        else
            add_to_expression!(new_expr, var * con_obj.func.aff.terms[var])
        end
    end

    for (i, pair) in enumerate(con_obj.func.terms.keys)
        var1 = pair.a
        var2 = pair.b

        expr_var1 = var1 in comp_vars ? var_copy_map[var1] : var1
        expr_var2 = var2 in comp_vars ? var_copy_map[var2] : var2
        
        add_to_expression!(new_expr, expr_var1 * expr_var2 * con_obj.func.terms[pair])
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
    con_obj::ScalarConstraint{GenericQuadExpr{Float64, V}, MOI.GreaterThan{Float64}},#ConstraintRef{Plasmo.OptiEdge{Plasmo.OptiGraph}, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}}},
    comp_vars,
    var_copy_map,
    next_object::T,
    next_node::N,
    link::Bool;
    slack::Bool = false,
    slack_up = nothing,
    slack_down = nothing
) where {V <: JuMP.AbstractVariableRef, T <: Plasmo.AbstractOptiGraph, N <: Union{Plasmo.OptiNode, Plasmo.RemoteNodeRef}}
    # Create empty expression (for a constraint)
    new_expr = GenericQuadExpr{Float64, V}()

    # Add variables or variable copy to the expression
    for (i, var) in enumerate(con_obj.func.aff.terms.keys)
        if var in comp_vars
            add_to_expression!(new_expr, var_copy_map[var] * con_obj.func.aff.terms[var])
        else
            add_to_expression!(new_expr, var * con_obj.func.aff.terms[var])
        end
    end

    for (i, pair) in enumerate(con_obj.func.terms.keys)
        var1 = pair.a
        var2 = pair.b

        expr_var1 = var1 in comp_vars ? var_copy_map[var1] : var1
        expr_var2 = var2 in comp_vars ? var_copy_map[var2] : var2
        
        add_to_expression!(new_expr, expr_var1 * expr_var2 * con_obj.func.terms[pair])
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


########## END QUADRATICS




# Define function to add theta to the objective
function _add_theta_to_objective!(
    obj_func::AffExpr,
    node::N,
    _theta::Vector{V}
) where {V <: JuMP.AbstractVariableRef, N <: Union{Plasmo.OptiNode, Plasmo.RemoteNodeRef}}
    for i in 1:length(_theta)
        add_to_expression!(obj_func, _theta[i] * 1)
    end
    JuMP.set_objective_function(node, obj_func)

end

function _add_theta_to_objective!(
    obj_func::V,
    node::N,
    _theta::Vector{V}
) where {V <: JuMP.AbstractVariableRef, N <: Union{Plasmo.OptiNode, Plasmo.RemoteNodeRef}}
    obj_func = obj_func + sum(_theta[i] for i in 1:length(_theta))
    JuMP.set_objective_function(node, obj_func)
end

function _set_theta_lower_bound!(object::T, M) where {T <: Plasmo.AbstractOptiGraph}
    theta = object[:_theta_node][:_theta]
    for i in 1:length(theta)
        JuMP.set_lower_bound(object[:_theta_node][:_theta][i], M)
    end
end

function _build_comp_vars_copies(node::Plasmo.OptiNode, num_var_copies::Int)
    return @variable(node, _comp_vars_copy[1:num_var_copies])
end

function _build_comp_vars_copies(rnode::Plasmo.RemoteNodeRef, num_var_copies::Int)
    rgraph = rnode.remote_graph

    darray = rgraph.graph
    pnode = Plasmo._convert_remote_to_proxy(rgraph, rnode)
    f = @spawnat rgraph.worker begin
        lgraph = Plasmo.local_graph(darray)
        lnode = Plasmo._convert_proxy_to_local(lgraph, pnode)
        vars = @variable(lnode, _comp_vars_copy[1:num_var_copies])
        [(var.index, Symbol(name(var))) for var in vars]
    end
    var_tuples = fetch(f)

    vars = [Plasmo.RemoteVariableRef(rnode, t[1], t[2]) for t in var_tuples]
    return vars
end

function _get_con_obj_var_sets(con_obj::ScalarConstraint{GenericAffExpr{Float64, V}, S}, next_object_nodes, last_object_nodes) where {V, S}
    vars = con_obj.func.terms.keys

    next_object_link_vars = [var for var in vars if JuMP.owner_model(var) in next_object_nodes]
    last_object_link_vars = [var for var in vars if JuMP.owner_model(var) in last_object_nodes]

    return vars, next_object_link_vars, last_object_link_vars
end

function _get_con_obj_var_sets(con_obj::ScalarConstraint{GenericQuadExpr{Float64, V}, S}, next_object_nodes, last_object_nodes) where {V, S}
    affvars = con_obj.func.aff.terms.keys
    quadterms = collect(con_obj.func.terms.keys)
    quadtermsa = [term.a for term in quadterms]
    quadtermsb = [term.b for term in quadterms]
    quadvars = unique(vcat(quadtermsa, quadtermsb))

    vars = union(affvars, quadvars)
    
    next_object_link_vars = [var for var in vars if JuMP.owner_model(var) in next_object_nodes]
    last_object_link_vars = [var for var in vars if JuMP.owner_model(var) in last_object_nodes]

    return vars, next_object_link_vars, last_object_link_vars
end
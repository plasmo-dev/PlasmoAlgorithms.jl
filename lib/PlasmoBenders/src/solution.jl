function _add_cut_constraint!(
    optimizer::BendersAlgorithm{Plasmo.OptiGraph},
    last_object::Plasmo.OptiGraph,
    theta_expr::V,
    rhs_expr::GenericAffExpr{Float64, Plasmo.NodeVariableRef}
) where V <: Union{GenericAffExpr{Float64, Plasmo.NodeVariableRef}, NodeVariableRef}

    if length(rhs_expr.terms) != 0
        #if length(unique(JuMP.owner_model.(keys(rhs_expr.terms)))) =< 1 #TODO: FIx this!
        #    @constraint(JuMP.owner_model(theta_expr), theta_expr >= rhs_expr)
        #else
        @linkconstraint(last_object, theta_expr >= rhs_expr)
        #end
    else
        @constraint(JuMP.owner_model(theta_expr), theta_expr >= rhs_expr)
    end

    #MOI.Utilities.reset_optimizer(last_object)
end

"""
    _add_Benders_cuts!(optimizer::BendersAlgorithm)

Add Benders cuts to each nested problem; uses results from the backward pass and forward
pass to create cuts.
"""
function _add_Benders_cuts!(optimizer::BendersAlgorithm)
    # Loop through each object; compile information and add Benders cut
    for i in 1:(length(optimizer.solve_order))
        last_object = optimizer.solve_order[i]
        next_objects = optimizer.solve_order_dict[last_object]

        if length(next_objects) > 0
            agg_rhs_expr = GenericAffExpr{Float64, Plasmo.NodeVariableRef}()

            for (j, object) in enumerate(next_objects)
                rhs_expr = GenericAffExpr{Float64, Plasmo.NodeVariableRef}()
                # Complicating variables are on previous object
                comp_vars = optimizer.comp_vars[object]

                # Phi is the solution of the previous object
                phis = optimizer.phis[object]
                next_phi = phis[length(phis)]

                # Values of the complicating variables at last iteration
                primal_iters = optimizer.primal_iters[object]
                last_primals = primal_iters[:, size(primal_iters, 2)]

                # Dual variables come from the linking of complicating variables to the next object
                dual_iters = optimizer.dual_iters[object]
                next_duals = dual_iters[:, size(dual_iters, 2)]
                add_to_expression!(rhs_expr, next_phi)
                for k in 1:length(comp_vars)
                    add_to_expression!(rhs_expr, next_duals[k] * (comp_vars[k] - last_primals[k]))
                end

                if get_multicut(optimizer)
                    theta_var = _get_theta(optimizer, last_object, j)
                    _add_cut_constraint!(optimizer, last_object, theta_var, rhs_expr)
                else
                    add_to_expression!(agg_rhs_expr, rhs_expr)
                end
            end

            if !(get_multicut(optimizer))
                theta_vars = _get_theta(optimizer, last_object)
                theta_expr = sum(theta_vars[k] for k in 1:length(theta_vars))
                _add_cut_constraint!(optimizer, last_object, theta_expr, rhs_expr)
            end
        end
    end
end

function _solve_for_strengthened_cuts(optimizer, i)
    next_object = optimizer.solve_order[i + 1]

    comp_vars = optimizer.comp_vars[next_object]
    var_copy_map = optimizer.var_copy_map[next_object]

    # Loop through each complicating variable and set them
    # Ensure that each complicating variable copy matches the original variable
    for (j, var) in enumerate(comp_vars)
        var_copy = var_copy_map[var]

        if JuMP.is_binary(var)
            JuMP.set_binary(var_copy)
        elseif JuMP.is_integer(var)
            JuMP.set_integer(var_copy)
        end

        if JuMP.has_lower_bound(var)
            JuMP.set_lower_bound(var_copy, JuMP.lower_bound(var))
        end
        if JuMP.has_upper_bound(var)
            JuMP.set_upper_bound(var_copy, JuMP.upper_bound(var))
        end
    end

    _update_objective_and_optimize(optimizer, next_object)

    # Unset binary/integer complicating variable copies
    # These get fixed, so they do not need to be integer
    for (j, var) in enumerate(comp_vars)
        var_copy = var_copy_map[var]

        if JuMP.is_binary(var)
            JuMP.unset_binary(var_copy)
        elseif JuMP.is_integer(var)
            JuMP.unset_integer(var_copy)
        end
    end
end

function _update_objective_and_optimize(optimizer, next_object)

    comp_vars = optimizer.comp_vars[next_object]
    var_copy_map = optimizer.var_copy_map[next_object]

    primal_iters = optimizer.primal_iters[next_object]
    last_primals = primal_iters[:, size(primal_iters, 2)]

    dual_iters = optimizer.dual_iters[next_object]
    next_duals = dual_iters[:, size(dual_iters, 2)]

    # Get the original objective function
    next_object_objective_function = JuMP.objective_function(next_object)

    # Ensure the objective function is an Expr, not a single variable
    if typeof(next_object_objective_function) == NodeVariableRef
        next_object_objective_function = AffExpr(0, next_object_objective_function => 1)
    end

    # Get the original constant
    original_constant = next_object_objective_function.constant

    # Add the lagrangean relaxation term to the objective
    for (j, var) in enumerate(comp_vars)
        add_to_expression!(next_object_objective_function, - next_duals[j] * var_copy_map[var])
        add_to_expression!(next_object_objective_function, next_duals[j] * last_primals[j])
    end

    # Set the new objective function (with lagrangean relaxation term)
    JuMP.set_objective_function(next_object, next_object_objective_function)
    #@objective(next_object, Min, next_object_objective_function)

    # Get the new solution
    optimize!(next_object)

    # Save the data from the new solution
    new_phi_LR = objective_value(next_object)
    push!(optimizer.phis_LR[next_object], new_phi_LR)

    # Reset objective
    for (j, var) in enumerate(comp_vars)
        add_to_expression!(next_object_objective_function, next_duals[j] * var_copy_map[var])
        add_to_expression!(next_object_objective_function, - next_duals[j] * last_primals[j])
    end

    JuMP.set_objective_function(next_object, next_object_objective_function)
    #@objective(next_object, Min, next_object_objective_function)

end

function _optimize_in_forward_pass!(optimizer, i, ub)
    next_object = optimizer.solve_order[i]

    comp_vars = optimizer.comp_vars[next_object]
    var_copy_map = optimizer.var_copy_map[next_object]
    primal_iters = optimizer.primal_iters[next_object]
    last_primals = primal_iters[:, size(primal_iters, 2)]

    # Fix primal solutions
    for (j, var) in enumerate(comp_vars)
        JuMP.fix(var_copy_map[var], last_primals[j], force = true)
    end

    # Optimize the next node
    optimize!(next_object)

    # Check termination status
    _check_termination_status(optimizer, next_object, i)

    # Add to the upper bound; if it's not the last object, subtract the cost-to-go from upper bound
    obj_val = JuMP.value(next_object, JuMP.objective_function(next_object))

    next_objects = optimizer.solve_order_dict[next_object]
    if length(next_objects) > 0
        theta_val = _theta_value(optimizer, next_object)
    else
        theta_val = 0
    end

    if get_regularize(optimizer)
        get_regularize_lbs(optimizer)[next_object] = obj_val

        _regularize_pass!(optimizer, next_object, ub)
    else
        # Save primal information to upcoming objects
        _add_to_upper_bound!(optimizer, next_object, ub)
        next_objects = optimizer.solve_order_dict[next_object]
        for object in next_objects
            next_comp_vars = optimizer.comp_vars[object]
            last_primals = [JuMP.value(next_object, var) for var in next_comp_vars]

            primal_iters = optimizer.primal_iters[object]
            optimizer.primal_iters[object] = hcat(primal_iters, last_primals)
        end

        if !optimizer.is_MIP
            next_duals = _get_next_duals(optimizer, next_object)
            dual_iters = optimizer.dual_iters[next_object]
            optimizer.dual_iters[next_object] = hcat(dual_iters, next_duals)

            next_phi = JuMP.objective_value(next_object)
            push!(optimizer.phis[next_object], next_phi)
        end

        # Save the solutions
        next_object_vars = JuMP.all_variables(next_object)
        optimizer.last_solutions[next_object] = [JuMP.value(next_object, var) for var in next_object_vars]
    end

    _check_fixed_slacks!(optimizer, next_object)

    for (j, var) in enumerate(comp_vars)
        JuMP.unfix(var_copy_map[var])
    end
end

function _optimize_in_backward_pass(optimizer, i)
    object = optimizer.solve_order[i]

    ### Unset binary/integer variables
    bin_vars_dict = optimizer.binary_map[object]
    int_vars_dict = optimizer.integer_map[object]

    bin_vars = collect(keys(bin_vars_dict))
    int_vars = collect(keys(int_vars_dict))

    # Unset binary/integer variables and set bounds
    JuMP.unset_binary.(bin_vars)
    bin_vars_not_fixed = (!).(JuMP.is_fixed.(bin_vars))
    JuMP.set_upper_bound.(bin_vars[bin_vars_not_fixed], 1)
    JuMP.set_lower_bound.(bin_vars[bin_vars_not_fixed], 0)
    JuMP.unset_integer.(int_vars)

    if i != 1
        comp_vars = optimizer.comp_vars[object]
        var_copy_map = optimizer.var_copy_map[object]
        primal_iters = optimizer.primal_iters[object]
        last_primals = primal_iters[:, size(primal_iters, 2)]

        # Fix primal solutions
        for (j, var) in enumerate(comp_vars)
            JuMP.fix(var_copy_map[var], last_primals[j], force = true)
        end
    end

    # Optimize the next node
    optimize!(object)

    # Check termination status
    _check_termination_status(optimizer, object, i)

    # If it's not the first node, save the dual and phi values
    if i != 1
        dual_iters = optimizer.dual_iters[object]
        next_duals = _get_next_duals(optimizer, object)

        next_phi = JuMP.value(object, JuMP.objective_function(object))

        optimizer.dual_iters[object] = hcat(dual_iters, next_duals)
        push!(optimizer.phis[object], next_phi)

        comp_vars = optimizer.comp_vars[object]
        var_copy_map = optimizer.var_copy_map[object]

        # Fix primal solutions
        for (j, var) in enumerate(comp_vars)
            JuMP.unfix(var_copy_map[var])
        end
    end

    _check_fixed_slacks!(optimizer, object)

    # Reset binary/integer variables
    bin_vars_with_lower_bound = JuMP.has_lower_bound.(bin_vars)
    bin_vars_with_upper_bound = JuMP.has_upper_bound.(bin_vars)
    JuMP.delete_lower_bound.(bin_vars[bin_vars_with_lower_bound])
    JuMP.delete_upper_bound.(bin_vars[bin_vars_with_upper_bound])
    JuMP.set_binary.(bin_vars)
    JuMP.set_integer.(int_vars) #TODO: fix bounds on integer vars?
end

"""
    _add_strengthened_cuts!(optimizer::BendersAlgorithm)

Add strengthened Benders cuts to each nested problem; Follows the process described by
Zou et al., 2019 (https://doi.org/10.1007/s10107-018-1249-5) where the cuts come from a
Lagrangian relaxation where the lagrange multipliers are the dual variables of the
backward pass.
"""
function _add_strengthened_cuts!(optimizer::BendersAlgorithm)
    # Loop through each object and add strengthened Benders cut
    Threads.@threads for i in 1:(length(optimizer.solve_order) - 1)
    #for i in 1:(length(optimizer.solve_order) - 1)
        _solve_for_strengthened_cuts(optimizer, i)
    end

    for i in 1:(length(optimizer.solve_order))
        last_object = optimizer.solve_order[i]
        next_objects = optimizer.solve_order_dict[last_object]

        if length(next_objects) > 0
            agg_rhs_expr = GenericAffExpr{Float64, Plasmo.NodeVariableRef}()
            agg_rhs_expr_LR = GenericAffExpr{Float64, Plasmo.NodeVariableRef}()

            for (j, object) in enumerate(next_objects)
                rhs_expr = GenericAffExpr{Float64, Plasmo.NodeVariableRef}()
                rhs_expr_LR = GenericAffExpr{Float64, Plasmo.NodeVariableRef}()

                # Complicating variables are on previous object
                comp_vars = optimizer.comp_vars[object]

                # Phi is the solution of the previous object
                phis = optimizer.phis[object]
                phis_LR = optimizer.phis_LR[object]

                next_phi = phis[length(phis)]
                next_phi_LR = phis_LR[length(phis_LR)]

                # Values of the complicating variables at last iteration
                primal_iters = optimizer.primal_iters[object]
                last_primals = primal_iters[:, size(primal_iters, 2)]

                # Dual variables come from the linking of complicating variables to the next object
                dual_iters = optimizer.dual_iters[object]
                next_duals = dual_iters[:, size(dual_iters, 2)]

                add_to_expression!(rhs_expr, next_phi)
                for k in 1:length(comp_vars)
                    add_to_expression!(rhs_expr, next_duals[k] * (comp_vars[k] - last_primals[k]))
                end

                add_to_expression!(rhs_expr_LR, next_phi_LR)
                for k in 1:length(comp_vars)
                    add_to_expression!(rhs_expr_LR, next_duals[k] * (comp_vars[k] - last_primals[k]))
                end
                if get_multicut(optimizer)
                    theta_var = _get_theta(optimizer, last_object, j)
                    if rhs_expr.constant > rhs_expr_LR.constant
                        _add_cut_constraint!(optimizer, last_object, theta_var, rhs_expr)
                    else
                        _add_cut_constraint!(optimizer, last_object, theta_var, rhs_expr_LR)
                    end
                else
                    add_to_expression!(agg_rhs_expr, rhs_expr)
                    add_to_expression!(agg_rhs_expr_LR, rhs_expr_LR)
                end

            end

            if !(get_multicut(optimizer))
                theta_vars = _get_theta(optimizer, last_object)
                theta_expr = sum(theta_vars[k] for k in 1:length(theta_vars))
                if agg_rhs_expr.constant > agg_rhs_expr_LR.constant
                    _add_cut_constraint!(optimizer, last_object, theta_expr, agg_rhs_expr)
                else
                    _add_cut_constraint!(optimizer, last_object, theta_expr, agg_rhs_expr_LR)
                end
            end
        end
    end
end

"""
    _add_initial_relaxed_cuts!(optimizer::BendersAlgorithm)
"""
function _add_initial_relaxed_cuts!(
    optimizer::BendersAlgorithm{T}
) where {T <: Union{Plasmo.OptiNode, Plasmo.OptiGraph}}
    link_vars_mapping = Dict{T, Vector{ConstraintRef}}()

    for i in 1:(length(optimizer.solve_order))
        last_object = optimizer.solve_order[i]

        # Unset binary variables
        bin_value_dict = optimizer.binary_map[last_object]
        int_value_dict = optimizer.integer_map[last_object]

        bin_vars = collect(keys(bin_value_dict))
        int_vars = collect(keys(int_value_dict))

        JuMP.unset_binary.(bin_vars)
        JuMP.unset_integer.(int_vars)
        bin_vars_not_fixed = (!).(JuMP.is_fixed.(bin_vars))
        JuMP.set_upper_bound.(bin_vars[bin_vars_not_fixed], 1)
        JuMP.set_lower_bound.(bin_vars[bin_vars_not_fixed], 0)

        # Add a constraint for linking variables on the linking constraint
        # These are needed for creating the Benders cuts
        next_objects = optimizer.solve_order_dict[last_object]
        for object in next_objects
            comp_vars = optimizer.comp_vars[object]
            var_copy_map = optimizer.var_copy_map[object]

            link_cons = @linkconstraint(optimizer.graph, [i = 1:length(comp_vars)], comp_vars[i] == var_copy_map[comp_vars[i]])
            link_vars_mapping[object] = link_cons
        end
    end

    optimize!(optimizer.graph)

    # Query dual, primal, and objective values, then add the Benders cut
    for i in 1:(length(optimizer.solve_order))
        last_object = optimizer.solve_order[i]
        next_objects = optimizer.solve_order_dict[last_object]

        agg_rhs_expr = GenericAffExpr{Float64, Plasmo.NodeVariableRef}()
        for (j, object) in enumerate(next_objects)
            rhs_expr = GenericAffExpr{Float64, Plasmo.NodeVariableRef}()

            link_cons = link_vars_mapping[object]
            comp_vars = optimizer.comp_vars[object]

            next_duals = JuMP.dual.(link_cons)
            last_primals = [JuMP.value(object, comp_var) for comp_var in comp_vars]
            next_phi = JuMP.value(object, JuMP.objective_function(object))

            dual_iters = optimizer.dual_iters[object]
            primal_iters = optimizer.primal_iters[object]

            push!(optimizer.phis[object], next_phi)
            optimizer.dual_iters[object] = hcat(dual_iters, next_duals)
            optimizer.primal_iters[object] = hcat(primal_iters, last_primals)

            add_to_expression!(rhs_expr, next_phi)
            for k in 1:length(comp_vars)
                add_to_expression!(rhs_expr, next_duals[k] * (comp_vars[k] - last_primals[k]))
            end

            if get_multicut(optimizer)
                theta_var = _get_theta(optimizer, last_object, j)
                _add_cut_constraint!(optimizer, last_object, theta_var, rhs_expr)
            else
                add_to_expression!(agg_rhs_expr, rhs_expr)
            end

        end

        # Add lower bound on theta as defined by the user
        if length(next_objects) > 0 && !(get_multicut(optimizer))
            theta_vars = _get_theta(optimizer, last_object)
            theta_expr = sum(theta_vars[k] for k in 1:length(theta_vars))
            _set_theta_lower_bound!(last_object, optimizer.M)
            _add_cut_constraint!(optimizer, last_object, theta_expr, agg_rhs_expr)
        end
    end

    # Reset variables as binary and integer
    for i in 1:(length(optimizer.solve_order))
        current_object = optimizer.solve_order[i]

        bin_value_dict = optimizer.binary_map[current_object]
        int_value_dict = optimizer.integer_map[current_object]

        bin_vars = collect(keys(bin_value_dict))
        int_vars = collect(keys(int_value_dict))

        bin_vars_with_lower_bound = JuMP.has_lower_bound.(bin_vars)
        bin_vars_with_upper_bound = JuMP.has_upper_bound.(bin_vars)
        JuMP.delete_lower_bound.(bin_vars[bin_vars_with_lower_bound])
        JuMP.delete_upper_bound.(bin_vars[bin_vars_with_upper_bound])
        JuMP.set_binary.(bin_vars)
        JuMP.set_integer.(int_vars)
    end

    optimizer.ext["link_var_mapping"] = link_vars_mapping
end

function _add_to_upper_bound!(optimizer::BendersAlgorithm, object::OptiGraph, ub)
    obj_val = JuMP.value(object, JuMP.objective_function(object))
    ub[1] += obj_val
    if length(optimizer.solve_order_dict[object]) > 0
        theta_val = _theta_value(optimizer, object)
        ub[1] -= theta_val
    end
end

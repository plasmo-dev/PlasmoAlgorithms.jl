_rd_fields = [
    :objective_function,
    :constraint,
    :slack,
    :ubs,
    :lbs,
    :best_ub
]

for field in _rd_fields
    method = Symbol("get_regularize_", field)
    @eval begin
        @doc """
            $($method)(optimizer::BendersAlgorithm)
        Return the value of $($(QuoteNode(field))) from the `regularize_data` field of the `BendersAlgorithm`
        """
        $method(optimizer::BendersAlgorithm) = getproperty(optimizer.regularize_data, $(QuoteNode(field)))
    end
    @eval export $method
end


function _construct_regularize!(
    optimizer::BendersAlgorithm{T}
) where {T <: Union{Plasmo.OptiGraph, Plasmo.OptiNode}}
    for (i, object) in enumerate(optimizer.solve_order)
        next_objects = optimizer.solve_order_dict[object]
        if length(next_objects) > 0
            object = optimizer.solve_order[i]

            obj_dict = get_regularize_objective_function(optimizer)
            obj_dict[object] = objective_function(object)
        end

        get_regularize_best_ub(optimizer)[object] = Inf
        get_regularize_ubs(optimizer)[object] = Inf
        get_regularize_lbs(optimizer)[object] = -Inf
    end
end

function _regularize_pass!(
    optimizer::BendersAlgorithm{T},
    object,
    ub
) where {T <: Union{Plasmo.OptiGraph, Plasmo.OptiNode}}
    next_objects = optimizer.solve_order_dict[object]
    if length(next_objects) > 0
        original_objective = get_regularize_objective_function(optimizer)[object]

        alpha = get_regularize_param(optimizer)
        LB = get_regularize_lbs(optimizer)[object]
        UB = get_regularize_best_ub(optimizer)[object]

        rhs =  maximum([LB + alpha * (UB - LB), LB])

        if rhs != Inf
            @linkconstraint(object, _reg_con, original_objective <= rhs)
        end

        @objective(object, Min, 0 * all_variables(object)[1])

        optimize!(object)
        if termination_status(object) != MOI.INFEASIBLE
            ub[1] += value(object, original_objective) - _theta_value(optimizer, object)
            get_regularize_ubs(optimizer)[object] = value(object, original_objective) - _theta_value(optimizer, object)
            for next_object in next_objects
                comp_vars = optimizer.comp_vars[next_object]
                last_primals = JuMP.value.(comp_vars)

                old_primals = optimizer.primal_iters[next_object]
                optimizer.primal_iters[next_object] = hcat(old_primals, last_primals)
            end
        else
            error(
                "Regularization failed with infeasible regularization subproblem." *
                "Regularization problem should be feasible; please turn off regularization and open an issue"
            )
        end

        if !optimizer.is_MIP
            if haskey(optimizer.parent_objects, object)
                next_duals = _get_next_duals(optimizer, object)
                dual_iters = optimizer.dual_iters[object]
                optimizer.dual_iters[object] = hcat(dual_iters, next_duals)

                next_phi = JuMP.objective_value(object)
                push!(optimizer.phis[object], next_phi)
            end
        end

        # Save the solutions
        object_vars = JuMP.all_variables(object)
        optimizer.last_solutions[object] = JuMP.value.(object_vars)

        @objective(object, Min, original_objective)

        if haskey(object, :_reg_con)
            MOI.delete(JuMP.owner_model(object[:_reg_con]), object[:_reg_con])
            delete!(object_dictionary(object), :_reg_con)
        end
    else
        _add_to_upper_bound!(optimizer, object, ub)
        get_regularize_ubs(optimizer)[object] = value(object, objective_function(object))

        if !optimizer.is_MIP
            next_duals = _get_next_duals(optimizer, object)
            dual_iters = optimizer.dual_iters[object]
            optimizer.dual_iters[object] = hcat(dual_iters, next_duals)

            next_phi = JuMP.objective_value(object)
            push!(optimizer.phis[object], next_phi)
        end

        # Save the solutions
        object_vars = JuMP.all_variables(object)
        optimizer.last_solutions[object] = JuMP.value.(object_vars)
    end
end

function _init_regularize_bounds!(optimizer::BendersAlgorithm)
    len_solve_order = length(optimizer.solve_order)
    for i in len_solve_order:-1:1
        object = optimizer.solve_order[i]

        new_ub = get_regularize_ubs(optimizer)[object]
        if new_ub < get_regularize_best_ub(optimizer)[object]
            get_regularize_best_ub(optimizer)[object] = new_ub
        end

        if i != 1
            parent_object = optimizer.parent_objects[object]
            get_regularize_ubs(optimizer)[parent_object] += get_regularize_ubs(optimizer)[object]
        end
    end
end

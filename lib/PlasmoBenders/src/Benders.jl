"""
    AbstractPBAlgorithm{T}
"""
abstract type AbstractPBAlgorithm{T, V} end

abstract type AbstractPBOptions end

abstract type PBData{T} end

"""
    BendersOptions

Options object for dual dynamic programming with graphs

Attributes include
 - `strengthened::Bool` - whether to use strengthened cuts
 - `multicut::Bool` - whether the problem should use aggregated or multiple cuts when
   there are multiple objects generating information in the forward pass
 - `feasibility_cuts::Bool` - whether to allow feasibility cuts if infeasible
 - `regularize::Bool` - whether to used regularization for getting next cuts
 - `parallelize_benders::Bool` - whether to parallelize Benders problems if applicable
 - `parallelize_forward::Bool` - whether to parallelize the forward pass if possible
 - `parallelize_backward::Bool` - whether to parallelize the backward pass if possible
 - `add_slacks::Bool` - whether to add slack variables to linking constraints
 - `fix_slacks::Bool` - whether to fix the slack variables to zero; slacks will be relaxed
   if the problem is infeasible
 - `warm_start::Bool` - whether to warm start the problem using the previous best solution
 - `relaxed_init_cuts::Bool` - whether to create some initial cuts by relaxing the problem
   and solving the full problem as an LP; only applies if the problem is a MIP
 - `slack_penalty::Real` - penalty for nonzero slacks; only applies for `add_slacks = true`
 - `regularize_param::Real` - parameter for regularization; essentially how far from optimal
   the solution can end up being for choosing the next iterates.
"""
mutable struct BendersOptions <: AbstractPBOptions
    strengthened::Bool
    multicut::Bool
    feasibility_cuts::Bool
    regularize::Bool
    parallelize_benders::Bool
    parallelize_forward::Bool
    parallelize_backward::Bool
    add_slacks::Bool
    fix_slacks::Bool
    warm_start::Bool
    relaxed_init_cuts::Bool

    slack_penalty::Real
    regularize_param::Real

    function BendersOptions()
        options = new()

        options.strengthened = false
        options.multicut = false
        options.feasibility_cuts = false
        options.regularize = false
        options.parallelize_benders = false
        options.parallelize_forward = false
        options.parallelize_backward = false
        options.add_slacks = false
        options.fix_slacks = false
        options.warm_start = true
        options.relaxed_init_cuts = false

        options.slack_penalty = 1.0e6
        options.regularize_param = 0.5

        return options
    end
end

"""
    RegularizeData

Data structure for storing regularization data

All attributes include a dictionary which is indexed by the subproblem objects of Benders
 - `objective_function` - the objective function of an object
 - `ubs` - value of the object's objective (without theta)
 - `lbs` - value of the object's objective (with theta)
 - `best_ub` - best upper bound
"""
mutable struct RegularizeData{T} <: PBData{T}
    objective_function::Dict
    ubs::Dict
    lbs::Dict
    best_ub::Dict

    function RegularizeData{T}() where {T <: Plasmo.AbstractOptiGraph}
        rd = new()

        rd.objective_function = Dict{T, Any}()
        rd.ubs = Dict{T, Float64}()
        rd.lbs = Dict{T, Float64}()
        rd.best_ub = Dict{T, Float64}()

        return rd
    end
end

"""
    BendersAlgorithm

Optimizer object for dual dynamic programming with graphs. Currently only implemented for linear tree structures.

Attributes include the following (where noted as dictionaries, these are mappings from the
Benders subproblems to the data described)
 - `graph` - Plasmo OptiGraph for appplying Benders
 - `root_object` - node or subgraph in `graph` indicating where to start the Benders algorithm
 - `is_MIP` = Boolean indicating if the problem is a MIP or not
 - `solve_order` - vector of OptiNodes in the order that they are solved by Benders
 - `solve_order_dict` - dictionary mapping to a vector of all the "next objects" for a given object
 - `parent_objects` - dictionary mapping to the previous subproblem
 - `max_iters` - maximum number of Benders iterations
 - `time_limit` - maximum time (seconds) allowed for running Benders
 - `tol` - (absolute) termination tolerance between upper and lower bounds
 - `current_iter` - current iteration of Benders algorithm
 - `M` - lower bound for the cost-to-go function at each iteration
 - `dual_iters` - dictionary mapping optinodes to the corresponding dual values
    at each iteration; dual values come from solution of following node
 - `primal_iters` - dictionary mapping optinodes to the primal solutions of the
    complicating variables on the given node
 - `phis` - dictionary mapping optinodes to the objective of the immediately following node
 - `phis_LR` - dictionary mapping optinodes to the objective of the lagrangean (used for strenghtened cuts)
    relaxation of the following node; used for strengthened Benders cuts in MIPs
 - `time_forward_pass` - time spent in the forward pass (seconds)
 - `time_backward_pass` - time spent in the backward pass (seconds)
 - `time_init` - time initializing the optimizer (seconds)
 - `time_iterations` - vector of the times spent in each Iteration
 - `comp_vars` - dictionary mapping the node to a vector of its complicating variables
 - `comp_var_map` - dictionary mapping the node to a dictionary of complicating variables
    mapped to their index in `comp_vars`
 - `var_copy_map` - dictionary mapping the node to a dictionary mapping the complicating
    variables to their copies on the following node
 - `objective_value` - the final objective value (from upper bound)
 - `lower_bounds` - vector of lower bounds at each iteration
 - `upper_bounds` - vector of upper bounds at each iteration
 - `regularize_data` - regularization data object
 - `binary_map` - dictionary mapping the nodes to a vector of binary variables on the node
 - `integer_map` - dictionary mapping the nodes to a vector of integer variables on the node
 - `last_solutions` - dictionary mapping the nodes to a vector of the last solutions
 - `var_solution_map` - dictionary mapping the variables to their index on the `last_solutions` vector
 - `var_to_graph_map` - dictionary mapping the variables to their owning subproblem subgraph
 - `feasibility_map` - dictionary matching `solve_order` for whether a problem was feasible or not
 - `options` - solver options for Benders algorithm
 - `ext` - Dictionary for extending certain procedures
"""
mutable struct BendersAlgorithm{T, V} <: AbstractPBAlgorithm{T, V}
    graph::T
    root_object::T

    is_MIP::Bool
    solve_order::Vector{T}
    solve_order_dict::Dict{T, Vector{T}}
    parent_objects::Dict{T, T}

    max_iters::Int
    time_limit::Real
    tol::Real
    current_iter::Int
    M::Real

    dual_iters::Dict{T, Matrix{Float64}}
    primal_iters::Dict{T, Matrix{Float64}}
    phis::Dict{T, Vector{Float64}}
    phis_LR::Dict{T, Vector{Float64}}

    time_forward_pass::Float64
    time_backward_pass::Float64
    time_init::Float64
    time_iterations::Vector{Float64}

    comp_vars::Dict{T, Vector{V}}
    comp_var_map::Dict{T, Dict{V, Int}}
    var_copy_map::Dict{T, Dict{V, V}}
    slack_vars::Dict{T, Vector{V}}

    objective_value::Float64
    status::MOI.TerminationStatusCode
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}

    regularize_data::RegularizeData{T}

    binary_map::Dict{T, Dict{V, Float64}}
    integer_map::Dict{T, Dict{V, Float64}}

    last_solutions::Dict{T, Vector{Float64}}
    var_solution_map::Dict{T, Dict{V, Int}}
    var_to_graph_map::Dict{V, T}
    best_solutions::Dict{T, Vector{Float64}}
    best_upper_bound::Float64

    feasibility_map::Dict{T, Bool}
    options::BendersOptions

    ext::Dict

    function BendersAlgorithm{T, V}() where {T <: Plasmo.AbstractOptiGraph, V <: JuMP.AbstractVariableRef}
        function build_graph(T)
            if T == Plasmo.OptiGraph
                return Plasmo.OptiGraph()
            elseif T == Plasmo.RemoteOptiGraph
                return Plasmo.RemoteOptiGraph()
            else
                error("$T must be a `Plasmo.OptiGraph` of `Plasmo.RemoteOptiGraph`")
            end
        end
        optimizer = new{T, V}()
        g = build_graph(T)
        optimizer.graph = build_graph(T)
        optimizer.is_MIP = false

        optimizer.root_object = T()
        optimizer.solve_order = Vector{T}()
        optimizer.solve_order_dict = Dict{T, Vector{T}}()
        optimizer.parent_objects = Dict{T, T}()

        optimizer.max_iters = 0
        optimizer.time_limit = Inf
        optimizer.tol = 0.
        optimizer.current_iter = 0
        optimizer.M = 0.

        optimizer.dual_iters = Dict{T, Matrix{Float64}}()
        optimizer.primal_iters = Dict{T, Matrix{Float64}}()
        optimizer.phis = Dict{T, Vector{Float64}}()
        optimizer.phis_LR = Dict{T, Vector{Float64}}()

        optimizer.time_forward_pass = 0.
        optimizer.time_backward_pass = 0.
        optimizer.time_init = 0.
        optimizer.time_iterations = Vector{Float64}()

        optimizer.comp_vars = Dict{T, Vector{V}}()
        optimizer.comp_var_map = Dict{T, Dict{V, Int}}()
        optimizer.var_copy_map = Dict{T, Dict{V, V}}()
        optimizer.slack_vars = Dict{T, Vector{V}}()

        optimizer.objective_value = Inf
        optimizer.status = MOI.OPTIMIZE_NOT_CALLED
        optimizer.lower_bounds = Vector{Float64}()
        optimizer.upper_bounds = Vector{Float64}()

        optimizer.regularize_data = RegularizeData{T}()

        optimizer.binary_map = Dict{T, Dict{V, Float64}}()
        optimizer.integer_map = Dict{T, Dict{V, Float64}}()

        optimizer.last_solutions = Dict{T, Vector{Float64}}()
        optimizer.var_solution_map = Dict{T, Dict{V, Int}}()
        optimizer.var_to_graph_map = Dict{V, Int}()
        optimizer.best_solutions = Dict{T, Vector{Float64}}()
        optimizer.best_upper_bound = Inf

        optimizer.feasibility_map = Dict{T, Bool}()
        optimizer.options = BendersOptions()

        optimizer.ext = Dict{String, Any}()

        return optimizer
    end
end

BendersAlgorithm() = BendersAlgorithm{Plasmo.OptiGraph, NodeVariableRef}()

"""
    BendersAlgorithm(graph, root_object; kwargs...)

Function for creating the BendersAlgorithm object from `graph`. `root_object` must be an
OptiNode or subgraph on `graph`. key ward arguments include the following
 - `max_iters = 100` - maximum number of iterations
 - `tol = 1e-7` - termination tolerance between upper and lower bounds
 - `M = 0` - lower bound on cost-to-go estimator for each node
 - `is_MIP = nothing` - indicates if the problem is a MIP. If it is passed as `nothing`,
    PlasmoBenders will check to see if it is a MIP, and will set `is_MIP` accordingly
 - `solver = nothing` - if defined, this solver will be set for all subproblems
 - `strengthened = false` - whether to use strengthened Benders cuts (see https://doi.org/10.1007/s10107-018-1249-5.)
 - `multicut = true` - whether to use multicuts (rather than aggregated cuts) when applicable
 - `feasibility_cuts = false` - whether to allow feasibility cuts
 - `regularize = false` - whether to regularize solution of next iterates
 - `parallelize_benders = false` - whether to parallelize subproblem solution when the problem
   has a Benders-type structure defined
 - `parallelize_forward = false` - whether to parallelize forward pass if possible; not yet supported
 - `parallelize_backward = false` - whether to parallelize backward pass
 - `add_slacks::Bool = false` - whether to add slack variables to the linking constraints to
    help ensure feasibility between solutions; slack variables are penalized in objective
 - `fix_slacks = false` - whether to fix the slack variables to zero they only relax if a problem
   is infeasible
 - `warm_start::Bool = true` - whether to set the previous iterations solutions as the
    starting values for the next iterations forward pass
 - `relaxed_init_cuts = false` - whether to generate initial cuts by relaxing the problem and
   solving the whole problem; only applies for MILPs
 - `slack_penalty = 1e6` - coefficient on slack variables in objective
 - `regularize_param = 0.5` - regularization parameter; must be between 0 and 1
"""
function BendersAlgorithm(
    graph::Plasmo.OptiGraph,
    root_object::T;
    solver = nothing,
    args...
) where {T <: Plasmo.OptiNode}

    nodes = all_nodes(graph)
    if !(root_object in nodes)
        error("root_object is not in the nodes of the graph")
    end
    if isnothing(solver)
        @warn("You have set subproblems as nodes. These will be restructured as graphs,
        but no solver has been set. Make sure to set a solver on the resulting subgraphs")
    end
    node_membership_vector = [i for i in 1:length(nodes)]
    for i in 1:length(nodes)
        if root_object == nodes[i]
            node_membership_vector[i] = 1
            node_membership_vector[1] = i
            break
        end
    end
    node_partition = Partition(graph, node_membership_vector)

    apply_partition!(graph, node_partition)
    start_graph = local_subgraphs(graph)[1]

    return BendersAlgorithm(graph, start_graph; solver = solver, args...)
end

function BendersAlgorithm(
    graph::T,
    root_object::T;
    max_iters = 100,
    time_limit = Inf,
    tol::Float64 = 1e-7,
    M = 0.,
    is_MIP = nothing,
    solver = nothing,
    strengthened::Bool = false,
    multicut::Bool = true,
    feasibility_cuts::Bool = false,
    regularize::Bool = false,
    parallelize_benders::Bool = false,
    parallelize_forward::Bool = false,
    parallelize_backward::Bool = false,
    add_slacks::Bool = false,
    fix_slacks::Bool = false,
    warm_start::Bool = true,
    relaxed_init_cuts::Bool = false,
    slack_penalty = 1e6,
    regularize_param::Real = 0.5
) where {T <: Plasmo.AbstractOptiGraph}

    if !(root_object in local_subgraphs(graph))
        error("root_object is not defined in the graph")
    end

    println("Initializing BendersAlgorithm...")
    flush(stdout)

    time_init = @elapsed begin

        # Initiailize optimizer and graph
        V = Plasmo.variable_type(graph)
        optimizer = BendersAlgorithm{T, V}()
        optimizer.graph = graph

        set_strengthened!(optimizer, strengthened)
        set_multicut!(optimizer, multicut)
        set_feasibility_cuts!(optimizer, feasibility_cuts)
        set_regularize!(optimizer, regularize)
        set_parallelize_benders!(optimizer, parallelize_benders)
        set_parallelize_forward!(optimizer, parallelize_forward)
        set_parallelize_backward!(optimizer, parallelize_backward)
        set_add_slacks!(optimizer, add_slacks)
        set_fix_slacks!(optimizer, fix_slacks)
        set_warm_start!(optimizer, warm_start)
        set_relaxed_init_cuts!(optimizer, relaxed_init_cuts)

        set_slack_penalty!(optimizer, slack_penalty)
        set_regularize_param!(optimizer, regularize_param)

        if parallelize_forward
            @warn("`parallelize_forward` is not yet supported. Benders will run, but the forward pass will not be parallelized")
        end
        if feasibility_cuts
            @warn(
                "`feasibility_cuts` have been implemented but are still under development. " *
                "If you experience bugs, please report them in the package's github issue tracker"
            )
        end
        if get_regularize(optimizer)
            if regularize_param < 0 || regularize_param > 1
                error("regularize_param is $regularize_param but needs to be in [0, 1]")
            end
        end
        if get_regularize(optimizer) && get_warm_start(optimizer)
            @warn("`regularize` and `warm_start` cannot both be true; setting `warm_start` to false")
            set_warm_start!(optimizer, false)
        end
        if !(get_multicut(optimizer)) && get_feasibility_cuts(optimizer)
            @warn("`multicut` cannot be false while `feasibility_cuts`` is true; setting `multicut` to true")
            set_multicut!(optimizer, true)
        end

        # Set initial data
        optimizer.root_object = root_object
        optimizer.max_iters = max_iters
        optimizer.time_limit = time_limit
        optimizer.tol = tol
        optimizer.M = M
        optimizer.feasibility_map[root_object] = true

        # Set is_MIP
        if isnothing(is_MIP)
            _set_is_MIP(optimizer)
        else
            optimizer.is_MIP = is_MIP
        end

        # Set solver if it is defined
        if !isnothing(solver)
            set_optimizer(optimizer.graph, solver)
            objects = _get_objects(graph)
            for i in 1:length(objects)
                object = objects[i]
                set_optimizer(object, solver)
            end
        end

        # Add start object
        push!(optimizer.solve_order, root_object)

        _init_ext!(optimizer)

        # Add second object to solve order
        _add_second_object!(optimizer, get_relaxed_init_cuts(optimizer))

        while length(optimizer.ext["search_next"]) > 0
            search_next = optimizer.ext["search_next"][1]
            parent_object = optimizer.parent_objects[search_next]
            ############### Add complicating variables ##############
            # Get the linking constraints between last and current node
            _add_complicating_variables!(optimizer, parent_object, search_next, get_add_slacks(optimizer), get_slack_penalty(optimizer))
            optimizer.feasibility_map[search_next] = true

            ################ Get the next object(s) in sequene ###################
            _add_next_object!(optimizer, parent_object, search_next, get_relaxed_init_cuts(optimizer))
        end

        _check_two_stage_tree(optimizer)

        # Check to make sure all objects have been added to the solve_order vector
        if length(intersect(_get_objects(optimizer), optimizer.solve_order)) != length(_get_objects(optimizer))
            error("Number of nodes/graphs being solved does not match the number of nodes/graphs in the overall graph")
        end

        if get_add_slacks(optimizer) && get_fix_slacks(optimizer)
            slack_vars = optimizer.slack_vars
            for i in 2:length(optimizer.solve_order)
                object = optimizer.solve_order[i]

                slacks = slack_vars[object]

                for var in slacks
                    JuMP.fix(var, 0, force = true)
                end
            end
        end

        var_solution_map = optimizer.var_solution_map
        var_to_graph_map = optimizer.var_to_graph_map
        # Create mapping of binary and integer variables on each node
        for (i, object) in enumerate(optimizer.solve_order)
            object_var_solution_map = Dict{V, Int}()
            all_vars = all_variables(object)

            bin_value_dict = Dict{V, Float64}()
            int_value_dict = Dict{V, Float64}()

            binary_bool_vec = _get_binary_bool_vector(object)
            integer_bool_vec = _get_integer_bool_vector(object)
            for i in 1:length(binary_bool_vec)
                if binary_bool_vec[i]
                    bin_value_dict[all_vars[i]] = 0
                end
            end
            for i in 1:length(integer_bool_vec)
                if integer_bool_vec[i]
                    int_value_dict[all_vars[i]] = 0
                end
            end

            for (j, var) in enumerate(all_vars)
                object_var_solution_map[var] = j
                var_to_graph_map[var] = object
            end

            # TODO: Maybe save solutions on the remote workers or each individual subgraph
            # the challenge here is that we need to make sure we are saving remotely
            optimizer.var_solution_map[object] = object_var_solution_map
            optimizer.binary_map[object] = bin_value_dict
            optimizer.integer_map[object] = int_value_dict

            optimizer.last_solutions[object] = zeros(Float64, length(all_vars))
        end

        Plasmo.set_to_node_objectives(optimizer.graph)
        if get_relaxed_init_cuts(optimizer) && optimizer.is_MIP
            _add_initial_relaxed_cuts!(optimizer)
        end

        #Plasmo._init_graph_backend(optimizer.graph)

        for i in 1:length(optimizer.solve_order)
            Plasmo.set_to_node_objectives(optimizer.solve_order[i])
        end
        Plasmo.set_to_node_objectives(optimizer.graph)
        # Construct regularization abilities if needed
        if get_regularize(optimizer)
            _construct_regularize!(optimizer)
        end
    end

    println("BendersAlgorithm Initialized!")

    optimizer.time_init += time_init

    return optimizer
end

"""
    run_algorithm!(optimizer::PB.BendersAlgorithm; output::Bool = true, run_gc::Bool = false)

Optimize the graph in BendersAlgorithm by using the Benders algorithm. Keyword argument `output`
indicates whether to print the upper/lower bounds and gap and each iteraiton. `run_gc` indicates
whether to run the garbage collector after each iteraiton. 
"""
function run_algorithm!(
    optimizer::BendersAlgorithm;
    output::Bool = true,
    run_gc::Bool = false
)

    println("################################################")
    println("Running BendersAlgorithm")
    println("################################################")
    println()
    println("Number of Variables: $(num_variables(optimizer.graph))")
    println("Number of Subproblems: $(length(optimizer.solve_order))")
    println()

    # Run the optimizer for up to max iterations
    while optimizer.current_iter < optimizer.max_iters
        if get_regularize(optimizer)
            _init_regularize_bounds!(optimizer)
        end
        # Get upper bound from forward pass
        time_forward_pass = @elapsed begin
            ub, lb = _forward_pass!(optimizer)
            push!(optimizer.upper_bounds, ub)

            # If the new solution is better than former solutions, save it
            if ub < optimizer.best_upper_bound || length(optimizer.upper_bounds) == 1
                optimizer.best_upper_bound = ub
                for i in 1:length(optimizer.solve_order)
                    object = optimizer.solve_order[i]
                    optimizer.best_solutions[object] = optimizer.last_solutions[object]
                end

                # Warm start by setting start values at the best solution
                if get_warm_start(optimizer)
                    for object in optimizer.solve_order
                        _warm_start(optimizer, object)
                    end
                end
            end
        end

        # Save forward pass time
        optimizer.time_forward_pass += time_forward_pass

        # Save the best upper bound
        best_upper_bound = optimizer.best_upper_bound

        # Save the lower bound
        push!(optimizer.lower_bounds, lb)

        # Get the lower bound from backward pass
        time_backward_pass = @elapsed begin
            if optimizer.is_MIP
                _backward_pass!(optimizer; strengthened = get_strengthened(optimizer))
            end
        end

        # Save backward pass time
        optimizer.time_backward_pass += time_backward_pass

        # Save iteration time
        push!(optimizer.time_iterations, time_forward_pass + time_backward_pass)

        # Define error
        if best_upper_bound != 0
            err = abs((best_upper_bound - lb) / best_upper_bound)
        else
            err = abs(best_upper_bound - lb)
        end


        # Update the iteration
        optimizer.current_iter += 1

        # Provide output
        if output
            if optimizer.current_iter % 20 == 0 || optimizer.current_iter == 1
                @printf "%4s | %8s | %8s | %8s | %8s" "Iter" "Gap" "LowerBound" "UpperBound" "Time (s)\n"
            end
            @printf(
                "%4i | %7.3g %% | %7.3e | %7.3e | %7.3e\n",
                optimizer.current_iter,
                round(err*100, digits = 3),
                optimizer.lower_bounds[end],
                optimizer.best_upper_bound,
                optimizer.time_iterations[end]
            )
            flush(stdout)
        end

        if run_gc
            GC.gc()
        end

        # End function if the error is within the tolerance
        if err < optimizer.tol
            optimizer.objective_value = optimizer.best_upper_bound
            println("Optimal Solution Found!")
            optimizer.status = MOI.OPTIMAL
            return nothing
        end

        if sum(optimizer.time_iterations) > optimizer.time_limit
            optimizer.status = MOI.TIME_LIMIT
            break
        end
    end

    # Upper/lower bounds
    final_lower_bound = optimizer.lower_bounds[end]
    final_upper_bound = optimizer.best_upper_bound

    # Get error
    err = abs(final_upper_bound - final_lower_bound)
    optimizer.objective_value = optimizer.best_upper_bound

    # Return maximum number of iterations reached
    if sum(optimizer.time_iterations) > optimizer.time_limit
        println("Maximum Time Limit of $(optimizer.time_limit) s Exceeded")
    else
        optimizer.status = MOI.ITERATION_LIMIT
        println("Maximum Number of Iterations Exceeded")
    end
    println("Optimality Gap at Termination was ", err)
    return nothing
end

"""
    _forward_pass!(optimizer::BendersAlgorithm)

Runs the forward pass for Benders. Follows the node order in `solve_order` attribute. Save the
primal information for complicating variables on each node. If it is not a MIP, also saves
the dual and objective value information and adds the Benders cut. Returns the upper and
lower bounds, where lower bound is only valid if the problem is not a MIP (otherwise the lower
bound comes from the backward pass)
"""
function _forward_pass!(optimizer::BendersAlgorithm)
    ########## Solve the first node ############
    root_object = optimizer.solve_order[1]

    # acs = all_constraints(root_object, include_variable_in_set_constraints = false)
    # if optimizer.current_iter >= 6
    #     offset = optimizer.current_iter - 6
    #     set_normalized_rhs(acs[77 + offset], -1e9)
    # end

    JuMP.optimize!(root_object)

    root_object_feasibility = _check_termination_status(optimizer, root_object, 1)

    # acs = all_constraints(root_object, include_variable_in_set_constraints = false)
    # lhs = [value(root_object, constraint_object(c).func) for c in acs[77:length(acs)]]
    # rhs = [constraint_object(c).set.lower for c in acs[77:length(acs)]]
    # comp_vals = [(lhs[i], rhs[i]) for i in 1:length(rhs)]
    # condition = x -> abs(x[1] - x[2]) < 1e-4
    # active_set_idxs = findall(condition, comp_vals)
    # println("active constraint refs: ", active_set_idxs, " out of $(length(acs)-76) constraints")

    # Get initial objective
    root_objective = JuMP.objective_value(root_object)

    # Initialize upper bound
    ub = [0.]
    lb = [0.]

    # Add objective (minus cost-to-go) to the upper bound value
    lb[1] += root_objective
    if get_regularize(optimizer)
        get_regularize_lbs(optimizer)[root_object] = lb[1]

        _regularize_pass!(optimizer, root_object, ub)
    else
        _add_to_upper_bound!(optimizer, root_object, ub)
        # Save primal information to upcoming objects
        next_objects = optimizer.solve_order_dict[root_object]
        for object in next_objects
            next_comp_vars = optimizer.comp_vars[object]
            last_primals = _get_variable_values(root_object, next_comp_vars)

            primal_iters = optimizer.primal_iters[object]
            optimizer.primal_iters[object] = hcat(primal_iters, last_primals)
        end

        # Save the solutions #FIX
        optimizer.last_solutions[root_object] = _get_object_last_solutions(root_object)
    end

    ############# Solve each successive object #################
    if get_parallelize_benders(optimizer)
        Threads.@threads for i in 2:(length(optimizer.solve_order))
            _optimize_in_forward_pass!(optimizer, i, ub)
        end
    else
        for i in 2:(length(optimizer.solve_order))
            _optimize_in_forward_pass!(optimizer, i, ub)
        end
    end

    ############### if it is not a MIP, do the backward pass now ###################
    if !optimizer.is_MIP
        _add_Benders_cuts!(optimizer)
    end

    return ub[1], lb[1]
end

"""
    _backward_pass!(optimizer::BendersAlgorithm)

Perform backward pass for a MIP problem. Backward pass can be done in parallel. Binary
and integer variables are relaxed and the problem is solved to get the dual and objective
values for producing Benders cuts. If `strengthened = true`, the stengthened cuts are ALSO
added using the approach of Zou et al. https://doi.org/10.1007/s10107-018-1249-5.
"""
function _backward_pass!(optimizer::BendersAlgorithm; strengthened::Bool = false)

    len_solve_order = length(optimizer.solve_order)

    # Perform backward pass in parallel
    #Threads.@threads for i in 1:len_solve_order
    if get_parallelize_backward(optimizer) || get_parallelize_benders(optimizer)
        Threads.@threads for i in 1:len_solve_order
            _optimize_in_backward_pass(optimizer, i)
        end
    else
        for i in 1:len_solve_order
            _optimize_in_backward_pass(optimizer, i)
        end
    end

    # Add constraint on the cost-to-go to last node Bender's cut
    if get_strengthened(optimizer)
        _add_strengthened_cuts!(optimizer)
    else
        _add_Benders_cuts!(optimizer)
    end

    return nothing
end

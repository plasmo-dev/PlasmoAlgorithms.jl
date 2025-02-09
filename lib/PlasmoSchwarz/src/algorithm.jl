"""
    Timers

Simple data structure to track the timing of various stages of the SchwarzAlgorithm.

Fields:
  - `start_time::Float64`: Time when the algorithm starts.
  - `initialize_time::Float64`: Time taken to initialize the algorithm.
  - `eval_objective_time::Float64`: Time taken to evaluate the objective function.
  - `eval_primal_feasibility_time::Float64`: Time taken to evaluate primal feasibility.
  - `eval_dual_feasibility_time::Float64`: Time taken to evaluate dual feasibility.
  - `communicate_time::Float64`: Time spent in inter-subproblem communication.
  - `update_subproblem_time::Float64`: Time taken to update subproblems.
  - `solve_subproblem_time::Float64`: Time taken to solve subproblems.
  - `total_time::Float64`: Total time for the entire process.
"""
@kwdef mutable struct Timers
    start_time::Float64 = 0.0
    initialize_time::Float64 = 0.0
    eval_objective_time::Float64 = 0.0
    eval_primal_feasibility_time::Float64 = 0.0
    eval_dual_feasibility_time::Float64 = 0.0
    communicate_time::Float64 = 0.0
    update_subproblem_time::Float64 = 0.0
    solve_subproblem_time::Float64 = 0.0
    total_time::Float64 = 0.0
end

"""
    Options

Configuration options for the algorithm.

Fields:
  - `tolerance::Float64`: Convergence tolerance for primal and dual errors (default: `1e-4`).
  - `max_iterations::Int64`: Maximum number of iterations to run (default: `1000`).
  - `mu::Float64`: Penalty parameter for the augmented Lagrangian (default: `1.0`).
  - `overlap_distance::Int64`: Distance for overlap in subproblem expansion (default: `1`).
  - `use_node_objectives::Bool`: Whether to use node-specific objectives instead of a global graph objective (default: `true`).
  - `subproblem_optimizer`: The optimizer to use for solving subproblems (default: `nothing`).
"""
@kwdef mutable struct Options
    tolerance::Float64 = 1e-4         # primal and dual tolerance measure
    max_iterations::Int64 = 1000      # maximum number of iterations
    mu::Float64 = 1.0                 # augmented lagrangian penalty
    overlap_distance::Int64 = 1
    use_node_objectives::Bool = true  # whether to ignore graph objective and use nodes
    subproblem_optimizer = nothing
end

"""
    SubProblemData{GT<:Plasmo.AbstractOptiGraph}

Holds data specific to a subproblem in the optimization graph.

Fields:
  - `restricted_subgraph::GT`: The subgraph representing the "restricted" subproblem. This is the subproblem before adding overlap.
  - `incident_variable_map::OrderedDict{NodeVariableRef,GT}`: Maps incident variables to their owning subgraphs.
  - `incident_constraint_map::OrderedDict{EdgeConstraintRef,GT}`: Maps incident constraints to their owning subgraphs.
  - `primal_values::OrderedDict{NodeVariableRef,Float64}`: Stores the current primal variable values.
  - `dual_values::OrderedDict{EdgeConstraintRef,Float64}`: Stores the current dual variable values.
  - `primal_parameters::OrderedDict{NodeVariableRef,NodeVariableRef}`: Maps primal variables to their parameter values.
  - `dual_parameters::OrderedDict{EdgeConstraintRef,NodeVariableRef}`: Maps dual variables to their parameter values.
  - `node_objectives::OrderedDict{OptiNode,Plasmo.AbstractJuMPScalar}`: Objective functions associated with nodes in the subproblem.
  - `objective_function::Union{Nothing,Plasmo.AbstractJuMPScalar}`: The total objective function for the subproblem, including penalties.
  - `last_termination_status::MOI.TerminationStatusCode`: The termination status from the last solve.
"""
mutable struct SubProblemData{GT<:Plasmo.AbstractOptiGraph}
    restricted_subgraph::GT

    # where a subproblem should look for primal and dual values
    incident_variable_map::OrderedDict{NodeVariableRef,GT}
    incident_constraint_map::OrderedDict{EdgeConstraintRef,GT}

    # the current primal and dual values
    primal_values::OrderedDict{NodeVariableRef,Float64}
    dual_values::OrderedDict{EdgeConstraintRef,Float64}

    # primal and dual parameters (used for updating the objective)
    primal_parameters::OrderedDict{NodeVariableRef,NodeVariableRef}
    dual_parameters::OrderedDict{EdgeConstraintRef,NodeVariableRef}

    # node objectives on subproblem
    node_objectives::OrderedDict{OptiNode,Plasmo.AbstractJuMPScalar}

    # the current total objective function (with penalties)
    objective_function::Union{Nothing,Plasmo.AbstractJuMPScalar}

    # subproblem termination status
    last_termination_status::MOI.TerminationStatusCode
end
function SubProblemData(restricted_subgraph::GT) where {GT<:Plasmo.AbstractOptiGraph}
    incident_variable_map = OrderedDict{NodeVariableRef,GT}()
    incident_constraint_map = OrderedDict{EdgeConstraintRef,GT}()
    primal_values = OrderedDict{NodeVariableRef,Float64}()
    dual_values = OrderedDict{EdgeConstraintRef,Float64}()
    primal_parameters = OrderedDict{NodeVariableRef,NodeVariableRef}()
    dual_parameters = OrderedDict{EdgeConstraintRef,NodeVariableRef}()
    node_objectives = OrderedDict{OptiNode,Plasmo.AbstractJuMPScalar}()
    return SubProblemData(
        restricted_subgraph,
        incident_variable_map,
        incident_constraint_map,
        primal_values,
        dual_values,
        primal_parameters,
        dual_parameters,
        node_objectives,
        nothing,
        MOI.OPTIMIZE_NOT_CALLED,
    )
end

"""
    SchwarzAlgorithm{GT<:Plasmo.AbstractOptiGraph}

Represents a Schwarz-based optimization algorithm applied to a partitioned graph.

Fields:
  - `graph::GT`: The global optimization graph passed to the algorithm.
  - `subproblems::Vector{GT}`: A list of subgraphs representing the (expanded) subproblems.
  - `element_subproblem_map::Dict{Plasmo.OptiElement,GT}`: Maps elements in the graph to their associated subproblems.
  - `objective_func::Plasmo.AbstractJuMPScalar`: The global objective function for the graph.
  - `options::Options`: Configuration options for the algorithm.
  - `initialized::Bool`: Whether the algorithm has been initialized.
  - `status::MOI.TerminationStatusCode`: The current status of the algorithm.
  - `err_pr::Union{Nothing,Float64}`: Current primal error.
  - `err_du::Union{Nothing,Float64}`: Current dual error.
  - `objective_value::Union{Nothing,Float64}`: Current objective value.
  - `iteration::Int64`: Current iteration count.
  - `primal_error_iters::Vector{Float64}`: History of primal errors per iteration.
  - `dual_error_iters::Vector{Float64}`: History of dual errors per iteration.
  - `objective_iters::Vector{Float64}`: History of objective values per iteration.
  - `solve_time::Float64`: Total solve time.
  - `timers::Timers`: Timers to measure performance metrics.
"""
mutable struct SchwarzAlgorithm{GT<:Plasmo.AbstractOptiGraph}
    graph::GT
    subproblems::Vector{GT}

    # NOTE: we have to track which elements map to which subproblems. this is the only
    # data that has to be centralized on the algorithm.
    element_subproblem_map::Dict{Plasmo.OptiElement,GT}

    # for evaluating total objective
    objective_func::Plasmo.AbstractJuMPScalar

    # algorithm options
    options::Options

    # current status
    initialized::Bool
    status::MOI.TerminationStatusCode
    err_pr::Union{Nothing,Float64}
    err_du::Union{Nothing,Float64}
    objective_value::Union{Nothing,Float64}

    # algorithm iteration outputs
    iteration::Int64
    primal_error_iters::Vector{Float64}
    dual_error_iters::Vector{Float64}
    objective_iters::Vector{Float64}
    solve_time::Float64

    # timers
    timers::Timers

    # inner constructor
    function SchwarzAlgorithm(graph::GT, options::Options) where {GT<:Plasmo.AbstractOptiGraph}
        subproblems = Vector{GT}()
        if options.use_node_objectives
            objective_func = sum(objective_function(node) for node in all_nodes(graph))
        else
            objective_func = Plasmo.objective_function(graph)
        end
        return new{GT}(
            graph,
            subproblems,
            Dict{Plasmo.OptiElement,SubProblemData{GT}}(),
            objective_func,
            options,
            false,
            MOI.OPTIMIZE_NOT_CALLED,
            nothing,
            nothing,
            nothing,
            0,
            Float64[],
            Float64[],
            Float64[],
            0.0,
            Timers(),
        )
    end
end
"""
    SchwarzAlgorithm(graph::OptiGraph, expanded_subgraphs::Vector{OptiGraph}; kwargs...)

Create an algorithm instance by providing the subproblems directly.


    SchwarzAlgorithm(graph::OptiGraph, partition::Plasmo.Partition; kwargs...)

Create an algorithm instance by providing a valid `Plasmo.Partition`.


    SchwarzAlgorithm(graph::OptiGraph, partition::Plasmo.Partition; kwargs...)

Create an algorithm instance and use Metis internally to create partitions.
"""
function SchwarzAlgorithm(graph::OptiGraph, expanded_subgraphs::Vector{OptiGraph}; kwargs...)
    options = Options(; kwargs...)
    algorithm = SchwarzAlgorithm(graph, options)
    algorithm.subproblems = expanded_subgraphs
    return algorithm
end
function SchwarzAlgorithm(graph::OptiGraph, partition::Plasmo.Partition; kwargs...)
    partitioned_graph = assemble_optigraph(partition)
    subgraphs = local_subgraphs(partitioned_graph)
    projection = hyper_projection(graph)
    expanded_subgraphs = expand.(projection, subgraphs, kwargs[:overlap_distance])
    algorithm = SchwarzAlgorithm(partitioned_graph, expanded_subgraphs; kwargs...)
    return algorithm
end
function SchwarzAlgorithm(graph::OptiGraph; n_partitions=4, kwargs...)
    clique_proj = Plasmo.clique_projection(graph)
    metis_vector = Int64.(Metis.partition(clique_proj.projected_graph, n_partitions))
    metis_partition = Partition(clique_proj, metis_vector)
    algorithm = SchwarzAlgorithm(graph, metis_partition; kwargs...)
    return algorithm
end

Base.broadcastable(algorithm::SchwarzAlgorithm) = Ref(algorithm)

# printing

function Base.string(algorithm::SchwarzAlgorithm)
    return @sprintf(
        """
        SchwarzAlgorithm
        %30s %9s
        %30s %9s
        %30s %9s
        %30s %9s
        """,
        "Number of subproblems:",
        length(algorithm.subproblems),
        "Number of variables:",
        Plasmo.num_variables(algorithm.graph),
        "Number of constraints:",
        Plasmo.num_constraints(algorithm.graph),
        "Number of linking constraints:",
        Plasmo.num_local_link_constraints(algorithm.graph)
    )
end
Base.print(io::IO, algorithm::SchwarzAlgorithm) = Base.print(io, Base.string(algorithm))
Base.show(io::IO, algorithm::SchwarzAlgorithm) = Base.print(io, algorithm)

"""
    check_valid_problem(algorithm::SchwarzAlgorithm) -> Bool

Validates the setup of an optimization problem for the specified algorithm. This function performs the following checks:

1. Ensures a subproblem optimizer is defined in the algorithm options.
2. Validates the consistency of subproblem graphs with their respective restricted subgraphs.
3. Confirms that the objective function is separable. Non-separable objectives are currently unsupported.
4. Verifies that the optimization graph does not exhibit a hierarchical structure, which is unsupported.

If some of these checks fail, the algorithm's `status` is set to `MOI.INVALID_MODEL`, and an appropriate error is thrown.

Notes:
- Additional validations planned for future implementations:
  - Checking for non-contiguous partitions.
  - Ensuring there is sufficient overlap (at least 1) in expanded graphs.
"""
function check_valid_problem(algorithm::SchwarzAlgorithm)
    graph = algorithm.graph
    restricted_subgraphs = local_subgraphs(graph)
    n_subproblems = length(algorithm.subproblems)
    if algorithm.options.subproblem_optimizer == nothing
        error(
            "No algorithm set for the subproblems.  Please provide an algorithm constructor to use to solve subproblem optigraphs",
        )
    end

    # check subgraphs are valid
    for i in 1:n_subproblems
        subproblem_graph = algorithm.subproblems[i]
        restricted_subgraph = restricted_subgraphs[i]
        if !(
            intersect(all_nodes(restricted_subgraph), all_nodes(subproblem_graph)) ==
            all_nodes(restricted_subgraph)
        )
            algorithm.status = MOI.INVALID_MODEL
            error("Invalid subproblems given to algorithm.")
        end
    end

    # TODO: check if objective is separable. need custom way to handle non-separable.
    if !(Plasmo.is_separable(objective_function(algorithm.graph)))
        algorithm.status = MOI.INVALID_MODEL
        error("SchwarzAlgorithm does not yet support non-separable objective functions.")
    end

    # TODO: check if graph is hierarchical.
    # TODO: come up with way to handle 'parent' nodes in Schwarz setting
    if _is_hierarchical(graph)
        algorithm.status = MOI.INVALID_MODEL
        error("SchwarzAlgorithm does not yet support optigraphs with a hierarchical structure.")
    end

    # TODO: raise warning for non-contiguous partitions
    # _check_contiguous_partitions

    # TODO: check for an overlap of at least 1 in expanded graphs.
    # _check_overlap

    return true
end

"""
    initialize!(algorithm::SchwarzAlgorithm)

Initializes the algorithm by validating the input problem and setting up subproblems. This includes:
  - Validating the global graph and subproblem consistency.
  - Assigning incident constraints and variables.
  - Setting up objective penalties for the subproblems.

Throws:
  - `ArgumentError` if the algorithm fails to initialize.
"""
function initialize!(algorithm::SchwarzAlgorithm)
    if algorithm.initialized
        error("SchwarzAlgorithm already initialized. Create a new instance of 
        `SchwarzAlgorithm` if you want to restart the algorithm completely.")
    end

    return algorithm.timers.initialize_time = @elapsed begin
        check_valid_problem(algorithm)
        graph = algorithm.graph
        n_subproblems = length(algorithm.subproblems)
        subgraphs = local_subgraphs(graph)

        # create subproblem data for each subproblem
        for i in 1:length(algorithm.subproblems)
            restricted_subgraph = subgraphs[i]
            subproblem = algorithm.subproblems[i]
            subproblem.ext[:subproblem_data] = SubProblemData(restricted_subgraph)
            for node in all_nodes(restricted_subgraph)
                algorithm.element_subproblem_map[node] = subproblem
                # node[:subproblem_graph] = expanded_subgraph
            end
            for edge in all_edges(restricted_subgraph)
                algorithm.element_subproblem_map[edge] = subproblem
                # edge[:subproblem_graph] = expanded_subgraph
            end
        end

        # assign edges that cross/couple subgraphs
        for link_edge in local_edges(algorithm.graph)
            # we assign edge to the first node
            assigned_node = link_edge.nodes[1]
            subproblem = algorithm.element_subproblem_map[assigned_node]
            algorithm.element_subproblem_map[link_edge] = subproblem
        end

        # gather incident constraints
        all_incident_edges = _find_boundary_edges(algorithm.graph, algorithm.subproblems)
        all_incident_constraints = _extract_constraints(all_incident_edges)
        @assert length(all_incident_constraints) == n_subproblems

        # setup subproblems
        for i in 1:n_subproblems
            subproblem = algorithm.subproblems[i]
            subproblem_data = subproblem.ext[:subproblem_data]
            subproblem_incident_edges = all_incident_edges[i]

            # create an auxilliary node to hold parameter variables
            parameter_node = add_node(subproblem)
            for linking_edge in subproblem_incident_edges

                # build primal information map
                linked_variables = all_variables(linking_edge)
                incident_variables = setdiff(linked_variables, all_variables(subproblem))
                for incident_variable in incident_variables
                    owning_node = get_node(incident_variable)
                    # owning_subproblem = owning_node[:subproblem_graph]
                    owning_subproblem = _get_subproblem(algorithm, owning_node)
                    subproblem_data.incident_variable_map[incident_variable] =
                        owning_subproblem
                    if Plasmo.start_value(incident_variable) == nothing
                        primal_start = 1.0
                    else
                        primal_start = Plasmo.start_value(incident_variable)
                    end
                    # add parameter for primal variable value
                    @variable(parameter_node, p in MOI.Parameter(primal_start))
                    subproblem_data.primal_parameters[incident_variable] = p
                    subproblem_data.primal_values[incident_variable] = primal_start
                end

                # build dual information map
                for link_constraint in all_constraints(linking_edge)
                    owning_edge = get_edge(link_constraint)
                    owning_subproblem = _get_subproblem(algorithm, owning_edge)
                    # owning_subproblem = owning_edge[:subproblem_graph]
                    subproblem_data.incident_constraint_map[link_constraint] =
                        owning_subproblem
                    if Plasmo.start_value(link_constraint) == nothing
                        dual_start = 1.0
                    else
                        dual_start = Plasmo.start_value(link_constraint)
                    end
                    # add parameter for dual variable value
                    @variable(parameter_node, d in MOI.Parameter(dual_start))
                    subproblem_data.dual_parameters[link_constraint] = d
                    subproblem_data.dual_values[link_constraint] = dual_start
                end
            end
        end

        _initialize_subproblem_objectives(algorithm)

        reset_iterations(algorithm)

        algorithm.initialized = true
    end
end

"""
    _get_subproblem(algorithm::SchwarzAlgorithm, element::Plasmo.OptiElement) -> GT

Retrieves the subproblem associated with a given element in the optimization graph.
"""
function _get_subproblem(algorithm::SchwarzAlgorithm, element::Plasmo.OptiElement)
    return algorithm.element_subproblem_map[element]
end

"""
    _initialize_subproblem_objectives(algorithm::SchwarzAlgorithm)

Initializes the objective functions for each subproblem in the algorithm. 

This function either sets node objectives directly, if specified, or extracts
and assigns the graph's separable objective terms to the subproblems. 

Additionally, augmented Lagrangian penalty terms are added to the objectives.
"""
function _initialize_subproblem_objectives(algorithm::SchwarzAlgorithm)
    if algorithm.options.use_node_objectives
        for expanded_subgraph in algorithm.subproblems
            Plasmo.set_to_node_objectives(expanded_subgraph)
        end
    else
        # extract objective terms from graph if we are not using the node objectives
        _extract_node_objectives(algorithm)

        # set the objective on each subproblem
        for expanded_subgraph in algorithm.subproblems
            node_objectives = expanded_subgraph.ext[:subproblem_data].node_objectives
            set_objective_function(
                expanded_subgraph,
                sum(node_objectives[node] for node in all_nodes(expanded_subgraph)),
            )
            set_objective_sense(expanded_subgraph, Plasmo.objective_sense(algorithm.graph))
        end
    end

    # add objective penalties
    for expanded_subgraph in algorithm.subproblems
        _formulate_objective_penalty(expanded_subgraph, algorithm.options.mu)
    end
    return nothing
end

"""
    _extract_node_objectives(algorithm::SchwarzAlgorithm)

Extracts separable objective terms from the global graph and assigns them to the nodes 
of each subproblem's restricted subgraph.
"""
function _extract_node_objectives(algorithm::SchwarzAlgorithm)
    objective_func = Plasmo.objective_function(algorithm.graph)
    for expanded_subgraph in algorithm.subproblems
        subproblem_data = expanded_subgraph.ext[:subproblem_data]
        restricted_subgraph = subproblem_data.restricted_subgraph

        # extract objective terms that are on this subgraph
        node_objectives = Plasmo.extract_separable_terms(
            objective_func, restricted_subgraph
        )
        for node in all_nodes(restricted_subgraph)
            subproblem_data.node_objectives[node] = sum(node_objectives[node])
        end
    end
    return nothing
end

"""
    _formulate_objective_penalty(expanded_subgraph::OptiGraph, mu::Float64)

Adds dual penalties and an augmented Lagrangian penalty term to the objective of the 
given subproblem.

Args:
  - `expanded_subgraph`: The subproblem graph with overlap regions.
  - `mu`: The penalty parameter for the augmented Lagrangian method.
"""
function _formulate_objective_penalty(expanded_subgraph::OptiGraph, mu::Float64)
    subproblem_data = expanded_subgraph.ext[:subproblem_data]
    # variable swap function for generating penalty terms
    incident_variables = keys(subproblem_data.incident_variable_map)
    function swap_variable_func(nvref::Plasmo.NodeVariableRef)
        if nvref in incident_variables
            return subproblem_data.primal_parameters[nvref]
        else
            return nvref
        end
    end

    # start with subproblem objective
    subproblem_objective = Plasmo.objective_function(expanded_subgraph)

    # add penalty for each incident linking constraint
    for link_reference in keys(subproblem_data.incident_constraint_map)
        link_constraint = Plasmo.constraint_object(link_reference)
        link_func = Plasmo.jump_function(link_constraint)
        link_rhs = Plasmo.moi_set(link_constraint).value

        # generate penalty term by swapping incident terms with their values
        penalty_term = Plasmo.value(swap_variable_func, link_func) - link_rhs
        augmented_penalty_term = penalty_term^2

        # update objective with dual term and augmented penalty
        link_dual_parameter = subproblem_data.dual_parameters[link_reference]
        subproblem_objective += *(-1, link_dual_parameter * penalty_term)
        subproblem_objective += *(0.5 * mu, augmented_penalty_term)

        # NOTE: `add_to_expression` does not work on nonlinear terms
        # Plasmo.add_to_expression!(obj, -1, link_dual*penalty_term)
        # Plasmo.add_to_expression!(obj, 0.5*mu, augmented_penalty_term)
    end
    set_objective_function(expanded_subgraph, subproblem_objective)
    return nothing
end

"""
    reset_iterations(algorithm::SchwarzAlgorithm)

Resets the algorithm's iteration counters and initializes the timers.
Also sets up optimizers for each subproblem.
"""
function reset_iterations(algorithm::SchwarzAlgorithm)
    algorithm.timers = Timers()
    algorithm.timers.start_time = time()

    # initialize subproblem algorithms
    for subgraph in algorithm.subproblems
        Plasmo.set_optimizer(subgraph, algorithm.options.subproblem_optimizer)
    end

    algorithm.err_pr = Inf
    algorithm.err_du = Inf
    algorithm.iteration = 0
    return nothing
end

"""
    do_iteration(algorithm::SchwarzAlgorithm)

Performs a single iteration of the algorithm. This involves:
1. Solving each subproblem.
2. Communicating primal and dual values between subproblems.
3. Updating subproblem objectives, penalties, and neighbor values.
"""
function do_iteration(algorithm::SchwarzAlgorithm)
    algorithm.timers.update_subproblem_time += @elapsed begin
        Threads.@threads for subproblem_graph in algorithm.subproblems
            _update_subproblem(subproblem_graph)
        end
    end

    # update based on current information
    algorithm.timers.solve_subproblem_time += @elapsed begin
        Threads.@threads for subproblem_graph in algorithm.subproblems
            _solve_subproblem(subproblem_graph)
        end
    end

    # communicate primal and dual values
    algorithm.timers.communicate_time += @elapsed begin
        Threads.@threads for subproblem_graph in algorithm.subproblems
            _retrieve_neighbor_values(subproblem_graph)
        end
    end

    return nothing
end

"""
    _solve_subproblem(subproblem_graph::OptiGraph)

Solves an individual subproblem and logs its termination status.
"""
function _solve_subproblem(subproblem_graph::OptiGraph)
    Plasmo.optimize!(subproblem_graph)
    term_status = Plasmo.termination_status(subproblem_graph)
    !(
        term_status in [
            MOI.TerminationStatusCode(4),
            MOI.TerminationStatusCode(1),
            MOI.TerminationStatusCode(10),
        ]
    ) && @warn("Suboptimal solution detected for subproblem with status $term_status")
    subproblem_graph.ext[:subproblem_data].last_termination_status = term_status
    return nothing
end

"""
    _retrieve_neighbor_values(subproblem_graph::OptiGraph)

Fetches primal and dual values from neighboring subproblems linked via incident
variables and constraints.
"""
function _retrieve_neighbor_values(subproblem_graph::OptiGraph)
    subproblem_data = subproblem_graph.ext[:subproblem_data]
    for (variable, remote_graph) in subproblem_data.incident_variable_map
        subproblem_data.primal_values[variable] = Plasmo.value(remote_graph, variable)
    end
    for (link_reference, remote_graph) in subproblem_data.incident_constraint_map
        subproblem_data.dual_values[link_reference] = Plasmo.dual(
            remote_graph, link_reference
        )
    end
    return nothing
end

"""
    _update_suproblem(subproblem_graph::OptiGraph)

Updates the objective penalty terms of the given subproblem based on the current
neighbor values.
"""
function _update_subproblem(subproblem_graph::OptiGraph)
    if Plasmo.primal_status(subproblem_graph; result=1) != MOI.NO_SOLUTION
        var_vals = value.(Ref(subproblem_graph), all_variables(subproblem_graph))
        for (i, var) in enumerate(all_variables(subproblem_graph))
            if !(Plasmo.is_parameter(var))
                Plasmo.set_start_value(var, var_vals[i])
            end
        end
    end
    _update_objective_penalty(subproblem_graph)
    return nothing
end

"""
    _update_objective_penalty(expanded_subgraph::OptiGraph)

Updates the terms in the subproblem's objective based on the current primal and dual values.
"""
function _update_objective_penalty(expanded_subgraph::OptiGraph)
    subproblem_data = expanded_subgraph.ext[:subproblem_data]
    for variable in keys(subproblem_data.incident_variable_map)
        variable_primal_value = subproblem_data.primal_values[variable]
        variable_parameter = subproblem_data.primal_parameters[variable]
        Plasmo.set_parameter_value(variable_parameter, variable_primal_value)
    end

    for link_reference in keys(subproblem_data.incident_constraint_map)
        link_dual_value = subproblem_data.dual_values[link_reference]
        link_dual_parameter = subproblem_data.dual_parameters[link_reference]
        Plasmo.set_parameter_value(link_dual_parameter, link_dual_value)
    end
    return nothing
end

"""
    calculate_objective_value(algorithm::SchwarzAlgorithm) -> Float64

Evaluate the objective value defined over the algorithm's graph. The function uses 
the current variable values from the subproblems to compute the overall objective.

Args:
  - `algorithm::SchwarzAlgorithm`: The algorithm object containing the objective function 
    and subproblem mappings.

Returns:
    The evaluated objective value as a `Float64`.
"""
function calculate_objective_value(algorithm::SchwarzAlgorithm)
    # grab variable values from corresponding subproblems
    vars = Plasmo.extract_variables(algorithm.objective_func)
    var_vals = Dict{NodeVariableRef,Float64}()
    for var in vars
        node = get_node(var)
        subproblem_graph = _get_subproblem(algorithm, node)
        #subproblem_graph = node[:subproblem_graph]
        var_vals[var] = Plasmo.value(subproblem_graph, var)
    end

    # evaluate using restricted solutions
    objective_value = Plasmo.value(i -> get(var_vals, i, 0.0), algorithm.objective_func)
    algorithm.objective_value = objective_value
    return objective_value
end

"""
    calculate_primal_feasibility(algorithm::SchwarzAlgorithm) -> Vector{Float64}

Evaluate the primal feasibility of the linking constraints defined over the 
algorithm's graph. This is done by checking the residuals between the linking 
constraints and their expected values based on the subproblem solutions.

Args:
  - `algorithm::SchwarzAlgorithm`: The algorithm object containing the graph and linking constraints.

Returns:
    A vector of primal residuals for the linking constraints.
"""
function calculate_primal_feasibility(algorithm::SchwarzAlgorithm)
    link_constraint_refs = local_link_constraints(algorithm.graph)
    n_link_constraints = length(link_constraint_refs)
    primal_residuals = Vector{Float64}(undef, n_link_constraints)
    for i in 1:n_link_constraints
        link_reference = link_constraint_refs[i]
        link_constraint = Plasmo.constraint_object(link_reference)
        func = Plasmo.jump_function(link_constraint)
        set = Plasmo.moi_set(link_constraint)

        # grab variable values from corresponding subproblems
        vars = Plasmo.extract_variables(func)
        var_vals = Dict{NodeVariableRef,Float64}()
        for var in vars
            node = get_node(var)
            subproblem = _get_subproblem(algorithm, node)
            # subproblem_graph = node[:subproblem_graph]
            var_vals[var] = Plasmo.value(subproblem, var)
        end

        # evaluate linking constraint using subproblem variable values
        constraint_value = Plasmo.value(i -> get(var_vals, i, 0.0), func)
        primal_residuals[i] = constraint_value - set.value
    end
    return primal_residuals
end

"""
    calculate_dual_feasibility(algorithm::SchwarzAlgorithm) -> Vector{Float64}

Evaluate the dual feasibility of the linking constraints defined over the 
algorithm's graph. For each linking constraint, the function calculates the 
difference between the maximum and minimum dual values across subproblems.

Args:
  - `algorithm::SchwarzAlgorithm`: The algorithm object containing the graph and linking constraints.

Returns:
    A vector of dual residuals for the linking constraints.
"""
function calculate_dual_feasibility(algorithm::SchwarzAlgorithm)
    graph = algorithm.graph
    link_constraint_refs = local_link_constraints(algorithm.graph)
    n_link_constraints = length(link_constraint_refs)
    dual_residuals = Vector{Float64}(undef, n_link_constraints)
    for i in 1:n_link_constraints
        link_reference = link_constraint_refs[i]
        edge = Plasmo.get_edge(link_reference)

        # get all subproblems that contain this edge
        graphs = Set{typeof(algorithm.graph)}()
        for node in all_nodes(edge)
            subproblem_graph = _get_subproblem(algorithm, node)
            #subproblem_graph = node[:subproblem_graph]
            push!(graphs, subproblem_graph)
        end

        # check each subproblem's dual value for this linkconstraint
        subproblem_duals = Vector{Float64}(undef, length(graphs))
        for (i, graph) in enumerate(collect(graphs))
            subproblem_duals[i] = dual(graph, link_reference)
        end
        dual_residuals[i] = abs(maximum(subproblem_duals) - minimum(subproblem_duals))
    end
    return dual_residuals
end

"""
    eval_iteration(algorithm::SchwarzAlgorithm; save_iteration::Bool=true)

Evaluate the current iterate of the algorithm by calculating primal and dual 
feasibility, as well as the objective value. Optionally, the iterate values 
can be stored internally for tracking progress.

Args:
  - `algorithm::SchwarzAlgorithm`: The algorithm object representing the current optimization state.
  - `save_iteration::Bool`: Whether to save the iterate values. Default is `true`.

Returns:
    Nothing. Updates the algorithm state with current feasibility and objective values.
"""
function eval_iteration(algorithm::SchwarzAlgorithm; save_iteration=true)
    # evaluate residuals
    algorithm.timers.eval_primal_feasibility_time += @elapsed prf = calculate_primal_feasibility(
        algorithm
    )
    algorithm.timers.eval_dual_feasibility_time += @elapsed duf = calculate_dual_feasibility(
        algorithm
    )
    algorithm.err_pr = norm(prf[:], Inf)
    algorithm.err_du = norm(duf[:], Inf)

    # eval objective
    algorithm.timers.eval_objective_time += @elapsed algorithm.objective_value = calculate_objective_value(
        algorithm
    )

    # save iteration data
    if save_iteration
        push!(algorithm.primal_error_iters, algorithm.err_pr)
        push!(algorithm.dual_error_iters, algorithm.err_du)
        push!(algorithm.objective_iters, algorithm.objective_value)
    end
    return nothing
end

"""
    _check_tolerance(algorithm::SchwarzAlgorithm) -> Bool

Check if the primal or dual residuals exceed the specified tolerance. This 
determines whether the algorithm should continue iterating.
"""
function _check_tolerance(algorithm::SchwarzAlgorithm)
    if algorithm.err_pr > algorithm.options.tolerance
        return true
    elseif algorithm.err_du > algorithm.options.tolerance
        return true
    end
    return false
end

"""
    run_algorithm!(algorithm::SchwarzAlgorithm)

Execute the optimization algorithm, iteratively solving subproblems and updating 
the solution until convergence or a stopping criterion is met.
"""
function run_algorithm!(algorithm::SchwarzAlgorithm)
    if !algorithm.initialized
        println("Initializing SchwarzAlgorithm...")
        initialize!(algorithm)
    end

    println("###########################################################")
    println("Running Algorithm: PlasmoSchwarz v0.1.0 using $(Threads.nthreads()) threads")
    println("###########################################################")
    println()
    println("Number of variables: $(num_variables(algorithm.graph))")
    println("Number of constraints: $(num_constraints(algorithm.graph))")
    println("Number of subproblems: $(length(algorithm.subproblems))")
    println("Subproblem variables:   $(num_variables.(algorithm.subproblems))")
    println("Subproblem constraints: $(num_constraints.(algorithm.subproblems))")
    println()

    reset_iterations(algorithm)

    while _check_tolerance(algorithm)
        algorithm.iteration += 1
        if algorithm.iteration > algorithm.options.max_iterations
            algorithm.status = MOI.ITERATION_LIMIT
            break
        end

        do_iteration(algorithm)

        eval_iteration(algorithm)

        # print iteration
        if algorithm.iteration % 20 == 0 || algorithm.iteration == 1
            @printf "%4s | %8s | %8s | %8s" "Iter" "Obj" "Prf" "Duf\n"
        end
        @printf(
            "%4i | %7.2e | %7.2e | %7.2e\n",
            algorithm.iteration,
            algorithm.objective_value,
            algorithm.err_pr,
            algorithm.err_du
        )

        # # update start values
        # algorithm.timers.update_subproblem_time += @elapsed begin
        #     for subproblem in algorithm.subproblems
        #         JuMP.set_start_value.(
        #             Ref(subproblem),
        #             all_variables(subproblem),
        #             value.(Ref(subproblem), all_variables(subproblem)),
        #         )
        #     end
        # end
    end

    algorithm.timers.total_time = time() - algorithm.timers.start_time
    if algorithm.status != MOI.ITERATION_LIMIT
        algorithm.status =
            algorithm.subproblems[1].ext[:subproblem_data].last_termination_status
    end

    println()
    println("Number of Iterations: ", algorithm.iteration)
    @printf "%8s | %8s | %8s" "Obj" "Prf" "Duf\n"
    @printf(
        "%7.2e | %7.2e | %7.2e\n",
        algorithm.objective_value,
        algorithm.err_pr,
        algorithm.err_du
    )
    println()
    println("Time spent in subproblems: ", algorithm.timers.solve_subproblem_time)
    println("Solution Time: ", algorithm.timers.total_time)
    println("EXIT: PlasmoSchwarz Finished with status: ", algorithm.status)
    return nothing
end

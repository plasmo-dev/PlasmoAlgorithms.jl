# Quickstart

This quickstart gives a brief overview of the functions needed for using PlasmoSchwarz to solve optimization problems defined in Plasmo. 

Once PlasmoSchwarz has been installed, you can use it from a Julia session using 

```julia 
julia> using PlasmoSchwarz
```

PlasmoSchwarz also requires defining a graph with Plasmo, and we will need a solver for the subproblems. We will use Ipopt. These packages will also need to be loaded
```julia
julia> using Plasmo, Ipopt
```

## PlasmoSchwarz Overview

PlasmoSchwarz is an iterative algorithm that solves overlapping subproblems and then shares information between the subproblems to help converge to the overall solution. The structure that PlasmoSchwarz "sees" is either defined by the user by partitioning into subgraphs or the partitioning can be done automatically using packages like `Metis.jl` or `KaHyPar`. The graph is passed to the `SchwarzAlgorithm` function which overlaps the graphs and updates the graphs as needed to apply Schwarz Decomposition. Schwarz Decomposition can then be applied by calling `run_algorithm!`. 

## Optimal Control Example

To show how PlasmoSchwarz can be applied for solving a problem, we will show a simple example of an optimal control problem over 200 time steps. The problem can be written as 

```math
\begin{align*} 
    \min &\; \sum_{t \in \mathcal{T}} x_t^2 + u_t^2\\
    \textrm{s.t.} &\; x_{t+1} = x_t + u_t + d_t\\
    &\; x_0 = \bar{x}_o \\
    &\; x_t \ge 0 \\
    &\; u_t \ge \underline{u}
\end{align*}
```

Here, $x$ are the states, $u$ are the inputs, and $d$ are disturbances. 

We will start by importing the required packages and defining problem parameters

```julia
using Plasmo, Ipopt
using PlasmoSchwarz

T = 200             # number of time points
d = sin.(1:T)       # a disturbance vector
imbalance = 0.1     # partition imbalance
distance = 2        # expansion distance
n_parts = 10        # number of partitions
```

Next, we will create the OptiGraph and nodes. Each state and input variable will be represented by a node.

```julia
# create the optigraph
graph = Plasmo.OptiGraph()
@optinode(graph, state[1:T])
@optinode(graph, control[1:(T - 1)])
for (i, node) in enumerate(state)
    @variable(node, x)
    @constraint(node, x >= 0)
    @objective(node, Min, x^2)
end
for node in control
    @variable(node, u)
    @constraint(node, u >= -1000)
    @objective(node, Min, u^2)
end
```

Next, we will set the initial value of the state and then link the state variables to the previous state and inputs.

```julia
# initial condition
n1 = state[1]
@constraint(n1, n1[:x] == 0)

# dynamics
for i in 1:(T - 1)
    @linkconstraint(graph, state[i + 1][:x] == state[i][:x] + control[i][:u] + d[i])
end
```

Now, we can define the sub-solver that will be used on the individual subproblem subgraphs, and we can then construct the `SchwarzAlgorithm` object.

```julia
# subproblem optimizer
sub_optimizer = Plasmo.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

# optimize using overlapping schwarz decomposition
optimizer = SchwarzAlgorithm(
    graph;
    n_partitions=n_parts,
    overlap_distance=1,
    subproblem_optimizer=sub_optimizer,
    max_iterations=100,
    mu=10.0,
)
```

Here, we have not defined the partitions to be used, so the `SchwarzAlgorithm` function will do the partitioning internally using `Metis.jl` and then overlap each subgraph by a distance of 1. Now, we can run the algorithm by calling 

```julia
run_algorithm!(optimizer)
```

## Querying Solutions

PlasmoSchwarz.jl provides API functions for querying solutions and information from the `SchwarzAlgorithm` object. To query the objective value, termination status, or solve time, we can use
```julia
Plasmo.objective_value(optimizer)
Plasmo.termination_status(optimizer)
Plasmo.solve_time(optimizer)
```
The primal and dual variables can be queried by calling `value` or `dual` with the first argument being the `SchwarzAlgorithm` object
```julia
Plasmo.value(optmizer, state[1][:x])
cons = all_constraints(graph)
Plasmo.dual(optimizer, cons[1])
```

You can also access the primal and dual feasibility vectors by viewing the algorithm methods.
```julia
prf = PlasmoSchwarz.calculate_primal_feasibility(optimizer)
duf = PlasmoSchwarz.calculate_dual_feasibility(optimizer)
```
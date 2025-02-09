# PlasmoSchwarz.jl

## Overview
PlasmoSchwarz.jl implements overlapping Schwarz decomposition for graph-structured optimization problems using the algorithm outlined in this [paper](https://arxiv.org/abs/1810.00491) and in Chapter 6 of [this PhD thesis](https://asset.library.wisc.edu/1711.dl/V2UHW7KSFIKBQ8Q/R/file-e04b2.pdf).
The package works with the graph-based algebraic modeling package [Plasmo.jl](https://github.com/plasmo-dev/Plasmo.jl) to formulate and solve problems.

## Installation
PlasmoSchwarz.jl can be installed using the following Julia Pkg command:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/plasmo-dev/PlasmoAlgorithms.jl/tree/main/lib/PlasmoSchwarz.git"))
```

## Simple Example
The following example solves a long-horizon optimal control problem where `x` are the states and `u` are the controls.

```julia
using Plasmo, Ipopt
using PlasmoSchwarz

T = 200             # number of time points
d = sin.(1:T)       # a disturbance vector
n_parts = 10        # number of partitions

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

# initial condition
n1 = state[1]
@constraint(n1, n1[:x] == 0)

# dynamics
for i in 1:(T - 1)
    @linkconstraint(graph, state[i + 1][:x] == state[i][:x] + control[i][:u] + d[i])
end

# subproblem optimizer
sub_optimizer = Plasmo.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

# generate algorithm instance. the algorithm will use Metis to partition internally.
optimizer = PlasmoSchwarz.Algorithm(
    graph;
    n_partitions=n_parts,
    subproblem_optimizer=sub_optimizer,
    max_iterations=100,
    mu=10.0,            # augmented lagrangian penalty
)

# run the optimizer
PlasmoSchwarz.run_algorithm!(optimizer)

# check termination status
@show Plasmo.termination_status(optimizer)

# check objective value
@show Plasmo.objective_value(optimizer)

# check first state and control values
@show Plasmo.value.(optimizer, state[1][:x])
@show Plasmo.value.(optimizer, control[1][:u])
```

### Providing Custom Partitions
It is also possible for users to provide custom partitions (in the form of a `Plasmo.Partition`). 

```julia
using KaHyPar
imbalance = 0.1             # partition imbalance
overlap_distance = 2        # expansion distance

# create a hypergraph projection and partition with KaHyPar
projection = Plasmo.hyper_projection(graph)
partition_vector = KaHyPar.partition(
    projection, 
    n_parts; 
    imbalance=imbalance, 
    configuration=:edge_cut
)

# create a `Plasmo.Partition` object using produced vector
partition = Plasmo.Partition(projection, partition_vector)

# run optimizer using provided partition (overlap will run internally)
optimizer = PlasmoSchwarz.Algorithm(
    graph,
    partition;
    overlap_distance=overlap_distance,
    subproblem_optimizer=sub_optimizer,
    max_iterations=100,
    mu=10.0,            # augmented lagrangian penalty
)

PlasmoSchwarz.run_algorithm!(optimizer)

# check termination status
@show Plasmo.termination_status(optimizer)
```

### Providing Custom Subproblems
Users may further provide their own custom overlap and provide  subproblems directly to the algorithm.

```julia
# create custom partitioned optigraph
partitioned_graph = Plasmo.assemble_optigraph(partition)

# generate custom subproblems using overlap distance
subgraphs = Plasmo.local_subgraphs(partitioned_graph)
expanded_subgraphs = Plasmo.expand.(projection, subgraphs, overlap_distance)

# run optimizer using provided subproblem
optimizer = PlasmoSchwarz.Algorithm(
    partitioned_graph,
    expanded_subgraphs;
    subproblem_optimizer=sub_optimizer,
    max_iterations=100,
    mu=10.0,
)

PlasmoSchwarz.run_algorithm!(optimizer)

# check termination status
@show Plasmo.termination_status(optimizer)
```

## Important Notes
- PlasmoSchwarz.jl does not yet perform automatic overlap improvement. This means the user needs to provide sufficient overlap to obtain convergence.
- PlasmoSchwarz.jl is not meant for problems with integer decision variables.
- Convergence may fail if the user provides non-contiguous subproblems (partitions), which means a subproblem contains distinct sets of unconnected nodes.

## Citing PlasmoSchwarz.jl

If you find Plasmo.jl useful for your work, you may cite the [manuscript](https://link.springer.com/article/10.1007/s12532-022-00223-3) as:
```
@article{Jalving2022,
  title={A Graph-Based Modeling Abstraction for Optimization: Concepts and Implementation in Plasmo.jl},
  author={Jordan Jalving and Sungho Shin and Victor M. Zavala},
  journal={Mathematical Programming Computation},
  year={2022},
  volume={14},
  pages={699 - 747},
  doi={10.1007/s12532-022-00223-3}
}
```
You can also access a freely available [pre-print](https://arxiv.org/abs/2006.05378).

There are also details on this algorithm in Prof. Sungho Shin's PhD dissertation
[here](https://asset.library.wisc.edu/1711.dl/V2UHW7KSFIKBQ8Q/R/file-e04b2.pdf):
``` sourceCode
@book{shin2021graph,
  title={Graph-Structured Nonlinear Programming: Properties and Algorithms},
  author={Shin, Sungho},
  year={2021},
  publisher={The University of Wisconsin-Madison}
}
```
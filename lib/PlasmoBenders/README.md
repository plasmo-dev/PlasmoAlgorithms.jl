# PlasmoBenders.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://plasmo-dev.github.io/PlasmoAlgorithms.jl/dev/PlasmoBenders/introduction/)


PlasmoBenders.jl implements Benders and Nested Benders decomposition to graph-based problems constructed in [Plasmo.jl](https://github.com/plasmo-dev/Plasmo.jl). The algorithm used by this package is outlined by the manuscript [Graph-Based Modeling and Decomposition of Hierarchical Optimization Problems](https://arxiv.org/pdf/2501.02098). After a user defines a graph in Plasmo.jl, they can pass the graph to PlasmoBenders.jl, set the root subgraph, and have PlasmoBenders apply Benders or Nested Benders decomposition to find a solution. These algorithms are applied *based on the user's model*, so algorithm performance is dependent in part on how the user constructed the problem in Plasmo.jl. Further, Benders and Nested Benders requires a specific structure (outlined below) which the user must supply. 

### Installation

PlasmoBenders.jl can be installed using the following Julia Pkg command: 

```julia
using Pkg
Pkg.add("PlasmoBenders")
```

### Overview

[Benders decomposition](https://en.wikipedia.org/wiki/Benders_decomposition) (BD) is a decomposition approach that breaks problems into a master problem and a subproblem(s) and is typically applied to linear and mixed integer linear programs. BD is an iterative algorithm that can be useful for problems where there are a set of complicating variables (in the master problem) that, once fixed, make the subproblem easier to solve. An iteration of BD generally includes 1) solving the master problem, 2) passing the solution of the master problem to the subproblem, 3) solving the subproblem with the master problem solution, and 4) passing primal and dual information from the subproblem to the master problem and forming cutting planes on the master problem. PlasmoBenders applies this approach to graph-based problems where each subgraph of an OptiGraph is the master problem or a subproblem. PD is applied to the graph based on the user-defined subproblems. 

Nested Benders Decomposition (NBD; also called [dual dynamic programming](https://www-sciencedirect-com.ezproxy.library.wisc.edu/science/article/pii/S0098135421000430)) uses similar ideas to BD but can have a sequence of subproblems (i.e., not all subproblems are connected to the original master problem, forming a nested structure). NBD can be applied to graphs with a tree structure, and NBD is likewise implemented in PlasmoBenders.

> [!NOTE]  
> PlasmoBenders requires Plasmo v0.6.2 or later

PlasmoBenders is built on a `BendersOptimizer` object which requires a user-defined graph and a subgraph of that graph as the root (master problem) graph. After the information is past, the `BendersOptimizer` constructor updates the graph to apply BD or NBD. `JuMP.optimize!` is extended so that the iterative BD/NBD algorithm is applied to find a solution. 


### Simple Exmaple
The following example shows how to apply BD/NBD to a graph-based problem using PlasmoBenders. Note that it requires using subgraphs as the subproblems, not nodes as subproblems. Below, we have a master problem with two subproblems. Each subproblem contains one node, but could have more than one.

```julia
using Plasmo, HiGHS, PlasmoBenders

g0 = OptiGraph()
g1 = OptiGraph()
g21 = OptiGraph()
g22 = OptiGraph()

@optinode(g1, n)
@variable(n, x[1:2] >= 0)
@constraint(n, 2 * x[1] + x[2] == 3)
@objective(n, Min, x[1] + 2 * x[2])

@optinode(g21, n)
@variable(n, x >= 0)
@objective(n, Min, 4 * x)

@optinode(g22, n)
@variable(n, x >= 1)
@objective(n, Min, x)

add_subgraph!(g0, g1)
add_subgraph!(g0, g21)
add_subgraph!(g0, g22)

set_to_node_objectives(g1)
set_to_node_objectives(g21)
set_to_node_objectives(g22)

@linkconstraint(g0, g1[:n][:x][2] + g21[:n][:x] == 3)
@linkconstraint(g0, g1[:n][:x][1] + g22[:n][:x] >= 1)

solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)

BendersOpt = BendersOptimizer(g0, g1; solver = solver)

optimize!(BendersOpt)
```

The `BendersOptimizer` constructor takes the overall graph, `g0` as the first object and then the root/master subgraph `g1` as the second argument. 

### Additional Functionality

PlasmoBenders includes additional functionality. The following keyword arguments can be passed to the `BendersOptimizer` constructor. These include the following key word arguments: 

 * `max_iters` - maximum number of iterations to use
 * `tol` - termination tolerance between upper and lower bounds
 * `solver` - the solver to use for the subproblems (the user can set these outside the constructor and not pass any `solver` object)
 * `strengthened` - whether to use [strengthened Benders cuts](https://link.springer.com/article/10.1007/s10107-018-1249-5) (for MIPs)
 * `multicut` - whether to use multicuts (as opposed to aggregated cuts) when adding cutting planes
 * `regularize` - whether to use a regularization scheme for choosing next solutions; currently only works for BD
 * `parallelize_benders` - whether to parallelize solutions of the subproblems in BD

### Citing PlasmoBenders.jl

If you find PlasmoBenders.jl useful for your work, you may cite the [preprint](https://arxiv.org/pdf/2501.02098) as:
```
@article{cole2025,
  title={Graph-Based Modeling and Decomposition of Hierarchical Optimization Problems},
  author={David L. Cole and Filippo Pecci and Omar J. Guerra and Harsha Gangammanavar and Jesse D. Jenkins and Victor M. Zavala},
  journal={arXiv preprint arXiv:2501.02098},
  year={2025}
}
```
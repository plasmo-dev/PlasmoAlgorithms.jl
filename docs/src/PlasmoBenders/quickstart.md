# Quickstart

This quickstart gives a brief overview of the functions needed for using PlasmoBenders to solve optimization problems defined in Plasmo. The example below will include a simple graph with three subproblems that we can solve with either Benders Decomposition (BD) or Nested Benders Decomposition (NBD). We will define the problem and then show how to create and solve the `BendersOptimizer` object.

Once PlasmoBenders has been installed, you can use it from a Julia session using 

```julia 
julia> using PlasmoBenders
```

PlasmoBenders also requires defining a graph with Plasmo, and we will need a solver for the subproblems. We will use HiGHS. These packages will also need to be loaded
```julia
julia> using Plasmo, HiGHS
```

## Required Structure

BD and NBD can exploit problem _structure_, and they require a specific graph structure in order to be applied. PlasmoBenders uses subgraphs as the subproblems, and the subgraphs must form a _tree structure_. This means that there are no cycles among the subgraphs. In other words, if you start at any subgraph $i$ and move along the edges of the graph to any other subgraph $j$, you cannot return to subgraph $i$ without traversing the same subgraphs as you did to reach subgraph $j$. PlasmoBenders, in constructing the `BendersOptimizer` object, will check structure and error out if the structure is not correct.

## Basic OptiGraph

First, we will define a simple graph with 3 subgraphs. Each subgraph will serve as a subproblem of BD/NBD, and each will contain one node. The code below creates this initial graph. 

```julia
using Plasmo, HiGHS, PlasmoBenders

# Define optigraphs
g = OptiGraph(); g1 = OptiGraph()
g2 = OptiGraph(); g3 = OptiGraph()

# Define problem on first OptiGraph
@optinode(g1, n)
@variable(n, x[1:2] >= 0)
@constraint(n, 2 * x[1] + x[2] == 3)
@objective(n, Min, x[1] + 2 * x[2])

# Define problem on second OptiGraph
@optinode(g2, n)
@variable(n, x >= 0)
@objective(n, Min, 4 * x)

# Define problem on third OptiGraph
@optinode(g3, n)
@variable(n, x >= 1)
@objective(n, Min, x)

# Add these OptiGraphs as subgraphs to the overall OptiGraph, g
add_subgraph!(g, g1); add_subgraph!(g, g2); add_subgraph!(g, g3)

# Set subgraph objectives to be their node's objectives
set_to_node_objectives(g1)
set_to_node_objectives(g2)
set_to_node_objectives(g3)

# Create links between the subgraphs
@linkconstraint(g0, g1[:n][:x][2] + g2[:n][:x] == 3)
@linkconstraint(g0, g1[:n][:x][1] + g3[:n][:x] >= 1)
```

Here, subgraph 1 is connected to subgraphs 2 and 3 by the link constraints. Subgraphs 2 and 3 are not connected to each other

!!! note
    PlasmoBenders currently only works with the `Min` objective sense for subgraphs. If your subgraphs have Maximizations for objectives, it is recommended that you reformulate your objective so that $\max f(x)$ becomes $\min -f(x)$.
    
## BendersOptimizer Object
The `BendersOptimizer` object updates the subgraphs to have the required cost-to-go variables and variable copies (fixed from previous stages). The constructor requires two arguments: 1) the overall graph and 2) the subgraph that will serve as the "root" subgraph. The root subgraph is the "master" problem for BD or the first stage problem for NBD. The `BendersOptimizer` constructor also supports several other key word arguments that are oulined in a separate page of the documentation. These include things like the maximum number of iterations, the convergence tolerance, or whether to use certain other options (e.g., regularization) within the solver. 

Before creating the `BendersOptimizer` object, we note that each subgraph requires a solver set on it. We will use HiGHS below, and we will set the output flag solver option to be `false` (this reduces the overall output and makes the BendersOptimizer output more readable). We can either set the solver manually, or we can pass the solver to the `BendersOptimizer` function using the keyward `solver`. To set the solver manually, we do
```julia
solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
set_optimizer(g1, solver)
set_optimizer(g2, solver)
set_optimizer(g3, solver)
```

We can now call the `BendersOptimizer` function (and if we have not set the solver, we can do so when calling this function).

```julia
solver = optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false)
benders_opt = BendersOptimizer(g, g1, solver = solver)
```

Here, `benders_opt` is the optimizer object. It includes several attributes, including the upper and lower bounds at each iteration, solver options, a mapping of complicating variables to their subproblems, and dual and primal information from each iteration. 

Because PlasmoBenders requires a tree structure, we can actually set _any_ subgraph as the root graph. Setting `g1` as the root graph results in two stages (`g1` in stage 1 and `g2` and `g3` in stage 2), but we could instead solve this problem with three stages by setting either `g2` or `g3` as the root graph FILL THIS IN.

## Optimizing the BendersOptimizer

The `BendersOptimizer` object can be solved directly after initialization by calling `JuMP.optimize!`. This function has been extended from JuMP to solve the `BendersOptimizer` object. The solution approach includes performing the forward and backward pass

```julia
JuMP.optimize!(benders_opt)
```

## Querying Solutions




BendersOptimizer constructor

solving the constructor

Note that PlasmoBenders expects a Minimization

querying solutions
# Graph Structure for PlasmoBenders

The graph structure defined by the user determines how the Benders Decomposition (BD) or Nested Benders Decomposition (NBD) algorithms are applied for the given problem. One of the strengths of Plasmo is that it provides flexibility in structuring problems and provides users with flexibility in building the graphs. Below, we will highlight some basic points about the graph structure chosen by the user. Further details can be found in our preprint [here](https://arxiv.org/abs/2501.02098).

## Subproblem Representation

PlasmoBenders applies BD and NBD based on the graph structure defined by the user, and we note here the objects that PlasmoBenders.jl uses for subproblems. Under previous versions of Plasmo (v0.5.4 and earlier), both nodes and subgraphs could be solved independently. Newer releases of Plasmo (v0.6 and on) are more sub-graph centric, and subgraphs are intended to be the "subproblem" objects. Plasmo supports nesting graphs within other graphs, such that you can have any "level" of nested graphs in a given OptiGraph. Any subgraph can be "aggregated" to form a node (i.e., the optimization problem on a graph can equivalently be represnted with a single node). For a given OptiGraph, $\mathcal{G}$, PlasmoBenders expects the graph to be partitioned into subgraphs (i.e., all nodes in the $\mathcal{G}$ are contained on a subgraph), and these subgraphs are treated as the subproblems. Note that each of these subgraphs could likewise have nested subgraphs within them, but PlasmoBenders does not operate on these additional nested subgraphs. Thus, in the `BendersAlgorithm`, there are two required arguments: the optigraph $\mathcal{G}$ and a subgraph of $\mathcal{G}$ which will serve as the root subgraph. 

It is possible that a user may want to use nodes as subproblems instead of subgraphs. PlasmoBenders provides a constructor function for this case and expects each node in the graph to be a subproblem. In this case, the user can pass a node of the graph to the `BendersAlgorithm` constructor function. However, PlasmoBenders is designed to only operate with subgraphs, so in this node-based case, the graph will be partitioned into subgraphs with each node being on a separate subgraph. 

## Required Graph Structure

Both BD and NBD require the subgraphs to form a tree structure. This means that there are no cycles among the subgraphs. In other words, if you start at any subgraph $i$ and move along the edges of the graph to any other subgraph $j$, you cannot return to subgraph $i$ without traversing the same subgraphs as you did to reach subgraph $j$. PlasmoBenders, in constructing the `BendersAlgorithm` object, will check structure and error out if the structure is not correct. The root subgraph specified by the user serves as a "starting point" for the graph structure, and each subgraph will then be placed into stages, with the root subgraph as the first stage. The second stage is all subgraphs connected by an edge to the root subgraph, the third stage is all subgraphs connected to subgraphs in stage two that does not include previous stages, and so forth. BD is a specific case of NBD where there are only two stages. These structures can be visualized as below, where the nodes aligned vertically can be considered in the same stage: 

![graph_structure](../figures/graph_structure.png)
## Building the Graph Structure

Some structures may initially look unfit for BD or NBD, but the flexibility in structuring and partitioning graphs provided by Plasmo can create the required structures. For instance, consider the following graph, where each "node" here can be considered its own subgraph.

![grid graph](../figures/grid_graph.png)

This problem can be partitioned into different subgraphs to form the required tree structure. For instance, you can see below different partitioning schemes that result in different tree structures that can then be used for PlasmoBenders. Note that the dotted lines represent new subgraphs formed by partitioning.

![grid graph partitioned](../figures/grid_graph_partitioning.png)

## Algorithm Performance

Importantly, how the algorithm performs will be dependent on the structure defined by a user. In the above example, each of the given partitionings could result in different algorithm performance. PlasmoBenders generalizes the BD and NBD algorithms to graph structure, but ther performance of these algorithms is dependent in large part on the mathematical problems that exist on the subgraphs. For instance, the choice of root graph will dictate which subgraphs receive cutting planes. Even in the case of two subgraphs connected by an edge, either subgraph could be treated as the root subgraph (or the master problem of BD). However, the master problem of BD is often deliberately chosen (e.g., it may contain all of the integer variables), and choosing one subgraph as the root graph over the other could drastically impact algorithm performance. 

The user is therefore advised to use prior knowledge about their problem in choosing the root graph. Understanding the full impacts of the root graph choice or developing methods for choosing the root graph is an open area that, to our knowledge, has not yet been explored. 
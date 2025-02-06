# PlasmoSchwarz.jl

PlasmoSchwarz.jl is a custom decomposition-based solver that exploits the graph structure defined for a Plasmo.jl OptiGraph and applies an overlapping Schwarz decomposition based on the user-defined structure. The algorithm is outlined in the paper ["A Graph-Based Modeling Abstraction for Optimization: Concepts and Implementation in Plasmo.jl"](https://arxiv.org/abs/2006.05378) and in Prof. Sungho Shin's PhD Dissertation, ["Graph-Structured Nonlinear Programming: Properties and Algorithms"](https://asset.library.wisc.edu/1711.dl/V2UHW7KSFIKBQ8Q/R/file-e04b2.pdf). 

Schwarz decomposition is an iterative algorithm takes advantage of Plasmo.jl's subgraph structuring capabilities by partitioning an optimization problem into separate subgraphs. The subgraphs can then be expanded to overlap with other subgraphs, and a single iteration of the algorithm solves each subgraph, which can be done in parallel. After optimization, solutions around the overlapping regions are shared between subgraphs and used within the algorithm to help the subgraphs converge to the solution of the overall problem. Schwarz decomposition does not support integer variables. Further, in practice, convergence of PlasmoSchwarz is dependent on problem formulation and the degree of overlap between problems. Convergence is not guaranteed, but primal and dual values can be used to detect whether a local solution has been reached. 

## Installation
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/plasmo-dev/PlasmoAlgorithms.jl/tree/main/lib/PlasmoSchwarz.git"))
```

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

## Contents

```@contents
Pages = [
    "quickstart.md"
    "algorithm.md"
    "api_docs.md"
    ]
Depth = 2
```

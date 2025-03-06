# PlasmoBenders.jl

PlasmoBenders.jl is a custom decomposition-based solver that exploits the graph structure defined for a Plasmo.jl OptiGraph and applies [Benders](https://en.wikipedia.org/wiki/Benders_decomposition) or Nested Benders Decomposition based on the user-defined structure (i.e., it detects the copmlicating variables based on the graph structure). The algorithm used by PlasmoBenders.jl is outlined in the manuscript [Graph-Based Modeling and Decomposition of Hierarchical Optimization Problems](https://arxiv.org/pdf/2501.02098). In this way, PlamsoBenders implements these decomposition approaches automatically for the user. 

## Installation
```julia
import Pkg
Pkg.add("PlasmoBenders")
```
or alternatively from the Julia package manager by performing the following:
```
pkg> add PlasmoBenders
```

## Citing PlasmoBenders.jl
If you find PlasmoBenders.jl useful for your work, you may cite the [preprint](https://arxiv.org/abs/2501.02098) as:
```
@article{cole2025,
  title={Graph-Based Modeling and Decomposition of Hierarchical Optimization Problems},
  author={David L. Cole and Filippo Pecci and Omar J. Guerra and Harsha Gangammanavar and Jesse D. Jenkins and Victor M. Zavala},
  journal={arXiv preprint arXiv:2501.02098},
  year={2025}
}
```

## Contents

```@contents
Pages = [
    "algorithm.md"
    "quickstart.md"
    "solver.md"
    "graph_structure.md"
    "api_docs.md"
    "storage_tutorial.md"
    ]
Depth = 2
```

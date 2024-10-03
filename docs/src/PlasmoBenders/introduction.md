# PlasmoBenders.jl

PlasmoBenders.jl is a custom decomposition-based solver that exploits the graph structure defined for a Plasmo.jl OptiGraph and applies [Benders](https://en.wikipedia.org/wiki/Benders_decomposition) or Nested Benders Decomposition based on the user-defined structure (i.e., it detects the copmlicating variables based on the graph structure). In this way, PlamsoBenders implements these decomposition approaches automatically for the user. 

## Installation
```julia
import Pkg
Pkg.add("PlasmoBenders")
```
or alternatively from the Julia package manager by performing the following:
```
pkg> add PlasmoBenders
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

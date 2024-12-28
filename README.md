# PlasmoAlgorithms.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://plasmo-dev.github.io/PlasmoAlgorithms.jl/dev/)

This repository is designed to be a location for meta-algorithms for OptiGraphs with Plasmo.jl. The `\lib` directory contains the algorithms implemented in Plasmo.jl, and these algorithms can be added as packages in Julia.

In the future, we hope to support several decomposition algorithms that exploit the graph structure of Plasmo.jl. Currently, only Benders decomposition is hosted in this directory (as PlasmoBenders.jl), but other algorithms are possible and will likely be added in the future. If you would like to develop a decomposition approach for Plasmo.jl OptiGraphs, we welcome any pull requests and contributions. 

Please report any issues and feature requests via the [Github issue tracker](https://github.com/plasmo-dev/PlasmoAlgorithms.jl/issues).

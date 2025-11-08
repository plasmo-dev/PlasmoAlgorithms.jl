#  Copyright 2024, David Cole, Jordan Jalving, Victor Zavala, and contributors
#  This Source Code Form is subject to the terms of the MIT License
#  This source code is adapted from that of Plasmo.jl which can be found at https://github.com/plasmo-dev/Plasmo.jl/blob/main/docs/make.jl

using Documenter, Plasmo, Suppressor, Graphs, PlasmoBenders, PlasmoSchwarz

DocMeta.setdocmeta!(Plasmo, :DocTestSetup, :(using Plasmo); recursive=true)
DocMeta.setdocmeta!(PlasmoBenders, :DocTestSetup, :(using PlasmoBenders); recursive=true)
DocMeta.setdocmeta!(PlasmoSchwarz, :DocTestSetup, :(using PlasmoSchwarz); recursive=true)

makedocs(;
    sitename="PlasmoAlgorithms.jl",
    modules=[PlasmoBenders],
    doctest=true,
    checkdocs=:export,
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    authors="Jordan Jalving and David Cole",
    pages=[
        "Introduction" => "index.md",
        "PlasmoBenders.jl" => [
            "Introduction" => "PlasmoBenders/introduction.md",
            "Algorithm" => "PlasmoBenders/algorithm.md",
            "Quickstart" => "PlasmoBenders/quickstart.md",
            "Solver Options" => "PlasmoBenders/solver.md",
            "Exploiting Graph Structure" => "PlasmoBenders/graph_structure.md",
            "API Documentation" => "PlasmoBenders/api_docs.md",
            "Tutorials" => [
                "Storage Operation" => "PlasmoBenders/storage_tutorial.md",
                "Equipment Sizing" => "PlasmoBenders/sizing_tutorial.md"
            ],
        ],
        "PlasmoSchwarz.jl" => [
            "Introduction" => "PlasmoSchwarz/introduction.md",
            "Quickstart" => "PlasmoSchwarz/quickstart.md",
            "Algorithm" => "PlasmoSchwarz/algorithm.md",
            "API Documentation" => "PlasmoSchwarz/api_docs.md",
        ]
    ],
)

deploydocs(; repo="github.com/plasmo-dev/PlasmoAlgorithms.jl.git")

#  Copyright 2024, David Cole, Jordan Jalving, Victor Zavala, and contributors
#  This Source Code Form is subject to the terms of the MIT License
#  This source code is adapted from that of Plasmo.jl which can be found at https://github.com/plasmo-dev/Plasmo.jl/blob/main/docs/make.jl

using Documenter, Plasmo, Suppressor, Graphs

doc_pkg = get(ENV, "DOC_PKG", "PlasmoBenders")  # default to PlasmoBenders

if doc_pkg == "PlasmoBenders"
    @info "Building docs for PlasmoBenders"
    using PlasmoBenders
    DocMeta.setdocmeta!(PlasmoBenders, :DocTestSetup, :(using PlasmoBenders); recursive=true)
elseif doc_pkg == "PlasmoSchwarz"
    @info "Building docs for PlasmoSchwarz"
    using PlasmoSchwarz
    DocMeta.setdocmeta!(PlasmoSchwarz, :DocTestSetup, :(using PlasmoSchwarz); recursive=true)
else
    error("Unsupported DOC_PKG=$(doc_pkg). Expected PlasmoBenders or PlasmoSchwarz.")
end

DocMeta.setdocmeta!(Plasmo, :DocTestSetup, :(using Plasmo); recursive=true)
DocMeta.setdocmeta!(Plasmo, :DocTestSetup, :(using Plasmo); recursive=true)

pages_benders = [
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
]

pages_schwarz = [
    "Introduction" => "index.md",
    "PlasmoSchwarz.jl" => [
        "Introduction" => "PlasmoSchwarz/introduction.md",
        "Quickstart" => "PlasmoSchwarz/quickstart.md",
        "Algorithm" => "PlasmoSchwarz/algorithm.md",
        "API Documentation" => "PlasmoSchwarz/api_docs.md",
    ]
]

modules_sel = doc_pkg == "PlasmoBenders" ? [PlasmoBenders] : [PlasmoSchwarz]
pages_sel = doc_pkg == "PlasmoBenders" ? pages_benders : pages_schwarz

makedocs(; 
    sitename="PlasmoAlgorithms.jl - $(doc_pkg)",
    modules=modules_sel,
    doctest=true,
    checkdocs=:export,
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    authors="Jordan Jalving and David Cole",
    pages=pages_sel,
)

devurl = doc_pkg == "PlasmoBenders" ? "benders" : "schwarz"
deploydocs(; repo="github.com/plasmo-dev/PlasmoAlgorithms.jl.git", devurl=devurl)

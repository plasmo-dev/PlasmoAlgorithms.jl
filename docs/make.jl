#  Copyright 2024, David Cole, Jordan Jalving, Victor Zavala, and contributors
#  This Source Code Form is subject to the terms of the MIT License
#  This source code is adapted from that of Plasmo.jl which can be found at https://github.com/plasmo-dev/Plasmo.jl/blob/main/docs/make.jl

using Documenter, Plasmo, Suppressor, Graphs
using Base: walkdir

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

# Build in a filtered temporary source directory to avoid processing pages
# of the other package when building docs for one package.
src_root = joinpath(@__DIR__, "src")
tmp_src = mktempdir()

# Always include index.md
mkpath(tmp_src)
cp(joinpath(src_root, "index.md"), joinpath(tmp_src, "index.md"); force=true)

if doc_pkg == "PlasmoBenders"
    # Copy only PlasmoBenders pages directory
    src_dir = joinpath(src_root, "PlasmoBenders")
    dest_dir = joinpath(tmp_src, "PlasmoBenders")
    mkpath(dest_dir)
    for (root, dirs, files) in walkdir(src_dir)
        rel = replace(root, src_dir => "")
        target_root = joinpath(dest_dir, rel)
        mkpath(target_root)
        for f in files
            cp(joinpath(root, f), joinpath(target_root, f); force=true)
        end
    end
else
    # Copy only PlasmoSchwarz pages directory
    src_dir = joinpath(src_root, "PlasmoSchwarz")
    dest_dir = joinpath(tmp_src, "PlasmoSchwarz")
    mkpath(dest_dir)
    for (root, dirs, files) in walkdir(src_dir)
        rel = replace(root, src_dir => "")
        target_root = joinpath(dest_dir, rel)
        mkpath(target_root)
        for f in files
            cp(joinpath(root, f), joinpath(target_root, f); force=true)
        end
    end
end

makedocs(; 
    sitename="PlasmoAlgorithms.jl - $(doc_pkg)",
    modules=modules_sel,
    doctest=true,
    checkdocs=:export,
    source=tmp_src,
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    authors="Jordan Jalving and David Cole",
    pages=pages_sel,
)

devurl = doc_pkg == "PlasmoBenders" ? "benders" : "schwarz"
deploydocs(; repo="github.com/plasmo-dev/PlasmoAlgorithms.jl.git", devurl=devurl)

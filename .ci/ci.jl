# This file is adapted from MadNLP's source code here:
# https://github.com/MadNLP/MadNLP.jl/blob/master/.ci/ci.jl

# Because the subdirectories include packages that need to be tested, this file
# will perform those tests and is called by the test.yaml file in the workflow

rm(joinpath(@__DIR__, "Project.toml"); force=true)
rm(joinpath(@__DIR__, "Manifest.toml"); force=true)

using Pkg

Pkg.activate(@__DIR__)

# if other packages are added to lib, this file can be updated with those package names
pkgs = ["PlasmoBenders", "PlasmoSchwarz"]

Pkg.develop.([PackageSpec(path=joinpath(@__DIR__,"..","lib",pkg)) for pkg in pkgs])
Pkg.build()

Pkg.test.(pkgs, coverage=true)

using Pkg

PA_DIR = pwd()

Pkg.develop(path=joinpath(PA_DIR, "lib", "PlasmoBenders"))
Pkg.develop(path=joinpath(PA_DIR, "lib", "PlasmoSchwarz"))
Pkg.instantiate()

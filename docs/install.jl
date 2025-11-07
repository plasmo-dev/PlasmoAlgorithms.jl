using Pkg

PA_DIR = pwd()

doc_pkg = get(ENV, "DOC_PKG", nothing)

if doc_pkg === nothing
	# Fallback: develop both (may conflict if dependencies are incompatible)
	@info "DOC_PKG not set; developing both PlasmoBenders and PlasmoSchwarz" 
	Pkg.develop(path=joinpath(PA_DIR, "lib", "PlasmoBenders"))
	Pkg.develop(path=joinpath(PA_DIR, "lib", "PlasmoSchwarz"))
else
	@info "Developing only $(doc_pkg) for docs build"
	Pkg.develop(path=joinpath(PA_DIR, "lib", doc_pkg))
end

Pkg.instantiate()

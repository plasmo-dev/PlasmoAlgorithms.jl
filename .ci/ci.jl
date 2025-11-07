"""
Continuous Integration test driver

This script runs tests for packages under `lib/` in isolated Julia environments
so that conflicting dependency requirements do not interfere with each other.

It is invoked by the GitHub Actions workflow and creates one environment per
package under `.ci/` before running tests. Dependency versions (including
Plasmo) are resolved automatically from each package's compat entries.
"""

using Pkg
using Dates

# Utility to (re)create a clean environment directory
function _clean_env!(envdir::AbstractString)
	if !isdir(envdir)
		mkpath(envdir)
	end
	for f in ("Project.toml", "Manifest.toml")
		fp = joinpath(envdir, f)
		if isfile(fp)
			rm(fp; force=true)
		end
	end
	return envdir
end

"""
	test_package(pkg::AbstractString)

Create a fresh isolated environment for `pkg`, develop the package from
`lib/pkg`, and run its tests with coverage enabled. All dependency versions are
resolved by Julia's package manager according to the package's compat bounds.
"""
function test_package(pkg::AbstractString)
	envdir = joinpath(@__DIR__, "env-" * pkg)
	_clean_env!(envdir)

	# Activate isolated env for this package
	Pkg.activate(envdir)

	# Bring the target package into this environment in development mode
	Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..", "lib", pkg)))

	# Resolve, build, and precompile before testing
	Pkg.instantiate()
	Pkg.build()
	try
		Pkg.precompile()
	catch err
		@warn "Precompile failed; proceeding to tests" exception=(err, catch_backtrace())
	end

	@info "Running tests for $(pkg) (auto-resolved dependencies)" timestamp=Dates.now()
	Pkg.test(pkg; coverage = true)
end

# Allow selection via environment variables for CI job splitting.
# If PKG is set, only test that package. Otherwise, test both sequentially.
pkg_sel = get(ENV, "PKG", nothing)

if pkg_sel !== nothing
	test_package(pkg_sel)
else
	# Run both if not explicitly selected
	test_package("PlasmoBenders")
	test_package("PlasmoSchwarz")
end

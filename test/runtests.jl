using TestItemRunner
using Wannier
using Documenter

# Get filter strings from command line arguments if provided
# Usage: julia --project test/runtests.jl "spread.jl" "bvector.jl"
filter_names = isempty(ARGS) ? nothing : ARGS

if isnothing(filter_names)
    println("Running all tests...")

    @run_package_tests verbose = true

    DocMeta.setdocmeta!(Wannier, :DocTestSetup, :(using Wannier); recursive = true)
    doctest(
        Wannier
        # fix=true,  # update all the output in `jldoctest`
    )
else
    println("Running specific tests: $(join(filter_names, ", "))")

    @run_package_tests verbose = true filter =
        ti -> any(name -> endswith(ti.filename, name), filter_names)
end

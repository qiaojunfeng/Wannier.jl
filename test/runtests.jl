using TestItemRunner
using Wannier
using Documenter

# Get filter strings from command line arguments if provided
# Usage: julia --project test/runtests.jl "spread.jl" "bvector.jl"
const selected_names = isempty(ARGS) ? nothing : ARGS

function filter_tests(test_item)
    if isnothing(selected_names)
        return true  # Run all tests
    else
        return any(selected_names) do sel
            if isdir(sel)
                occursin(sel, test_item.filename)
            else
                endswith(test_item.filename, sel)
            end
        end
    end
end

if isnothing(selected_names)
    println("Running all tests...")
    @run_package_tests verbose = true

    DocMeta.setdocmeta!(Wannier, :DocTestSetup, :(using Wannier); recursive = true)
    doctest(
        Wannier
        # fix=true,  # update all the output in `jldoctest`
    )
else
    println("Running specific tests: $(join(selected_names, ", "))")
    @run_package_tests verbose = true filter = filter_tests
end

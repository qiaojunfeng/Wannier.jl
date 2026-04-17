using TestItemRunner

# Optional flags:
#   --all        run full suite (including plot tests and doctests)
#   --plot       include plotting tests under test/plot/
#   --doctest    run doctests
# Usage examples:
#   julia --project=test test/runtests.jl
#   julia --project=test test/runtests.jl --all
#   julia --project=test test/runtests.jl --plot --doctest
#   julia --project=test test/runtests.jl spread.jl bvector.jl

const run_all = "--all" in ARGS
const run_plot_tests = run_all || ("--plot" in ARGS)
const run_doctest = run_all || ("--doctest" in ARGS)
const selected_names = filter(arg -> !startswith(arg, "--"), ARGS)
const has_selected_names = !run_all && !isempty(selected_names)

function is_plot_test(test_item)
    return occursin("/test/plot/", test_item.filename) ||
           occursin("\\test\\plot\\", test_item.filename)
end

function filter_tests(test_item)
    if has_selected_names
        return any(selected_names) do sel
            if isdir(sel)
                occursin(sel, test_item.filename)
            else
                endswith(test_item.filename, sel)
            end
        end
    end

    # Fast default run excludes plotting tests unless explicitly requested.
    return run_plot_tests || !is_plot_test(test_item)
end

if has_selected_names
    println("Running specific tests: $(join(selected_names, ", "))")
    @run_package_tests verbose = true filter = filter_tests
else
    if run_all
        println("Running full test suite (including plot tests and doctests)...")
    else
        println("Running default test suite (plot tests: $(run_plot_tests ? "enabled" : "disabled"))...")
    end
    @run_package_tests verbose = true filter = filter_tests

    if run_doctest
        println("Running doctests...")
        include("doctest.jl")
    else
        println("Skipping doctests (pass --doctest to enable).")
    end
end

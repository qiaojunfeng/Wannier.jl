# using Aqua
# only test ambiguities in current module, otherwise fails due to ambiguities in other packages
#   https://github.com/JuliaTesting/Aqua.jl/issues/77
# disable project_extras since we don't use julia < 1.2
# Aqua.test_all(Wannier; ambiguities=false, project_extras=false)
# Aqua.test_ambiguities(Wannier)

using TestItemRunner

# TODO enable everything once the refactor is done
# @run_package_tests verbose = true filter = ti -> occursin("exchanges", ti.filename)
@run_package_tests verbose = true filter = ti -> occursin("interpolation/", ti.filename) || occursin("exchanges", ti.filename)

using TestItemRunner

# Filter test files
@run_package_tests verbose = true filter = ti -> occursin("interpolation", ti.filename)

# @run_package_tests verbose = true

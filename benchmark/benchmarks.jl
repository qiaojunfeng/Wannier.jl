using BenchmarkTools

const FIXTURE_PATH = joinpath(@__DIR__, "..", "test", "fixtures")

SUITE = BenchmarkGroup()

for file in readdir(@__DIR__)
    if startswith(file, "bench_") && endswith(file, ".jl")
        name = file[(length("bench_") + 1):(end - length(".jl"))]
        SUITE[name] = include(file)
    end
end

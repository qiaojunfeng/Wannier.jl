using Wannier
using Test

const TEST_PATH = @__DIR__
const FIXTURE_PATH = "$TEST_PATH/fixtures"


@testset "Wannier.jl" begin

    # include("io/w90.jl")

    # include("bvector.jl")

    include("spread.jl")
end

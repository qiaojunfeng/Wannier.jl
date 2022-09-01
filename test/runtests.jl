using Wannier
using Test

using Aqua
# only test ambiguities in current module, otherwise fails due to ambiguities in other packages
#   https://github.com/JuliaTesting/Aqua.jl/issues/77
# disable project_extras since we don't use julia < 1.2
Aqua.test_all(Wannier; ambiguities=false, project_extras=false)
Aqua.test_ambiguities(Wannier)

const TEST_PATH = @__DIR__
const FIXTURE_PATH = joinpath(TEST_PATH, "fixtures")

@testset "Wannier.jl" begin
    include("io/w90.jl")
    include("io/model.jl")
    include("bvector.jl")
    include("spread.jl")
    include("util/kpoint.jl")
    include("util/io.jl")
    include("util/element.jl")
    include("util/center.jl")

    include("wannierize/max_localize.jl")
    include("wannierize/disentangle.jl")
    include("wannierize/opt_rotate.jl")
    include("wannierize/split.jl")
    include("wannierize/parallel_transport.jl")
    include("wannierize/constraint_center/max_localize.jl")
    include("wannierize/constraint_center/disentangle.jl")

    include("cli/help.jl")

    include("interpolate/rvector.jl")
    include("interpolate/fourier.jl")
    include("interpolate/kpath.jl")
    include("interpolate/hamiltonian.jl")
    include("realspace.jl")
end

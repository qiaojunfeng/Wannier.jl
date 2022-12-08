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

mat2vec(A::AbstractMatrix) = [v for v in eachcol(A)]

@testset "Wannier.jl" begin
    include("io/w90.jl")
    include("io/model.jl")
    include("bvector.jl")
    include("spread.jl")
    include("util/kpoint.jl")
    include("util/element.jl")
    include("util/center.jl")
    include("util/kpath.jl")

    include("wannierize/max_localize.jl")
    include("wannierize/disentangle.jl")
    include("wannierize/opt_rotate.jl")
    include("wannierize/split.jl")
    include("wannierize/parallel_transport.jl")
    include("wannierize/constrain_center/max_localize.jl")
    include("wannierize/constrain_center/disentangle.jl")

    include("cli/help.jl")

    include("interp/rvector.jl")
    include("interp/fourier.jl")
    include("interp/hamiltonian.jl")
    include("interp/fermisurf.jl")
    include("interp/derivative.jl")
    include("realspace.jl")
end

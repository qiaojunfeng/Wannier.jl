module Wannier

using Printf: @printf, @sprintf
using LinearAlgebra
using OrderedCollections
using DocStringExtensions
using StaticArrays

using Reexport: @reexport

@reexport using CrystalBase
@reexport using WannierIO

include("common/const.jl")
include("common/size.jl")
include("common/rgrid.jl")

include("utils/printing.jl")
include("utils/linalg.jl")
include("utils/center.jl")

include("defaults.jl")

## Wannierization
include("kpoints/kpoint.jl")
include("kpoints/kstencil_shell.jl")
include("kpoints/kstencil.jl")

include("model.jl")
include("spread.jl")

include("io/w90/nnkp.jl")
include("io/w90/amn.jl")
include("io/w90/mmn.jl")
include("io/w90/chk.jl")
include("io/w90/model.jl")
include("io/truncate.jl")
include("io/interface/mud.jl")
include("io/interface/orbital_magnetization.jl")

include("Datasets.jl")

include("localization/gauge.jl")
include("localization/max_localize.jl")
include("localization/disentangle.jl")
include("localization/opt_rotate.jl")
include("localization/parallel_transport/parallel_transport.jl")
include("localization/split.jl")
include("localization/coopt.jl")
include("localization/constrain_center/coopt.jl")

include("realspace/wavefunction.jl")
include("realspace/moment.jl")

## Wannier interpolation
include("interpolation/Rspace.jl")
include("interpolation/operator.jl")

# these are files for interpolation only, put them here since
# they need structs defined in previous files
include("io/w90/tb.jl")
include("io/w90/hr.jl")
include("io/w90/spn.jl")
include("io/volume/xsf.jl")
include("io/volume/cube.jl")
include("io/volume/bxsf.jl")

include("interpolation/fourier.jl")
include("interpolation/hamiltonian.jl")
include("interpolation/position.jl")
include("interpolation/hamiltonian_gradient.jl")
include("interpolation/hamiltonian_hessian.jl")
include("interpolation/spin.jl")
include("interpolation/berry_curvature.jl")
include("interpolation/orbital_magnetization.jl")
include("interpolation/fermi_energy.jl")

# include("interpolation/magmom.jl")

include("symmetry.jl")

# Some convenience functions for users
include("tools/Tools.jl")

# Plotting functions, requires Makie
include("plot.jl")

function __init__()
    return Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        if exc.f in [bandplot, bandplot!, get_bandplot, projbandplot, projbandplot!, get_projbandplot]
            if isempty(methods(exc.f))
                print(io, "\n$(exc.f) has no methods, yet. Makie has to be loaded for the plotting extension to be activated. Run `using Makie`, `using CairoMakie`, `using GLMakie` or any other package that also loads Makie.")
            end
        end
    end
end

end

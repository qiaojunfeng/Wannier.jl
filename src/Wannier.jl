module Wannier

using Printf: @printf, @sprintf
using LinearAlgebra
using StaticArrays
using OrderedCollections
using DocStringExtensions
using Reexport: @reexport

@reexport using CrystalBase
@reexport using WannierIO

include("common/const.jl")
include("common/size.jl")
include("common/rgrid.jl")
include("common/kernels.jl")
include("common/layouts.jl")

include("utils/printing.jl")
include("utils/linalg.jl")
include("utils/center.jl")

include("defaults.jl")

# Wannierization
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

include("localization/method.jl")
include("localization/gauge.jl")
include("localization/disentangle.jl")
include("localization/parallel_transport/parallel_transport.jl")
include("localization/split.jl")
include("localization/coopt.jl")
include("localization/objective.jl")
include("localization/solver.jl")
include("localization/localize.jl")

include("realspace/wavefunction.jl")
include("realspace/moment.jl")

# Wannier interpolation
include("interpolation/Rspace.jl")
include("interpolation/operator.jl")

include("io/volume/xsf.jl")
include("io/volume/cube.jl")
include("io/volume/bxsf.jl")

include("interpolation/fourier.jl")
include("interpolation/hamiltonian.jl")
include("interpolation/types.jl")
include("interpolation/real_space.jl")
include("interpolation/operators/berry_connection.jl")
include("interpolation/operators/orbital_magnetization.jl")
include("interpolation/construction.jl")
include("io/w90/interpolation.jl")
include("io/w90/tb.jl")
include("io/w90/hr.jl")
include("io/w90/spn.jl")
include("interpolation/fourier_kernel.jl")
include("interpolation/diagonalization.jl")
include("interpolation/planning.jl")
include("interpolation/observables/band_energy.jl")
include("interpolation/observables/band_velocity.jl")
include("interpolation/observables/spin_expectation.jl")
include("interpolation/workflows/band_structure.jl")
include("interpolation/position.jl")
include("interpolation/hamiltonian_gradient.jl")
include("interpolation/hamiltonian_hessian.jl")
include("interpolation/spin.jl")
include("interpolation/berry_curvature.jl")
include("interpolation/observables/berry_curvature.jl")
include("interpolation/orbital_magnetization.jl")
include("interpolation/observables/orbital_magnetization.jl")
include("interpolation/fermi_energy.jl")

# include("interpolation/magmom.jl")

include("symmetry/operations.jl")
include("interpolation/symmetry.jl")
include("symmetry/localization.jl")
include("symmetry/model.jl")

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

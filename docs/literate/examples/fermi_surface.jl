# # Interpolation of Fermi surface

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In this tutorial, we will run the disentanglement algorithm on
the copper and then compute the Fermi surface.

1. generate the `amn`, `mmn`, and `eig` files by using `Quantum ESPRESSO` (QE)
2. construct a [`Model`](@ref) for `Wannier.jl`, by reading the `win`, `amn`, `mmn`, and `eig` files
3. run `Wannier.jl` [`localize`](@ref) on the `Model` to minimize the spread
4. write the maximal localized gauge to a new `amn` file
=#

# ## Preparation
# Load the package
using WannierIO
using Wannier
using Wannier.Datasets
using WGLMakie

#=
!!! tip

    Use the `run.sh` script which automate the scf, nscf, pw2wannier90 steps.
=#

#=
## Model construction

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct a [`Model`](@ref) that abstracts the calculation
=#
model = load_dataset("Cu")

#=
!!! tip

    The [`read_w90`](@ref) function will parse the `win` file and set the frozen window for the `Model` according to
    the `dis_froz_min` and `dis_froz_max` parameters in the `win` file. However, you can also change these parameters
    by calling the [`set_frozen_win!`](@ref) function.
=#

#=
## Disentanglement and maximal localization

The [`localize`](@ref) function
will disentangle and maximally localize the spread
functional, and returns the gauge matrices `U`,
=#
U = localize(model);

# The initial spread is
spread(model)

# The final spread is
spread(model, U)

# save the new gauge to the model
model.gauges .= U;

#=
!!! note

    The convergence thresholds is determined by the
    keyword arguments of [`localize`](@ref), e.g., `f_tol` for the tolerance on spread,
    and `g_tol` for the tolerance on the norm of spread gradient, etc. You can use stricter thresholds
    to further minimize a bit the spread.
=#

# load QE band structure from bands.x output
using QuantumEspressoIO
qe_bands = QuantumEspressoIO.read_band_dat(dataset"Cu/outputs/qe_bands.dat");
# the Fermi energy from scf calculation
εF = 16.8985

#=
## Band structure interpolation
=#
# Auto-generate kpath from Cu crystal structure (requires Spglib + Brillouin).
# Cu.win only has a Monkhorst-Pack grid, no explicit kpoint_path block.
using Spglib
kpath = KPath(KSegment(model), default_w90_kpath_num_points())

interpolation_model = InterpolationModel(model)

# interpolate band structure; E_mat is n_bands × n_kpoints
E_mat = interpolate(
    interpolation_model, collect(kpath), BandEnergy()
).band_energy;
E = collect(eachcol(E_mat))

# plot band difference against QE reference
fig, ax, plt = get_bandplot(kpath, qe_bands.eigenvalues, E;
    kwargs1 = (label = "QE",),
    kwargs2 = (label = "Wannier.jl", linestyle = :dash),
    fermi_energy = εF,
)
fig

#=
## Fermi surface

Interpolate eigenvalues on a uniform ``30 \times 30 \times 30`` mesh.
`Wannier.fermisurf` handles the endpoint convention (bxsf needs the last kpoint
to be the periodic image of the first, so the actual grid is ``31^3``).
=#
Wannier.Tools.fermisurf(interpolation_model; nk = 30, ef = εF, outprefix = "Cu")

#=
The output `Cu.bxsf` can be visualised with e.g. FermiSurfer or VESTA.
=#
bxsf = read_bxsf("Cu.bxsf")

#=
!!! note

    To render the Fermi surface interactively, load a Makie backend and call
    the `bandplot` recipe on the bxsf grid, or use an external tool such as
    [FermiSurfer](https://fermisurfer.osdn.jp/).
=#

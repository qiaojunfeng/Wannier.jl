# # 8. Interpolation of Fermi surface

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
3. run `Wannier.jl` [`disentangle`](@ref) on the `Model` to minimize the spread
4. write the maximal localized gauge to a new `amn` file

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`tutorial.ipynb`](./tutorial.ipynb)
    - Julia script: [`tutorial.jl`](./tutorial.jl)
=#

# ## Preparation
# Load the package
using WannierIO
using Wannier
using WannierPlots

#=
!!! tip

    Use the `run.sh` script which automate the scf, nscf, pw2wannier90 steps.
=#

#=
## Model generation

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct a [`Model`](@ref) that abstracts the calculation
=#
model = read_w90("cu")

#=
!!! tip

    The [`read_w90`](@ref) function will parse the `win` file and set the frozen window for the `Model` according to
    the `dis_froz_min` and `dis_froz_max` parameters in the `win` file. However, you can also change these parameters
    by calling the [`set_frozen_win!`](@ref) function.
=#

#=
## Disentanglement and maximal localization

The [`disentangle`](@ref) function
will disentangle and maximally localize the spread
functional, and returns the gauge matrices `U`,
=#
U = disentangle(model);

# The initial spread is
omega(model)

# The final spread is
omega(model, U)

# save the new gauge to the model
model.U .= U;

#=
!!! note

    The convergence thresholds is determined by the
    keyword arguments of [`disentangle`](@ref), e.g., `f_tol` for the tolerance on spread,
    and `g_tol` for the tolerance on the norm of spread gradient, etc. You can use stricter thresholds
    to further minimize a bit the spread.
=#

# load QE band structure
kpoints_qe, E_qe = WannierIO.read_qe_band("qe_bands.dat");
# the Fermi energy from scf calculation
ef = 16.8985

#=
## Generate an [`InterpModel`](@ref)
=#
# Force using `kpoint_path` in `win` file
win = read_win("cu.win")
kpath = Wannier.get_kpath(win.unit_cell, win.kpoint_path)

interp_model = Wannier.InterpModel(model; kpath=kpath)

# interpolate band structure
kpi, E = Wannier.interpolate(interp_model)

# plot band difference
P = plot_band_diff(kpi, E_qe, E; fermi_energy=ef)
Main.HTMLPlot(P, 500)  # hide

#=
then interpolate the Fermi surface on a ``30 \times 30 \times 30`` mesh
=#
kpoints, E_fs = Wannier.fermi_surface(interp_model; n_k=30);

#=
save to a `bxsf` file
=#
# origin of the grid, always zeros
origin = zeros(Float64, 3)
WannierIO.write_bxsf("cu.bxsf", ef, origin, interp_model.recip_lattice, E_fs)

# show the Brillouin zone
using Brillouin
using PlotlyJS

# primitive reciprocal basis associated with k-path
bxsf = Wannier.read_bxsf("cu.bxsf")
fig = WannierPlots.plot_fermisurf_plotly(bxsf.rgrid, bxsf.fermi_energy, bxsf.E; kpath=kpath)
fig.layout.width = 500
fig.layout.height = 500
fig.layout.autosize = false
fig
Main.HTMLPlot(fig, 500)  # hide

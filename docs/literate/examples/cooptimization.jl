# # 10. Co-optimizing spin-up and spin-down WFs

#=
```@meta
CurrentModule = Wannier
```
=#

#=
Usually, for spin-polarized systems, we run two independent Wannierizations
for the spin-up and spin-down channels. However, in some cases, e.g., computing
magnetic exchange constants within an Heisenberg model, we need to make sure
the spin-up and spin-down WFs have as similar as possible real-space shapes,
and are centered on each atoms.

In this tutorial, we will show how to co-optimize the spin-up and spin-down
WFs, with two constraints:
1. WF center constraint to construct atom-centered WFs
2. spin-up and down overlap constraint to make sure each pair of
    spin-up and spin-down WFs are as similar as possible

We will Wannierize a 2D ``CrI_3`` system, using pseudo-atomic-orbital projections
as the starting guess (computed by QE).

## Outline

1. plot QE band structure as a reference
2. run two independent Wannierizations of spin-up and spin-down channels
3. construct a [`MagModel`](@ref) that merges the two spin channels
4. disentangle, with overlap constraint
5. disentangle, with both WF center and overlap constraints

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`tutorial.ipynb`](./tutorial.ipynb)
    - Julia script: [`tutorial.jl`](./tutorial.jl)
=#

# ## Preparation
# Load the packages
using WannierIO
using Wannier
using WannierPlots

#=
## Plot QE band structure

It's always good to check band-structure-interpolation qualities
of WFs against DFT bands.
To do this, let's first load QE band structure
=#
qe = WannierIO.read_qe_xml("qe_bands.xml")

#=
Note that I used Wannier90 generated kpath for QE bands calculation, to be
comparable with Wannier-interpolated bands.
To generate a kpath equivalent to Wannier90, here we use the [`get_kpath`](@ref)
function, with `unit_cell` and `kpoint_path` parsed from `win` file
by [`read_win`](@ref) function.
=#
win = read_win("up/cri3_up.win")
kpath = Wannier.get_kpath(win.unit_cell, win.kpoint_path)

#=
then we construct an [`KPathInterpolant`](@ref) object which stores the exact
kpoint coordinates to be interpolated, using 100 points in the 1st kpath segment
(equivalent to Wannier90 input `bands_num_points`).
=#
kpi = Wannier.interpolate_w90(kpath, 100)

# Now we can plot the QE bands for two spin channels
P = plot_band_diff(kpi, qe.E_up, qe.E_dn; fermi_energy=qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

#=
## Model construction

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct two [`Model`](@ref)s
for spin-up and spin-down channels.
Note the frozen windows for the two channels are set independently
according to the `dis_froz_max` inside the two `win` files, respectively;
in our case, they are both `dis_froz_max = -2 eV` since there are gaps
in both spin channels.
=#
model_up = read_w90("up/cri3_up")
model_dn = read_w90("dn/cri3_dn")

#=
## Projection-only WFs

The projection-only WFs, i.e., without any disentanglement or maximal
localization, are usually centered on each atom. This can be checked
by computing WF centers and spreads
=#
omega(model_up)

# and for spin-down channel
omega(model_dn)

#=
However, their band interpolations are often not good, so they shouldn't
be used for magnetic exchange constants calculations.

In our case, the projection-only WFs are `Cr:4s,3d` and `I:5s,5p` orbitals
centered on atoms. We will use Wannier interpolation for band structures,
and compare them with QE bands.

We first construct two [`InterpModel`](@ref)s for spin-up and spin-down channels,
using the `kpoint_path` from `win` file (otherwise, by default the `InterpModel`
will use `Brillouin.jl` to auto generate a kpath, which might be different
from user's input)
=#
interpModel_up = Wannier.InterpModel(model_up; kpath=kpath)
interpModel_dn = Wannier.InterpModel(model_dn; kpath=kpath)

# then interpolate eigenvalues
E_up_projonly = Wannier.interpolate(interpModel_up, kpi)
E_dn_projonly = Wannier.interpolate(interpModel_dn, kpi)

# and plot the spin-up bands compared with QE
P = plot_band_diff(kpi, qe.E_up, E_up_projonly; fermi_energy=qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

# and the spin-down bands
P = plot_band_diff(kpi, qe.E_dn, E_dn_projonly; fermi_energy=qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

#=
As can be seen from the above two figures, the projection-only WFs do not
reproduce DFT bands, i.e., they do not correctly describe the electronic
structure of the system, thus should not be used for physical property
calculations.

As a side node, we can also plot the Wannier-interpolated spin-up and down bands
in one figure,
=#
P = plot_band_diff(kpi, E_up_projonly, E_dn_projonly; fermi_energy=qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

#=
Save the MLWF gauge into `chk` files, to allow other codes, e.g., Wannier90,
to restart from the generated gauge.
Note we need to explicitly provide the `exclude_bands`, so as to be consistent
with `win` file.

!!! note

    The written `chk` are Fortran formatted files, which you can convert to
    binary file using Wannier90 `chk2chk.x` executable.
    There is an additional keyword argument `binary` for [`write_chk`](@ref),
    which is `false` by default, but can be set to `true` to write binary `chk`;
    however, Fortran binary format is compiler-dependent, so it's not gaurantted
    to be able to restart from the written binary `chk` file for other codes.
=#
exclude_bands = collect(1:8)
Wannier.write_chk("up/wjl_up_projonly.chk", model_up; exclude_bands)
Wannier.write_chk("dn/wjl_dn_projonly.chk", model_dn; exclude_bands)

#=
## Independent Wannierizations of two spin channels

Previous band structures show that projection-only WFs are not good enough.
Now, we will run disentanglement (and maximal localization) independently for
the two spin channels, to construct WFs that can accurately interpolate
band structures.
Note we just run 100 iterations here, this only converge to the order of
5e-2, but is already good enough for band interpolation.
=#
U_up_mlwf = disentangle(model_up);

# and WF centers and spreads
omega(model_up, U_up_mlwf)

# For band interpolations, we explicitly construct [`InterpModel`](@ref)s by
# reusing previous ``\bm{R}`` vectors, and compute Hamiltonian ``H(\bm{R})``.
# This skips ``\bm{R}`` vector generation, a bit faster
interpModel_up_mlwf = Wannier.InterpModel(
    interpModel_up.kRvectors,
    interpModel_up.kpath,
    fourier(interpModel_up.kRvectors, Wannier.get_Hk(model_up.E, U_up_mlwf)),
)
E_up_mlwf = Wannier.interpolate(interpModel_up_mlwf, kpi)

# Now the MLWF bands are very accurate, much better than projection-only
P = plot_band_diff(kpi, qe.E_up, E_up_mlwf; fermi_energy=qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

# similarly, for spin-down channel
U_dn_mlwf = disentangle(model_dn);
# and WF centers and spreads
omega(model_dn, U_dn_mlwf)
# and the interpolated bands
interpModel_dn_mlwf = Wannier.InterpModel(
    interpModel_dn.kRvectors,
    interpModel_dn.kpath,
    fourier(interpModel_dn.kRvectors, Wannier.get_Hk(model_dn.E, U_dn_mlwf)),
)
E_dn_mlwf = Wannier.interpolate(interpModel_dn_mlwf, kpi)
P = plot_band_diff(kpi, qe.E_dn, E_dn_mlwf; fermi_energy=qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

# Save the gauge into `chk` files
Wannier.write_chk("up/wjl_up_mlwf.chk", model_up, U_up_mlwf; exclude_bands)
Wannier.write_chk("dn/wjl_dn_mlwf.chk", model_dn, U_dn_mlwf; exclude_bands)

#=
Although the two sets of MLWFs accurately reproduce the band structures,
their centers and spreads are often different: in simple cases,
they might exchange orders; in more complex cases, they can hybridize
into different orbitals. This causes problems for magnetic exchange constants
calculations.
=#

#=
## Disentanglement with spin-up-down overlap constraint

In this section, we will impose spin-up and spin-down WF overlap constraint,
to disentangle simeultaneously the two sets of WFs, so that they can accurately
interpolate band structures, and more importantly, have the same centers and
spreads.

We first construct a [`MagModel`](@ref) for the two spin channels,
using the previous two `Model`s and an additional overlap matrix.

The spin-up and down overlap matrices is written in the same format as `amn`
=#
Mud = read_amn("updn/cri3_updn.mud");

# then assemble into a [`MagModel`](@ref)
model = Wannier.MagModel(model_up, model_dn, Mud)

#=
Now let's disentangle with spin overlap constraint.
Here `λs` is the Lagrange multiplier for the constraint.
=#
λs = 10.0
U_up, U_dn = disentangle(model, λs);
#=
The resulting spin-up and spin-down WFs have very similar centers and spreads,
however, their centers drift from the original positions which were centered
on atoms.
=#
omega(model, U_up, U_dn, λs)

# as a comparison, the spreads of independent Wannierizations are
omega(model, U_up_mlwf, U_dn_mlwf, λs)

# Save the gauge into `chk` files
Wannier.write_chk("up/wjl_up_cowf.chk", model_up, U_up; exclude_bands)
Wannier.write_chk("dn/wjl_dn_cowf.chk", model_dn, U_dn; exclude_bands)

#=
As can be seen from the above spin-up-down overlaps `<↑|↓>`,
the two sets of WFs are sort of randomly distributed, some overlaps are 0,
thus the two Hamiltonians are not using the same real-space basis functions.
=#

#=
## Disentanglement with both WF center and overlap constraints

Now, we impose both WF center constraints and spin-up-down WF overlap
constraints, to disentangle simeultaneously the two sets of WFs.
Thus, our final WFs can accurately interpolate band structures,
have the same WF centers and spreads, and are atom-centered.
These WFs satisfy the assumptions for magnetic exchange constants
calculations, i.e., constructing a Heisenberg model with
1. atom-centered orbitals with localized magnetic moments;
2. the spin-up and spin-down orbitals can map one-by-one to each other
    so they have the same real-space basis;
3. both spin channels accurately describe the electronic structure of
    the system.

Since we want atom-centered WFs, our target centers are just atom positions.
We store our target WF centers in a column-wise matrix,
=#
r₀ = zeros(3, model.up.n_wann)
# the first 6 WFs are the `4s,3d` orbitals of the 1st `Cr` atom
r₀[:, 1:6] .= model_up.atom_positions[:, 1]
# the next 6 WFs are the `4s,3d` orbitals of the 2nd `Cr` atom
r₀[:, 7:12] .= model_up.atom_positions[:, 2]
# the next 4 WFs are the `5s,5p` orbitals of `I` atom
r₀[:, 13:16] .= model_up.atom_positions[:, 3]
# and similarly for the remaining 5 `I` atoms
r₀[:, 17:20] .= model_up.atom_positions[:, 4]
r₀[:, 21:24] .= model_up.atom_positions[:, 5]
r₀[:, 25:28] .= model_up.atom_positions[:, 6]
r₀[:, 29:32] .= model_up.atom_positions[:, 7]
r₀[:, 33:36] .= model_up.atom_positions[:, 8]

# convert to Cartesian coordinates
r₀ = model.up.lattice * r₀

#=
Now we need to choose the Lagrange multiplier factors for the two constraints.
Here the `λc` is the Lagrange multiplier for the WF center constraint,
and `λs` is for the spin overlap constraint.

We use `10.0` for both `λc` and `λs` here, but you can try different values.

On output, the second last column, ωc, is the penalty of the WF center
constraint, i.e., the larger the ωc, the farther the WF center is from the
target position r₀.
=#
λc = 10.0
λs = 10.0
U_up, U_dn = disentangle_center(model, r₀, λc, λs);

# the final centers and spreads are
omega(model, U_up, U_dn, r₀, λc, λs)

# as a comparison, the spreads of independent Wannierizations are
omega(model, U_up_mlwf, U_dn_mlwf, r₀, λc, λs)

#=
Notice the 2nd last column of ωc, which is the penalty for WF centers:
in the case of independent Wannierizations, some WF centers are far away
from atom positions, so we don't have a atom-centered Hamiltonian model.
Moreover, as mentioned above, some spin overlaps are 0, so the spin up and
down have different real-space basis functions.

In comparison, the two constraints work properly: the resulting WFs
are well-centered on atoms (all penalties < 4e-3), and each spin-up WF
map one-by-one to each spin-down WF (all `<↑|↓>` > 0.98).

Finally, we check the interpolated bands.
=#
interpModel_up = Wannier.InterpModel(
    interpModel_up.kRvectors,
    interpModel_up.kpath,
    fourier(interpModel_up.kRvectors, Wannier.get_Hk(model_up.E, U_up)),
)
interpModel_dn = Wannier.InterpModel(
    interpModel_dn.kRvectors,
    interpModel_dn.kpath,
    fourier(interpModel_dn.kRvectors, Wannier.get_Hk(model_dn.E, U_dn)),
)
E_up = Wannier.interpolate(interpModel_up, kpi);
E_dn = Wannier.interpolate(interpModel_dn, kpi);

# and compare the up bands against QE
P = plot_band_diff(kpi, qe.E_up, E_up; fermi_energy=qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

# and the down bands
P = plot_band_diff(kpi, qe.E_dn, E_dn; fermi_energy=qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

# Save the gauge into `chk` files
Wannier.write_chk("up/wjl_up_cowf_c.chk", model_up, U_up; exclude_bands)
Wannier.write_chk("dn/wjl_dn_cowf_c.chk", model_dn, U_dn; exclude_bands)

#=
From the above two figures, we know that our WFs accurately reproduce the
``CrI_3`` Hamiltonian.

Finally, we have two sets of WFs that are atom-centered,
have similar centers and spreads,
and can accurately interpolate band structures!
We can now pass them to magnetic exchange constant calculations.
=#

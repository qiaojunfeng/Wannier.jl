# # 10. Collinear spin system

#=
```@meta
CurrentModule = Wannier
```
=#

#=
Usually, for spin-polarized systems, we run two independent Wannierizations
for the spin-up and spin-down channels.

In this tutorial, we will show how to co-optimize the spin-up and spin-down
WFs, with two constraints:

We will Wannierize a iron crystal, using pseudo-atomic-orbital projections
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
using Wannier.Datasets
# using WannierPlots

#=
## Plot QE band structure

It's always good to check band-structure-interpolation qualities
of WFs against DFT bands.
To do this, let's first load QE band structure
=#
qe = WannierIO.read_qe_xml(dataset"Fe_collinear/reference/qe_bands.xml")

#=
Note that I used Wannier90 generated kpath for QE bands calculation, to be
comparable with Wannier-interpolated bands.
To generate a kpath equivalent to Wannier90, here we use the [`get_kpath`](@ref)
function, with `unit_cell` and `kpoint_path` parsed from `win` file
by [`read_win`](@ref) function.
=#
win = read_win(dataset"Fe_collinear/Fe_up.win")
kpath = Wannier.generate_kpath(win.unit_cell_cart, win.kpoint_path)

#=
then we construct an [`KPathInterpolant`](@ref) object which stores the exact
kpoint coordinates to be interpolated, using 100 points in the 1st kpath segment
(equivalent to Wannier90 input `bands_num_points`).
=#
kpi = Wannier.generate_w90_kpoint_path(kpath)

using PlotlyJS
# Now we can plot the QE bands for two spin channels
P = plot_band_diff(kpi, qe.eigenvalues_up, qe.eigenvalues_dn; fermi_energy=qe.fermi_energy)
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
model_up = read_w90(dataset"Fe_collinear/Fe_up")
model_dn = read_w90(dataset"Fe_collinear/Fe_dn")

using LinearAlgebra
model = Wannier.Model(
    model_up.lattice,
    model_up.atom_positions,
    model_up.atom_labels,
    model_up.kstencil,
    map(model_up.overlaps, model_dn.overlaps) do Mb_up, Mb_dn
        map(Mb_up, Mb_dn) do M_up, M_dn
            [
                M_up zeros(size(M_up, 1), size(M_dn, 2))
                zeros(size(M_dn, 1), size(M_up, 2)) M_dn
            ]
        end
    end,
    map(model_up.gauges, model_dn.gauges) do U_up, U_dn
        [
            U_up zeros(size(U_up, 1), size(U_dn, 2))
            zeros(size(U_dn, 1), size(U_up, 2)) U_dn
        ]
    end,
    map(model_up.eigenvalues, model_dn.eigenvalues) do E_up, E_dn
        vcat(E_up, E_dn)
    end,
    map(model_up.frozen_bands, model_dn.frozen_bands) do f_up, f_dn
        vcat(f_up, f_dn)
    end,
    map(model_up.entangled_bands, model_dn.entangled_bands) do e_up, e_dn
        vcat(e_up, e_dn)
    end,
)

#=
## Projection-only WFs

The projection-only WFs, i.e., without any disentanglement or maximal
localization, are usually centered on each atom. This can be checked
by computing WF centers and spreads
=#
omega(model_up)

# and for spin-down channel
omega(model_dn)

omega(model)

U = disentangle(model);

omega(model, U)

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
Ham_0 = Wannier.TBHamiltonian(model)
Ham_1 = Wannier.TBHamiltonian(model, U)

interp_0 = Wannier.HamiltonianInterpolator(Ham_0)
interp_1 = Wannier.HamiltonianInterpolator(Ham_1)

bands_0 = interp_0(kpi)[1];
bands_1 = interp_1(kpi)[1];

# and compare the up bands against QE
bands_qe = map(qe.eigenvalues_up, qe.eigenvalues_dn) do E_up, E_dn
    vcat(E_up, E_dn)
end
P = plot_band_diff(kpi, bands_qe, bands_0; fermi_energy=qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

# and the down bands
P = plot_band_diff(kpi, bands_qe, bands_1; fermi_energy=qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

s_z = map(1:n_kpoints(model)) do ik
    [
        diagm(ones(n_bands(model_up))) zeros(n_bands(model_up), n_bands(model_dn))
        zeros(n_bands(model_dn), n_bands(model_up)) diagm(-ones(n_bands(model_dn)))
    ]
end
s_x = map(1:length(s_z)) do ik
    zeros(size(s_z[1]))
end
s_y = map(1:length(s_z)) do ik
    zeros(size(s_z[1]))
end
spin_vecs = map(zip(s_x, s_y, s_z)) do (x, y, z)
    Wannier.MVec3.(x, y, z)
end

Rspace = generate_Rspace(model)
spin = Wannier.TBSpin(Rspace, model.kpoints, spin_vecs, U)

θ = 0
ϕ = 0
spin_interp = Wannier.SpinProjectionInterpolator(Ham_1, spin, θ, ϕ)
S_band = spin_interp(kpi);

# Save the gauge into `chk` files
Wannier.plot_band(kpi, bands_1; fermi_energy=qe.fermi_energy, color=S_band)

#=
From the above two figures, we know that our WFs accurately reproduce the
``CrI_3`` Hamiltonian.

Finally, we have two sets of WFs that are atom-centered,
have similar centers and spreads,
and can accurately interpolate band structures!
We can now pass them to magnetic exchange constant calculations.
=#

# Use adaptive kmesh refinement to converge Fermi energy
using Printf
using LinearAlgebra
using Wannier
using Wannier.Datasets

prefix = dataset"Fe_soc/reference/MDRS/Fe"
tb = read_w90_tb(prefix)
interp = HamiltonianInterpolator(tb.hamiltonian)

lattice = real_lattice(interp)
recip_lattice = reciprocal_lattice(lattice)
kdistance = 0.08
kgrid = round.(Int, [norm(b) for b in eachcol(recip_lattice)] ./ kdistance)

kbT = 0.01  # eV
# For NoneSmearing, usually need a larger tolerance, otherwise cannot bisect Fermi energy
# smearing = Wannier.NoneSmearing()
# tol_n_electrons = 1e-4
# εF, adpt_kgrid = Wannier.compute_fermi_energy(kgrid, interp, num_electrons, kbT, smearing, tol_n_electrons)

n_electrons = 8
# This is a TB Hamiltonian from a SOC calculation
prefactor = 1

# smearing = Wannier.FermiDiracSmearing()
smearing = Wannier.ColdSmearing()
tol_εF = 5e-3  # convergence tolerance for Fermi energy in eV

# You can directly compute Fermi energy by
εF = Wannier.compute_fermi_energy(
    kgrid, interp, n_electrons, kbT, smearing; prefactor, tol_εF
)

# But here we will separate it into two steps, 1st construct the adaptive kgrid,
# which is actually a uniformlly-spaced grid on input, but the [`AdaptiveKgrid`](@ref)
# data structure allows it to be iteratively refined on selected kpoints later on.
adpt_kgrid = Wannier.AdaptiveKgrid(kgrid, interp)
# then compute Fermi energy
εF = Wannier.compute_fermi_energy!(
    adpt_kgrid, interp, n_electrons, kbT, smearing; prefactor, tol_εF
)
# The `Wannier.compute_fermi_energy!` function is a lower-level function,
# will modify the input `adpt_kgrid`, and from `adpt_kgrid` we can retrieve
# all the interpolated eigenvalues that are used to compute Fermi energy.
# For example, the fractional coordinates of the 1st kpoint is
adpt_kgrid.kvoxels[1].point
# and the interpolated eigenvalues are
adpt_kgrid.vals[1]

# actual output is: `Fermi Energy after interpolation: 17.60847657 eV`
@printf("Fermi Energy after interpolation: %.8f eV\n", εF)

# As a comparison, the Fermi energy from pw.x scf calculation is 17.6132 eV,
# which is computed on a 8x8x8 kgrid with 0.02 Ry cold smearing.

using ProgressMeter
# Compute n_electrons w.r.t. a range of Fermi energy
εF_range = range(εF - 1, εF + 1; step=10e-3)
n_electrons_range = @showprogress map(εF_range) do εi
    occs = Wannier.occupation(adpt_kgrid, εi, kbT, smearing; prefactor)
    kweights = Wannier.default_kweights(adpt_kgrid)
    Wannier.compute_n_electrons(occs, kweights)
end

# find VBM, CBM
vbm = Wannier.find_vbm(adpt_kgrid.vals, εF)[1]
cbm = Wannier.find_cbm(adpt_kgrid.vals, εF)[1]

using Plots
plot(εF_range, n_electrons_range)
xlabel!("Fermi energy (eV)")
ylabel!("n_electrons")
vline!([εF]; label="Fermi energy")
vline!([vbm]; label="VBM")
vline!([cbm]; label="CBM")
xlims!(εF - 0.5, εF + 0.5)
# gui()

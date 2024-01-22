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

εF_scf = read_win(dataset"Fe_soc/Fe.win").fermi_energy

# You can directly compute Fermi energy by
εF = Wannier.compute_fermi_energy(
    kgrid, interp, n_electrons, kbT, smearing; prefactor, tol_εF
)

# But here we will separate it into two steps, 1st construct the adaptive kgrid,
# which is actually a uniformlly-spaced grid on input, but the [`AdaptiveKgrid`](@ref)
# data structure allows it to be iteratively refined on selected kpoints later on.
kpoints = get_kpoints(kgrid)
eigenvalues = interp(kpoints)[1]
adpt_kgrid = Wannier.AdaptiveKgrid(kpoints, eigenvalues)

# note one can also store the initial uniform grid to a `bxsf` file
"""
Add replica of the 1st kpoint to the last one, so the kgrid is periodic and
can be saved to a `bxsf` file.

!!! warning

    This function assumes the kpoints of the input `adpt_kgrid` are constructed
    by `get_kpoints(kgrid_size; endpoint=false)`, i.e.,
    - uniformly spaced
    - do not contain kpoints having coordinate equal 1
    - the z coordinate increases the fastest
"""
function to_periodic_vals(adpt_kgrid::Wannier.AdaptiveKgrid)
    nks = Wannier.guess_kgrid_size([kv.point for kv in adpt_kgrid.kvoxels])
    nbands = length(adpt_kgrid.vals[1])
    E = zeros(Float64, nbands, (nks .+ 1)...)
    counter = 1
    for i in 1:nks[1]
        for j in 1:nks[2]
            for k in 1:nks[3]
                E[:, i, j, k] .= adpt_kgrid.vals[counter]
                counter += 1
            end
        end
    end
    E[:, :, :, end] .= E[:, :, :, 1]
    E[:, :, end, :] .= E[:, :, 1, :]
    E[:, end, :, :] .= E[:, 1, :, :]
    return E
end
origin = zeros(3)
E = to_periodic_vals(adpt_kgrid)
WannierIO.write_bxsf("Fe_soc.bxsf", εF, origin, recip_lattice, E)
# or construct the initial uniform grid from a `bxsf` file
"""
Remove the replica of the last kpoint, and construct an `AdaptiveKgrid` that
is uniformlly-spaced.

!!! warning

    This function assumes the input eigenvalues `E` are read from a `bxsf` file.
"""
function to_adpt_kgrid(E::Array{Float64, 4})
    nks = collect(size(E)[2:end] .- 1)
    kpoints = Wannier.get_kpoints(nks)
    nbands = size(E, 1)
    eigenvalues = [zeros(nbands) for _ in 1:length(kpoints)]
    counter = 1
    for i in 1:nks[1]
        for j in 1:nks[2]
            for k in 1:nks[3]
                eigenvalues[counter] .= E[:, i, j, k]
                counter += 1
            end
        end
    end
    return Wannier.AdaptiveKgrid(kpoints, eigenvalues)
end
E = WannierIO.read_bxsf("Fe_soc.bxsf").E
adpt_kgrid_bxsf = to_adpt_kgrid(E)

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

import Base.Broadcast.broadcastable
using Roots
using SpecialFunctions: erf

abstract type SmearingFunction end

Base.Broadcast.broadcastable(S::SmearingFunction) = Ref(S)

struct FermiDiracSmearing <: SmearingFunction end

occupation(x, ::FermiDiracSmearing) = 1 / (1 + exp(x))

struct ColdSmearing <: SmearingFunction end

function occupation(x::T, ::ColdSmearing) where {T}
    return (
        -erf(x + 1 / sqrt(T(2))) / 2 +
        1 / sqrt(2 * T(π)) * exp(-(-x - 1 / sqrt(T(2)))^2) +
        1 / T(2)
    )
end

struct NoneSmearing <: SmearingFunction end

occupation(x, ::NoneSmearing) = x > 0 ? zero(x) : one(x)

default_occupation_prefactor() = 2

@doc raw"""
Compute occupation given eigenvalues and Fermi energy.

# Arguments
- `eigenvalues`: eigenvalues in eV
- `εF`: Fermi energy in eV

# Keyword arguments
- `kBT`: temperature in the same unit as `E`, i.e., $k_B T$ in eV
- `prefactor`: 1 for collinear calculation, 2 for spinless
"""
function occupation(
    eigenvalues::AbstractVector,
    εF::Real,
    kBT::Real,
    smearing::SmearingFunction;
    prefactor::Real=default_occupation_prefactor(),
)
    T = promote_type(eltype(eltype(eigenvalues)), typeof(εF), typeof(kBT))
    inv_kBT = iszero(kBT) ? T(Inf) : 1 / kBT

    occ = map(eigenvalues) do εk
        prefactor * occupation.((εk .- εF) .* inv_kBT, smearing)
    end
    return occ
end

"""
    $(SIGNATURES)

Compute the number of electrons with a given density of states.

# Arguments
- `energy`: Vector of energy, in eV unit
- `dos`: density of states on energy grid, in states/eV unit
- `εF`: Fermi energy, in eV unit

# Return
- `n_electrons`: number of electrons
"""
function compute_n_electrons(energy::AbstractVector, dos::AbstractVector, εF::Number)
    dE = energy[2] - energy[1]
    cum_dos = cumsum(dos) * dE

    idx = argmin(abs.(energy .- εF))
    tol_energy = 1e-5
    if abs(energy[idx] - εF) > tol_energy
        error("Fermi energy not found in energy grid")
    end

    return cum_dos[idx]
end

default_kweights(eigenvalues::AbstractVector) = 1 / length(eigenvalues)

function compute_n_electrons(
    occupation::AbstractVector, kweights=default_kweights(occupation)
)
    return sum(kweights .* sum.(occupation))
end

function compute_fermi_energy(
    eigenvalues::AbstractVector,
    n_electrons::Real,
    kBT::Real,
    smearing::SmearingFunction;
    prefactor::Real=default_occupation_prefactor(),
    kweights=default_kweights(eigenvalues),
    tol_n_electrons::Real=1e-6,
)
    # Get rough bounds to bracket εF
    min_ε = minimum(minimum, eigenvalues) - 1
    max_ε = maximum(maximum, eigenvalues) + 1

    excess(εF) = begin
        occ = occupation(eigenvalues, εF, kBT, smearing; prefactor)
        compute_n_electrons(occ, kweights) - n_electrons
    end
    @assert excess(min_ε) <= 0 <= excess(max_ε) "Fermi energy not bracketed $(excess(min_ε)) $(excess(max_ε))"

    εF = Roots.find_zero(excess, (min_ε, max_ε), Roots.Bisection(); atol=tol_n_electrons)
    Δn_elec = excess(εF)
    abs(Δn_elec) > tol_n_electrons &&
        error("Failed to find Fermi energy within tolerance, Δn_elec = $Δn_elec")

    return εF
end

struct Kvoxel{T,VT<:AbstractVector{T}}
    """fractional coordinates of kpoint"""
    point::VT

    """length of the kvoxel along three dimensions"""
    dv::VT

    """weight of the kpoint"""
    weight::T
end

struct AdaptiveKgrid{KV<:Kvoxel,VT}
    kvoxels::Vector{KV}
    vals::Vector{VT}
end

Base.length(ag::AdaptiveKgrid) = length(ag.kvoxels)

function occupation(
    adpt_grid::AdaptiveKgrid,
    εF::Real,
    kBT::Real,
    smearing::SmearingFunction;
    prefactor::Real=default_occupation_prefactor(),
)
    T = promote_type(eltype(adpt_grid.vals), typeof(εF), typeof(kBT))
    inv_kBT = iszero(kBT) ? T(Inf) : 1 / kBT

    occ = map(adpt_grid.vals) do εk
        prefactor * occupation.((εk .- εF) .* inv_kBT, smearing)
    end
    return occ
end

default_kweights(adpt_grid::AdaptiveKgrid) = [kv.weight for kv in adpt_grid.kvoxels]

"""
Refine the kgrid by splitting the kvoxels into subvoxels.

# Arguments
- `ag`: `AdaptiveKgrid`
- `iks`: indices of kvoxels to be refined

# Keyword arguments
- `n_subvoxels`: number of subvoxels along each dimension. 2 -> split into 8 subvoxels
"""
function refine!(ag::AdaptiveKgrid, iks::AbstractVector, interp::Function; n_subvoxels=2)
    new_kvoxels = eltype(ag.kvoxels)[]

    # split the current kvoxel into 8 sub kvoxels, so 7 new kvoxels are added
    range_subs = 0:(n_subvoxels - 1)
    add_points = [Vec3(i, j, k) for i in range_subs for j in range_subs for k in range_subs]
    deleteat!(add_points, 1)

    for ik in iks
        # split the current kvoxel into 8 sub kvoxels
        vx0 = ag.kvoxels[ik]
        voxel = Kvoxel(vx0.point, vx0.dv ./ n_subvoxels, vx0.weight / n_subvoxels^3)
        ag.kvoxels[ik] = voxel
        sub_voxels = map(add_points) do pt
            Kvoxel(voxel.point + pt .* voxel.dv, voxel.dv, voxel.weight)
        end
        append!(new_kvoxels, sub_voxels)
    end

    new_vals = interp([v.point for v in new_kvoxels])
    append!(ag.kvoxels, new_kvoxels)
    append!(ag.vals, new_vals)
    return nothing
end

function AdaptiveKgrid(kgrid::AbstractVector, interp::HamiltonianInterpolator)
    kpoints = get_kpoints(kgrid)
    eigenvals, _ = interp(kpoints)

    dv = Vec3(1 ./ kgrid)
    kweight = default_kweights(eigenvals)
    kvoxels = map(kpoints) do kpt
        Kvoxel(kpt, dv, kweight)
    end
    return AdaptiveKgrid(kvoxels, eigenvals)
end

"""
Compute Fermi energy by recursively refining the kgrid when interpolating the Hamiltonian.

# Return
- `εF`: Fermi energy
"""
function compute_fermi_energy(
    kgrid::AbstractVector,
    interp::HamiltonianInterpolator,
    n_electrons::Real,
    kBT::Real,
    smearing::SmearingFunction;
    kwargs...,
)
    adpt_kgrid = AdaptiveKgrid(kgrid, interp)
    return compute_fermi_energy!(adpt_kgrid, interp, n_electrons, kBT, smearing; kwargs...)
end

"""
Compute Fermi energy by recursively refining the kgrid when interpolating the Hamiltonian.

On output, the input variable `adpt_kgrid` contains the refined kgrid and eigenvalues.

# Return
- `εF`: Fermi energy
"""
function compute_fermi_energy!(
    adpt_kgrid::AdaptiveKgrid,
    interp::HamiltonianInterpolator,
    n_electrons::Real,
    kBT::Real,
    smearing::SmearingFunction;
    prefactor::Real=default_occupation_prefactor(),
    tol_n_electrons::Real=1e-6,
    tol_εF::Real=5e-3,
    max_refine::Integer=10,
)
    # the initial guessing Fermi energy
    εF = compute_fermi_energy(
        adpt_kgrid.vals,
        n_electrons,
        kBT,
        smearing;
        prefactor,
        kweights=default_kweights(adpt_kgrid),
        tol_n_electrons,
    )
    @printf("εF on input kgrid   : %15.9f eV, n_kpoints = %8d\n", εF, length(kpoints))

    εF_prev = εF - 1
    iter = 1
    # search range
    width_εF = 0.5
    while abs(εF - εF_prev) > tol_εF && iter <= max_refine
        refine_iks = filter(1:length(adpt_kgrid)) do ik
            any(abs.(adpt_kgrid.vals[ik] .- εF) .<= width_εF)
        end
        # alternate between even and odd refinement, so it works for the
        # K/K' point of graphene as well
        # I should iterate odd grid 1st, otherwise it seems the graphene
        # case could still stuck at wrong εF with [8, 8, 1] kgrid
        n_subvoxels = iter % 2 == 0 ? 2 : 3
        refine!(adpt_kgrid, refine_iks, x -> interp(x)[1]; n_subvoxels)

        εF_prev = εF
        εF = compute_fermi_energy(
            adpt_kgrid.vals,
            n_electrons,
            kBT,
            smearing;
            prefactor,
            kweights=default_kweights(adpt_kgrid),
            tol_n_electrons,
        )
        # gradually reduce width_εF to save computation
        ΔεF = εF - εF_prev
        @printf(
            "εF at iteration %3d : %15.9f eV, n_kpoints = %8d, ΔεF = %16.9e eV\n",
            iter,
            εF,
            length(adpt_kgrid),
            ΔεF,
        )
        iter += 1
        # after 10 iters, the width is mutiplied by 0.8^10 ≈ 0.107
        # width_εF *= 0.8
        # set next search range according to ΔεF
        width_εF = min(width_εF, abs(ΔεF) * 5)
    end
    return εF
end

"""
Find the valence band maximum.

# Arguments
- `eigenvalues`: eigenvalues
- `εF`: Fermi energy

# Return
- `vbm`: valence band maximum
- `ik`: index of kpoint for the vbm
- `n`: index of band for the vbm
"""
function find_vbm(eigenvalues::AbstractVector, εF::Real)
    # convert to a n_bands x n_kpoints matrix
    E = reduce(hcat, eigenvalues)
    # mask the conduction bands
    E[E .> εF] .= -Inf
    n, ik = argmax(E).I
    vbm = E[n, ik]
    return vbm, ik, n
end

"""
Find the conduction band minimum.

# Arguments
- `eigenvalues`: eigenvalues
- `εF`: Fermi energy

# Return
- `cbm`: conduction band minimum
- `ik`: index of kpoint for the cbm
- `n`: index of band for the cbm
"""
function find_cbm(eigenvals::AbstractVector, εF::Real)
    # convert to a n_bands x n_kpoints matrix
    E = reduce(hcat, eigenvals)
    # mask the valence bands
    E[E .< εF] .= Inf
    n, ik = argmin(E).I
    cbm = E[n, ik]
    return cbm, ik, n
end

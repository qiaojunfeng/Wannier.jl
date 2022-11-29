export interpolate

"""
    get_Hk(E, A)

Construct k space Hamiltonian Hᵏ.

```math
H_{\\bm{k}} = A_{\\bm{k}}^\\dagger [\\epsilon_{n \\bm{k}}] A_{\\bm{k}},
```
where ``[\\epsilon_{n \\bm{k}}]`` is a diagonal matrix with
``\\epsilon_{n \\bm{k}}`` as the diagonal elements.

# Arguments
- `E`: `n_bands * n_kpts`, energy eigenvalue
- `A`: `n_bands * n_wann * n_kpts`, gauge matrices
"""
function get_Hk(E::Matrix{T}, A::Array{U,3}) where {T<:Number,U<:Number}
    n_bands, n_wann, n_kpts = size(A)
    size(E) != (n_bands, n_kpts) && error("size(E) != (n_bands, n_kpts)")

    Hᵏ = zeros(U, n_wann, n_wann, n_kpts)
    for ik in 1:n_kpts
        Hᵏ[:, :, ik] = A[:, :, ik]' * Diagonal(E[:, ik]) * A[:, :, ik]
    end

    return Hᵏ
end

"""
    interpolate(model::InterpModel{T}, kpoints::Matrix{T}) where {T<:Real}

Interpolate energy eigenvalues at `kpoints`.

# Arguments
- `model`: `InterpModel`
- `kpoints`: `3 * n_kpts`, kpoints to be interpolated, in fractional coordinates,
    can be nonuniform.
"""
function interpolate(model::InterpModel{T}, kpoints::Matrix{T}) where {T<:Real}
    Hᵏ = invfourier(model.kRvectors, model.H, kpoints)
    # diagonalize
    Eᵏ, _ = diag_Hk(Hᵏ)

    return Eᵏ
end

"""
    diag_Hk(Hᵏ:: AbstractArray{T, 3}) where {T<:Complex}

Diagonalize k space Hamiltonian `H`.

# Arguments
- `H`: `n_wann * n_wann * n_kpts`, k space Hamiltonian

# Return
- `E`: `n_wann * n_kpts`, energy eigen values
- `V`: `n_wann * n_wann * n_kpts`, `V[:, i, ik]` is the i-th eigen vector at `ik`-th kpoint
"""
function diag_Hk(H::AbstractArray{T,3}) where {T<:Complex}
    n_wann = size(H, 1)
    n_kpts = size(H, 3)
    size(H) == (n_wann, n_wann, n_kpts) || error("size(H) != (n_wann, n_wann, n_kpts)")

    E = zeros(real(T), n_wann, n_kpts)
    V = similar(H)
    for ik in axes(H, 3)
        Hᵏ = @view H[:, :, ik]
        # @assert ishermitian(Hᵏ) norm(Hᵏ - Hᵏ')
        @assert norm(Hᵏ - Hᵏ') < 1e-10
        # Hᵏ = 0.5 * (Hᵏ + Hᵏ')
        ϵ, v = eigen(Hᵏ)
        E[:, ik] = real.(ϵ)
        H[:, :, ik] = v
    end

    return E, H
end

"""
    interpolate(model::InterpModel, kpi::KPathInterpolant)

Interpolate band structure along the given kpath.

# Arguments
- `model`: `InterpModel`
- `kpi`: `KPathInterpolant`
"""
function interpolate(model::InterpModel, kpi::KPathInterpolant)
    kpoints = get_kpoints(kpi)
    return interpolate(model, kpoints)
end

"""
    interpolate(model::InterpModel)

Interpolate band structure along kpath.

The `model.kpath` will be used.

# Arguments
- `model`: `InterpModel`

!!! note

    The kpath has the same density as `Wannier90`'s default, i.e.,
    `kpath_num_points = 100`.
"""
function interpolate(model::InterpModel)
    kpi = interpolate_w90(model.kpath, 100)
    return kpi, interpolate(model, kpi)
end

"""
    sort_hamiltonian_by_norm(Rvectors::RVectorsMDRS{T}, Hᴿ::Array{Complex{T},3}) where {T<:Real}

Sort the Hamiltonian ``H(\bm{R})`` by its norm in descending order.

# Arguments
- `Rvectors`: `RVectorsMDRS`
- `Hᴿ`: `n_wann * n_wann * n_rvecs`, Hamiltonian in real space.

# Return
- `idx`: the index to sort the Hamiltonian, i.e. `Hᴿ[:, :, idx]` is the sorted Hamiltonian.
- `normR`: the norm of R vectors, `normR[idx]` is the sorted `R` vectors in descending order.
- `normH`: the norm of H, `normH[idx]` is the sorted `Hᴿ` in descending order.

# Example
```julia
R, H, pos = read_w90_tb("mos2")
# Expand the H(R) to R̃, I assume the H(R) is read from w90 tb.dat file
H1 = Wannier.mdrs_v1tov2(R, H)
idx, normR, normH = sort_hamiltonian_by_norm(R, H1)
```
"""
function sort_hamiltonian_by_norm(
    Rvectors::RVectorsMDRS{T}, Hᴿ::Array{Complex{T},3}
) where {T<:Real}
    @assert size(Hᴿ, 3) == Rvectors.n_r̃vecs
    # Use the mdrsv2 R̃ vectors
    R = Rvectors.R̃vectors.R

    normH = [norm(Hᴿ[:, :, i]) for i in axes(Hᴿ, 3)]
    R_cart = Rvectors.lattice * R
    normR = [norm(i) for i in eachcol(R_cart)]

    idx = sortperm(normR)

    @printf("# idx    ||R||    ||H(R)||\n")
    for i in 1:20
        @printf("%2d    %.4f    %.4f\n", i - 1, normR[idx[i]], normH[idx[i]])
    end

    return idx, normR, normH
end

"""
    cut_hamiltonian(Rvectors::RVectorsMDRS{T}, Hᴿ::Array{Complex{T},3}, Rcut::T) where {T<:Real}

Cut real space Hamiltonian `Hᴿ` by the given cutoff radius.

# Arguments
- `Rvectors`: `RVectorsMDRS`
- `Hᴿ`: `n_wann * n_wann * n_rvecs`, Hamiltonian in real space.
- `Rcut`: cutoff radius in angstrom, the Hamiltonian `H(R)` corresponds to the R vectors with norm larger than `Rcut` will be set to 0.

# Return
- `H`: `n_wann * n_wann * n_rvecs`, cutted Hamiltonian in real space.

# Example
```julia
R, H, pos = read_w90_tb("mos2")
# Expand the H(R) to R̃, I assume the H(R) is read from w90 tb.dat file
H1 = Wannier.mdrs_v1tov2(R, H)
H2 = cut_hamiltonian(R, H1, 7.0)
```
"""
function cut_hamiltonian(
    Rvectors::RVectorsMDRS{T}, Hᴿ::Array{Complex{T},3}, Rcut::T
) where {T<:Real}
    @assert size(Hᴿ, 3) == Rvectors.n_r̃vecs

    idx, normR, normH = sort_hamiltonian_by_norm(Rvectors, Hᴿ)
    idx_cut = normR .> Rcut

    H = deepcopy(Hᴿ)
    H[:, :, idx_cut] .= 0

    return H
end
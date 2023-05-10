export interpolate

"""
    get_Hk(E, U)

Construct k space Hamiltonian Hᵏ.

```math
H_{\\bm{k}} = U_{\\bm{k}}^\\dagger [\\epsilon_{n \\bm{k}}] U_{\\bm{k}},
```
where ``[\\epsilon_{n \\bm{k}}]`` is a diagonal matrix with
``\\epsilon_{n \\bm{k}}`` as the diagonal elements.

# Arguments
- `E`: `n_bands * n_kpts`, energy eigenvalue
- `U`: `n_bands * n_wann * n_kpts`, gauge matrices
"""
function get_Hk(E::Matrix{T}, U::Array{S,3}) where {T<:Number,S<:Number}
    n_bands, n_wann, n_kpts = size(U)
    size(E) != (n_bands, n_kpts) && error("size(E) != (n_bands, n_kpts)")

    Hᵏ = [zeros(S, n_wann, n_wann) for k = 1:n_kpts]
    Threads.@threads for ik in 1:n_kpts
        # I need to force Hermiticity here, otherwise in some cases,
        # especially degenerate eigenvalues, the eigenvectors of Hᵏ,
        #   F = eigen(Hᵏ)
        # does not satisfy unitarity,
        #   F.vectors' ≈ F.vectors
        # and this leads to
        #   norm(F.vectors * diagm(F.values) * F.vectors' - Hᵏ) ≈ 1e-1
        # If I compute explicitly its inverse,
        #   norm(F.vectors * diagm(F.values) * inv(F.vectors) - Hᵏ) ≈ 1e-14
        # However, replacing all the `'` by `inv` is not a good idea,
        # since gauge rotation is used a lot throughout the code;
        # so I enforce Hermiticity here.
        # See also
        # https://discourse.julialang.org/t/a-b-a-is-not-hermitian-even-when-b-is/70611
        Uk = @view U[:, :, ik]
        Ek = @view E[:, ik]
        Hᵏ[ik] .= Hermitian(Uk' * Diagonal(Ek) * Uk)
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

@doc raw"""
    diag_Hk(Hᵏ:: AbstractArray{T, 3}) where {T<:Complex}

Diagonalize k space Hamiltonian `H`.

```math
H = V E V^{-1}
```

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
        ϵ, v = eigen(Hermitian(Hᵏ))
        E[:, ik] = ϵ
        V[:, :, ik] = v
    end

    return E, V
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

## NEW
abstract type AbstractKGrid{T} end

core_kgrid(x::AbstractKGrid) = x.core

k_cryst(x::AbstractKGrid) = core_kgrid(x).k_cryst
k_cryst(x::Vec3) = x

Base.length(kgrid::AbstractKGrid) = length(core_kgrid(kgrid))

struct CoreKGrid{T} <: AbstractKGrid{T}
    k_cryst::Vector{Vec3{T}}
end

core_kgrid(x::CoreKGrid) = x

Base.length(x::CoreKGrid) = length(x.k_cryst)

struct HamiltonianKGrid{T,MT<:AbstractMatrix{Complex{T}},VT<:AbstractVector{T}} <:
       AbstractKGrid{T}
       
    core::CoreKGrid{T}
    
    Hk::Vector{MT}
    eigvals::Vector{VT}
    eigvecs::Vector{MT}
end
function HamiltonianKGrid(kpoints::Vector{<:Vec3}, args...)
    return HamiltonianKGrid(CoreKGrid(kpoints), args...)
end

@doc raw"""
	HamiltonianKGrid(hami::TBHamiltonian{T}, nk, H_function_k::Function = x -> nothing) where T
	HamiltonianKGrid(hami::TBHamiltonian{T}, k_grid, H_function_k::Function = x -> nothing) where T

Takes a k grid, calculates Hk for each of them and diagonalizes. Only the eigenvectors and eigenvalues of Hk are stored,
the `H_function_k` function is called on the intermediate Hk. 
"""
function HamiltonianKGrid(hami::TBHamiltonian{T}, kpoints::Vector{<:Vec3},
                          Hk_function::Function = x -> nothing) where {T}
                          
    n_eigvals = max(blocksize(hami)...)
    eigvals = [zeros(T, n_eigvals) for k in kpoints]
    # eigvals = hami[1].block isa AbstractMagneticMatrix ?
    #           [MagneticVector(zeros(T, n_eigvals)) for k in kpoints] :
    #           [zeros(T, n_eigvals) for k in kpoints]
              
    kgrid = HamiltonianKGrid(kpoints, [zeros_block(hami) for k in kpoints], eigvals,
                             [zeros_block(hami) for k in kpoints])
    nk = length(kpoints)
    calc_caches = [HermitianEigenWs(block(hami[1])) for i in 1:Threads.nthreads()]
    p = Progress(nk, 1, "Calculating H(k)...")
    
    Threads.@threads for i in 1:nk
        tid = Threads.threadid()
        
        Hk!(kgrid.eigvecs[i], hami, k_cryst(kgrid)[i])
        
        copy!(kgrid.Hk[i], copy(kgrid.eigvecs[i]))
        Hk_function(kgrid.Hk[i])
        eigen!(kgrid.eigvals[i], kgrid.eigvecs[i], calc_caches[tid])
        next!(p)
    end
    
    return kgrid
end

n_wann(h::HamiltonianKGrid) = size(h.Hk[1],1)

function Hk!(out::AbstractMatrix, tbhami::TBHamiltonian, kpoint::Vec3)
    fill!(out, zero(eltype(out)))
    
    invfourier(tbhami, kpoint) do i, iR, R_cart, b, fac
        @inbounds out[i] += fac * b.block[i]
    end
    
end

"""
    Hk(hamiltonian::TBHamiltonian, kpoint::Vec3)
    Hk!(hk::AbstractMatrix, hamiltonian::TBHamiltonian, kpoint::Vec3)

Constructs the reciprocal Hamiltonian at a given _k_-point.  
"""
function Hk(tbhami::TBHamiltonian, kpoint::Vec3)
    out = similar(tbhami[1].block)
    Hk!(out, tbhami, kpoint)
    return out
end

Hk(g::HamiltonianKGrid) = g.Hk
eigvecs(g::HamiltonianKGrid) = g.eigvecs
eigvals(g::HamiltonianKGrid) = g.eigvals

"""
    interpolate(model::InterpModel{T}, kpoints::Matrix{T}) where {T<:Real}

Interpolate energy eigenvalues at `kpoints`.

# Arguments
- `model`: `InterpModel`
- `kpoints`: `3 * n_kpts`, kpoints to be interpolated, in fractional coordinates,
    can be nonuniform.
"""
function interpolate(hami::TBHamiltonian{T}, kpoints) where {T<:Real}
    grid = HamiltonianKGrid(hami, kpoints)
    # diagonalize
    Eᵏ = [grid.eigvals[ik][iw] for iw in 1:n_wann(grid), ik in 1:length(grid)]
    return Eᵏ
end

@inline function LinearAlgebra.eigen!(vecs::AbstractMatrix, ws::HermitianEigenWs)
    return Eigen(decompose!(ws, 'V', 'A', 'U', vecs, 0., 0., 0, 0, 1e-16)...)
end
@inline function LinearAlgebra.eigen!(vals, vecs::AbstractMatrix, ws::HermitianEigenWs)
    ws.w = vals
    te = Eigen(decompose!(ws, 'V', 'A', 'U', vecs, 0., 0., 0, 0, 1e-16)...)
    return te
end


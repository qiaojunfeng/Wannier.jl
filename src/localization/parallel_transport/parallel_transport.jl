export parallel_transport

include("contraction.jl")

"""
Löwdin-orthonormalized obstruction matrix between kpoints `k1` and `k2` that are
separated by `dk` (fractional coordinates).

`b` is recomputed from the kpoint coordinates rather than assumed to be a unit
vector, so that grids with negative coordinates (e.g. `-0.25` instead of `0.75`,
as produced by `qe/open_grid.x`) are handled correctly.
"""
function overlap_obstruction(U, M, bvectors, k1::Integer, k2::Integer, dk)
    b = round.(Int, bvectors.kpoints[k1] + dk - bvectors.kpoints[k2])
    ib = index_bvector(bvectors, k1, k2, b)
    Nᵏᵇ = view(U, :, :, k1)' * view(M, :, :, ib, k1) * view(U, :, :, k2)
    return orthonorm_lowdin(Nᵏᵇ)
end

"""
Right-multiply the gauge at kpoint `ik` by the pulled-back corner obstruction
`exp(t · log O) = V · diag(exp(t · logd)) · V'`, where `(V, logd)` come from
[`eig_log`](@ref).
"""
function pullback!(U, ik::Integer, V, logd, t::Real)
    Oₖ = V * Diagonal(exp.(t .* logd)) * V'
    view(U, :, :, ik) .= view(U, :, :, ik) * Oₖ
    return nothing
end

"""
Factor the U(1) determinant winding out of an edge obstruction path
`O_path[:, :, i]`, so that the SU(N) contraction ([`matrix_transport`](@ref) or
[`powm`](@ref)) only deals with the special-unitary part. Return
`(O_normalized, logD)`, where `logD` is the continuous (unwrapped) phase of
`det(O_path)` along the path.

The winding is spread *equally over all `n_col` bands* (an `n_col`-th root),
rather than dumped onto a single band, and is reattached afterwards with
[`det_winding_phase`](@ref). This is applied identically to every edge (xy, xz,
yz), so no edge is treated specially.
"""
function factor_det_winding(O_path::AbstractArray{Complex{T}, 3}) where {T <: Real}
    n_col = size(O_path, 1)
    n_k = size(O_path, 3)
    logD = [imag(log(det(view(O_path, :, :, i)))) for i in 1:n_k]
    for i in 2:n_k
        kmin = argmin([abs(logD[i] + 2π * k - logD[i - 1]) for k in -1:1])
        logD[i] += (kmin - 2) * 2π
    end
    O_norm = similar(O_path)
    for i in 1:n_k
        O_norm[:, :, i] = exp(-im * logD[i] / n_col) * view(O_path, :, :, i)
    end
    return O_norm, logD
end

"""
Scalar phase reattaching the U(1) determinant winding removed by
[`factor_det_winding`](@ref), spread equally over the `n_col` bands, for
unwrapped determinant phase `logD_i` and interpolation parameter `t`.
"""
det_winding_phase(logD_i::Real, t::Real, n_col::Integer) = exp(im * logD_i * t / n_col)

"""
    parallel_transport(model::Model{T}; use_U=false, log_interp=false)

Parallel transport the gauge from the first kpoint to all other kpoints.

Assumptions:
- the kpoints are contained in a `N1 * N2 * N3` cartesian grid
- the neighbor list must contain the six cartesian neighbors along x, y, z directions

# Arguments
- `model`: model

# Keyword arguments
- `use_U`: use the gauge `U` instead of random matrix
- `log_interp`: use logarithmic interpolation method
"""
function parallel_transport(
        model::Model{T}; use_U::Bool = false, log_interp::Bool = false
    ) where {T <: Real}
    if log_interp
        println("log interpolation")
    else
        println("parallel transport")
    end

    n_kx, n_ky, n_kz = kgrid_size(model)
    n_kpts = n_kpoints(model)
    n_wann = n_wannier(model)
    kpts = kpoints(model)

    # start from 0
    tx = collect(0:(n_kx - 1)) / n_kx
    ty = collect(0:(n_ky - 1)) / n_ky
    tz = collect(0:(n_kz - 1)) / n_kz

    k_xyz, xyz_k = get_kpoint_mappings(kpts, kgrid_size(model))

    # for overlap matrices
    M = model.overlaps
    bvectors = model.kstencil

    # the new gauge
    if use_U
        U = deepcopy(model.gauges)
    else
        U = identity_gauge(Complex{T}, n_kpts, n_wann)
    end

    # 1. propagate along kx
    @info "Filling (kx,0,0)"
    dk = [1 / n_kx, 0.0, 0.0]
    propagate!(U, [xyz_k[i, 1, 1] for i in 1:n_kx], dk, M, bvectors)

    # Corner obstruction for dimension d = 1.
    # In GLS2019 paper, ũ(1) = ( τ₁ ũ(0) ) Vₒ, where Vₒ is the obstruction matrix.
    # Here we compute Vₒ = U(nkx)' * M(nkx,1) * U(0), which is the inverse of the
    # Vₒ in the paper, so no minus sign is needed when pulling it back.
    O1 = overlap_obstruction(U, M, bvectors, xyz_k[end, 1, 1], xyz_k[1, 1, 1], dk)
    V, logd = eig_log(O1)
    for i in 1:n_kx
        pullback!(U, xyz_k[i, 1, 1], V, logd, tx[i])
    end

    # 2. propagate along ky
    @info "Filling (kx,ky,0)"
    dk = [0.0, 1 / n_ky, 0.0]
    for ik in 1:n_kx
        propagate!(U, [xyz_k[ik, j, 1] for j in 1:n_ky], dk, M, bvectors)
    end

    # corner obstruction
    O2 = overlap_obstruction(U, M, bvectors, xyz_k[1, end, 1], xyz_k[1, 1, 1], dk)
    V, logd = eig_log(O2)
    for i in 1:n_kx, j in 1:n_ky
        pullback!(U, xyz_k[i, j, 1], V, logd, ty[j])
    end

    # Line obstruction, at ky = 1 along kx = 0 -> 1
    Oxy = zeros(Complex{T}, n_wann, n_wann, n_kx)
    for i in 1:n_kx
        Oxy[:, :, i] = overlap_obstruction(U, M, bvectors, xyz_k[i, n_ky, 1], xyz_k[i, 1, 1], dk)
    end
    Oxy, logDxy = factor_det_winding(Oxy)

    Uxy = log_interp ? zeros(Complex{T}, n_wann, n_wann, n_kx, n_ky) : matrix_transport(Oxy, ty)
    for i in 1:n_kx, j in 1:n_ky
        base = log_interp ? powm(Oxy[:, :, i], ty[j]) : Uxy[:, :, i, j]
        ik = xyz_k[i, j, 1]
        view(U, :, :, ik) .= view(U, :, :, ik) * (det_winding_phase(logDxy[i], ty[j], n_wann) * base)
    end

    # 3. Propagate along the third dimension
    @info "Filling (k1,k2,k3)"
    dk = [0.0, 0.0, 1 / n_kz]
    for i in 1:n_kx, j in 1:n_ky
        propagate!(U, [xyz_k[i, j, k] for k in 1:n_kz], dk, M, bvectors)
    end

    # Fix corner
    O4 = overlap_obstruction(U, M, bvectors, xyz_k[1, 1, n_kz], xyz_k[1, 1, 1], dk)
    V, logd = eig_log(O4)
    for k in 1:n_kz, i in 1:n_kx, j in 1:n_ky
        pullback!(U, xyz_k[i, j, k], V, logd, tz[k])
    end

    # Fix first edge, in x-z plane, at kz = 1 along kx = 0 -> 1
    Oxz = zeros(Complex{T}, n_wann, n_wann, n_kx)
    for i in 1:n_kx
        Oxz[:, :, i] = overlap_obstruction(U, M, bvectors, xyz_k[i, 1, n_kz], xyz_k[i, 1, 1], dk)
    end
    Oxz, logDxz = factor_det_winding(Oxz)
    Uxz = log_interp ? zeros(Complex{T}, n_wann, n_wann, n_kx, n_kz) : matrix_transport(Oxz, tz)
    for i in 1:n_kx, k in 1:n_kz
        base = log_interp ? powm(Oxz[:, :, i], tz[k]) : Uxz[:, :, i, k]
        W = det_winding_phase(logDxz[i], tz[k], n_wann) * base
        for j in 1:n_ky
            ik = xyz_k[i, j, k]
            view(U, :, :, ik) .= view(U, :, :, ik) * W
        end
    end

    # Fix second edge, in y-z plane, at kz = 1 along ky = 0 -> 1
    Oyz = zeros(Complex{T}, n_wann, n_wann, n_ky)
    for j in 1:n_ky
        Oyz[:, :, j] = overlap_obstruction(U, M, bvectors, xyz_k[1, j, n_kz], xyz_k[1, j, 1], dk)
    end
    Oyz, logDyz = factor_det_winding(Oyz)
    Uyz = log_interp ? zeros(Complex{T}, n_wann, n_wann, n_ky, n_kz) : matrix_transport(Oyz, tz)
    for j in 1:n_ky, k in 1:n_kz
        base = log_interp ? powm(Oyz[:, :, j], tz[k]) : Uyz[:, :, j, k]
        W = det_winding_phase(logDyz[j], tz[k], n_wann) * base
        for i in 1:n_kx
            ik = xyz_k[i, j, k]
            view(U, :, :, ik) .= view(U, :, :, ik) * W
        end
    end

    # Fix whole surface
    for i in 1:n_kx, j in 1:n_ky
        O = overlap_obstruction(U, M, bvectors, xyz_k[i, j, n_kz], xyz_k[i, j, 1], dk)
        for k in 1:n_kz
            ik = xyz_k[i, j, k]
            view(U, :, :, ik) .= view(U, :, :, ik) * powm(O, tz[k])
        end
    end

    compute_error(model, U)

    return U, nothing
end

"""
    compute_error(model, U::Array{Complex{T},3})

Compute the smoothness error of the gauge.
"""
function compute_error(model::Model{T}, U::AbstractArray{Complex{T}, 3}) where {T <: Real}
    # initial (ϵ0) and final (ϵ1) error
    ϵ0 = 0.0
    ϵ1 = 0.0

    n_kx, n_ky, n_kz = kgrid_size(model)
    kpts = kpoints(model)
    k_xyz, xyz_k = get_kpoint_mappings(kpts, kgrid_size(model))

    M = model.overlaps
    U0 = model.gauges
    bvectors = model.kstencil

    epsilon(k1, k2, dk, B) = norm(overlap_obstruction(B, M, bvectors, k1, k2, dk) - I)^2

    dkx = [1 / n_kx, 0.0, 0.0]
    dky = [0.0, 1 / n_ky, 0.0]
    dkz = [0.0, 0.0, 1 / n_kz]

    for i in 1:n_kx, j in 1:n_ky, k in 1:n_kz
        k1 = xyz_k[i, j, k]

        k2 = xyz_k[i == n_kx ? 1 : i + 1, j, k]
        ϵ0 += epsilon(k1, k2, dkx, U0)
        ϵ1 += epsilon(k1, k2, dkx, U)

        k2 = xyz_k[i, j == n_ky ? 1 : j + 1, k]
        ϵ0 += epsilon(k1, k2, dky, U0)
        ϵ1 += epsilon(k1, k2, dky, U)

        k2 = xyz_k[i, j, k == n_kz ? 1 : k + 1]
        ϵ0 += epsilon(k1, k2, dkz, U0)
        ϵ1 += epsilon(k1, k2, dkz, U)
    end

    ϵ0 = sqrt(ϵ0) / n_kpoints(model)
    ϵ1 = sqrt(ϵ1) / n_kpoints(model)

    println("initial error = ", round(ϵ0; digits = 4))
    println("final error   = ", round(ϵ1; digits = 4))
    println()

    return ϵ0, ϵ1
end

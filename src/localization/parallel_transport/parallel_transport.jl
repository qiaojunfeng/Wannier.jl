export parallel_transport

include("contraction.jl")

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
    model::Model{T}; use_U::Bool=false, log_interp::Bool=false
) where {T<:Real}
    if log_interp
        println("log interpolation")
    else
        println("parallel transport")
    end

    n_kx, n_ky, n_kz = model.kgrid
    n_kpts = model.n_kpts
    n_wann = model.n_wann
    kpoints = model.kpoints

    # start from 0
    tx = collect(0:(n_kx - 1)) / n_kx
    ty = collect(0:(n_ky - 1)) / n_ky
    tz = collect(0:(n_kz - 1)) / n_kz

    k_xyz, xyz_k = get_kpoint_mappings(kpoints, model.kgrid)

    # for overlap matrices
    M = model.M
    bvectors = model.kstencil

    # the new gauge
    if use_U
        U = deepcopy(model.U)
    else
        U = identity_gauge(Complex{T}, n_kpts, n_wann)
    end

    # 1. propagate along kx
    @info "Filling (kx,0,0)"
    kpts = [xyz_k[i, 1, 1] for i in 1:n_kx]
    dk = [1 / n_kx, 0.0, 0.0]
    propagate!(U, kpts, dk, M, bvectors)

    # compute obstruction matrix for dimension d = 1
    # In GLS2019 paper, ũ(1) = ( τ₁ ũ(0) ) Vₒ, where Vₒ is the obstruction matrix.
    # However, here we compute Vₒ = U(nkx)' * M(nkx,1) * U(0),
    # since kpoints are discretized 1...nkx, we need approximations:
    #     ũ(1) ≈ |ψ(nkx)> * U(nkx),
    #     τ₁ ũ(0) = |ψ(1+nkx)> * U(1),
    # so our Vₒ is actually the inverse of the Vₒ in the paper.
    k1 = xyz_k[end, 1, 1]
    k2 = xyz_k[1, 1, 1]

    # Compute shifting b vector, for regular MP grid, b = [1, 0, 0].
    # However, for grids having negative cooridnates (i.e, -0.25 instead of 0.75),
    # we need to recompute b vector.
    # For MP grid, k1 + dk == k2 + [1, 0, 0], so b = [1, 0, 0].
    # For other grids, b = k1 + dk - k2.
    b = round.(Int, kpoints[k1] + dk - kpoints[k2])
    ib = index_bvector(bvectors, k1, k2, b)
    Nᵏᵇ = U[k1]' * M[k1][ib] * U[k2]
    O1 = orthonorm_lowdin(Nᵏᵇ)
    @debug "Obstruction matrix =" V = O1

    d, V = eigen(O1)
    logd = log.(d)
    # Hack to avoid separating eigenvalues at -1. TODO understand that
    for i in 1:n_wann
        if imag(logd[i]) < -π + 0.01
            logd[i] += 2π * im
        end
    end
    @debug "log(d) =" logd

    # and pull it back
    for i in 1:n_kx
        ik = xyz_k[i, 1, 1]
        # Since our Vₒ is the inverse of the Vₒ in the paper,
        # we don't need a minus sign here.
        Oₖ = V * diagm(0 => exp.(tx[i] * logd)) * V'
        U[ik] *= Oₖ
    end

    # 2. propagate along ky
    @info "Filling (kx,ky,0)"
    dk = [0.0, 1 / n_ky, 0.0]
    for ik in 1:n_kx
        kpts = [xyz_k[ik, j, 1] for j in 1:n_ky]
        propagate!(U, kpts, dk, M, bvectors)
    end

    # corner obstruction
    k1 = xyz_k[1, end, 1]
    k2 = xyz_k[1, 1, 1]
    # For MP grid, b = [0, 1, 0]
    # For other grids,
    b = round.(Int, kpoints[k1] + dk - kpoints[k2])
    ib = index_bvector(bvectors, k1, k2, b)
    Nᵏᵇ = U[k1]' * M[k1][ib] * U[k2]
    O2 = orthonorm_lowdin(Nᵏᵇ)

    d, V = eigen(O2)
    logd = log.(d)
    # Hack to avoid separating eigenvalues at -1. TODO understand that
    for i in 1:n_wann
        if imag(logd[i]) < -π + 0.01
            logd[i] += 2π * im
        end
    end

    # pull it back
    for i in 1:n_kx
        for j in 1:n_ky
            ik = xyz_k[i, j, 1]
            # no need a minus sign here
            Oₖ = V * diagm(0 => exp.(ty[j] * logd)) * V'
            U[ik] *= Oₖ
        end
    end

    # pull back the line obstruction, at ky = 1 along kx = 0 -> 1
    Oxy = zeros(Complex{T}, n_wann, n_wann, n_kx)
    detO3 = zeros(Complex{T}, n_kx)
    # eigs = zeros(Complex{T}, n_wann, n_kx)
    for i in 1:n_kx
        k1 = xyz_k[i, n_ky, 1]
        k2 = xyz_k[i, 1, 1]
        # For MP grid, b = [0, 1, 0]
        # For other grids,
        b = round.(Int, kpoints[k1] + dk - kpoints[k2])
        ib = index_bvector(bvectors, k1, k2, b)
        Nᵏᵇ = U[k1]' * M[k1][ib] * U[k2]
        Oxy[:, :, i] = orthonorm_lowdin(Nᵏᵇ)
        detO3[i] = det(Oxy[:, :, i])
    end

    # find a continuous log of the determinant
    logD = imag(log.(detO3))
    for i in 2:n_kx
        kmin = argmin([abs(logD[i] + 2π * k - logD[i - 1]) for k in -1:1])
        logD[i] += (kmin - 2) * 2π
    end
    for i in 1:n_kx
        Oxy[:, :, i] = exp(-im * logD[i] / n_wann) * Oxy[:, :, i]
    end
    # Interpolate the line obstruction
    Uxy = zeros(Complex{T}, n_wann, n_wann, n_kx, n_ky)

    if !log_interp
        Uxy = matrix_transport(Oxy, ty)
    end
    for i in 1:n_kx
        O = Oxy[:, :, i]
        d, V = eigen(O)

        if log_interp
            for j in 1:n_ky
                Uxy[:, :, i, j] = powm(O, ty[j])
            end
        end

        for j in 1:n_ky
            ik = xyz_k[i, j, 1]

            U[ik] *= exp(im * logD[i] * ty[j] / 2)
            U[ik] *= Uxy[:, :, i, j]
        end
    end

    # Propagate along the third dimension
    @info "Filling (k1,k2,k3)"
    dk = [0.0, 0.0, 1 / n_kz]
    for i in 1:n_kx, j in 1:n_ky
        kpts = [xyz_k[i, j, k] for k in 1:n_kz]
        propagate!(U, kpts, dk, M, bvectors)
    end

    # Fix corner
    k1 = xyz_k[1, 1, n_kz]
    k2 = xyz_k[1, 1, 1]
    # For MP grid, b = [0, 0, 1]
    # For other grids,
    b = round.(Int, kpoints[k1] + dk - kpoints[k2])
    ib = index_bvector(bvectors, k1, k2, b)
    Nᵏᵇ = U[k1]' * M[k1][ib] * U[k2]
    O4 = orthonorm_lowdin(Nᵏᵇ)
    d, V = eigen(O4)
    logd = log.(d)

    # Hack to avoid separating eigenvalues at -1. TODO understand that
    for i in 1:n_wann
        if imag(logd[i]) < -π + 0.01
            logd[i] += 2π * im
        end
    end

    for k in 1:n_kz
        # fixer = powm(O4, t3[k])
        W = V * diagm(0 => exp.(tz[k] * logd)) * V'

        for i in 1:n_kx, j in 1:n_ky
            ik = xyz_k[i, j, k]
            U[ik] *= W
        end
    end

    # Fix first edge, in x-z plane, at kz = 1 along kx = 0 -> 1
    Oxz = zeros(Complex{T}, n_wann, n_wann, n_kx)
    Uxz = zeros(Complex{T}, n_wann, n_wann, n_kx, n_kz)

    for i in 1:n_kx
        k1 = xyz_k[i, 1, n_kz]
        k2 = xyz_k[i, 1, 1]
        # For MP grid, b = [0, 0, 1]
        # For other grids,
        b = round.(Int, kpoints[k1] + dk - kpoints[k2])
        ib = index_bvector(bvectors, k1, k2, b)
        Nᵏᵇ = U[k1]' * M[k1][ib] * U[k2]
        Oxz[:, :, i] = orthonorm_lowdin(Nᵏᵇ)
    end

    if !log_interp
        Uxz = matrix_transport(Oxz, tz)
    end

    for i in 1:n_kx
        for k in 1:n_kz
            if log_interp
                W = powm(Oxz[:, :, i], tz[k])

                for j in 1:n_ky
                    ik = xyz_k[i, j, k]
                    U[ik] *= W
                end
            else
                for j in 1:n_ky
                    ik = xyz_k[i, j, k]
                    U[ik] *= Uxz[:, :, i, k]
                end
            end
        end
    end

    # Fix second edge, in y-z plane, at kz = 1 along ky = 0 -> 1
    Oyz = zeros(Complex{T}, n_wann, n_wann, n_ky)
    Uyz = zeros(Complex{T}, n_wann, n_wann, n_ky, n_kz)

    for j in 1:n_ky
        k1 = xyz_k[1, j, n_kz]
        k2 = xyz_k[1, j, 1]
        # For MP grid, b = [0, 0, 1]
        # For other grids,
        b = round.(Int, kpoints[k1] + dk - kpoints[k2])
        ib = index_bvector(bvectors, k1, k2, b)
        Nᵏᵇ = U[k1]' * M[k1][ib] * U[k2]
        Oyz[:, :, j] = orthonorm_lowdin(Nᵏᵇ)
    end

    if !log_interp
        Uyz = matrix_transport(Oyz, tz)
    end

    for j in 1:n_ky
        for k in 1:n_kz
            if log_interp
                W = powm(Oyz[:, :, j], tz[k])

                for i in 1:n_kx
                    ik = xyz_k[i, j, k]
                    U[ik] *= W
                end
            else
                for i in 1:n_kx
                    ik = xyz_k[i, j, k]
                    U[ik] *= Uyz[:, :, j, k]
                end
            end
        end
    end

    # Fix whole surface
    for i in 1:n_kx, j in 1:n_ky
        k1 = xyz_k[i, j, n_kz]
        k2 = xyz_k[i, j, 1]
        # For MP grid, b = [0, 0, 1]
        # For other grids,
        b = round.(Int, kpoints[k1] + dk - kpoints[k2])
        ib = index_bvector(bvectors, k1, k2, b)
        Nᵏᵇ = U[k1]' * M[k1][ib] * U[k2]
        O = orthonorm_lowdin(Nᵏᵇ)

        for k in 1:n_kz
            ik = xyz_k[i, j, k]
            U[ik] *= powm(O, tz[k])
        end
    end

    compute_error(model, U)

    obs = Obstruction(Oxy, Oxz, Oyz, Uxy, Uxz, Uyz)
    return U, obs
end

"""
    compute_error(model, U::Array{Complex{T},3})

Compute the smoothness error of the gauge.
"""
function compute_error(model::Model{T}, U::Vector{Matrix{Complex{T}}}) where {T<:Real}
    # initial error
    ϵ0 = 0.0
    # final error
    ϵ1 = 0.0

    n_kx, n_ky, n_kz = model.kgrid
    kpoints = model.kpoints
    k_xyz, xyz_k = get_kpoint_mappings(kpoints, model.kgrid)

    M = model.M
    U0 = model.U

    epsilon(i, j, b, B) = begin
        ib = index_bvector(model.kstencil, i, j, b)
        Nᵏᵇ = B[i]' * M[i][ib] * B[j]
        norm(orthonorm_lowdin(Nᵏᵇ) - I)^2
    end

    dkx = [1 / n_kx, 0.0, 0.0]
    dky = [0.0, 1 / n_ky, 0.0]
    dkz = [0.0, 0.0, 1 / n_kz]

    for i in 1:n_kx, j in 1:n_ky, k in 1:n_kz
        k1 = xyz_k[i, j, k]
        if i == n_kx
            k2 = xyz_k[1, j, k]
            # For MP grid
            # b = [1, 0, 0]
        else
            k2 = xyz_k[i + 1, j, k]
            # For MP grid
            # b = [0, 0, 0]
        end
        # For other grids
        b = round.(Int, kpoints[k1] + dkx - kpoints[k2])

        ϵ0 += epsilon(k1, k2, b, U0)
        ϵ1 += epsilon(k1, k2, b, U)

        if j == n_ky
            k2 = xyz_k[i, 1, k]
            # For MP grid
            # b = [0, 1, 0]
        else
            k2 = xyz_k[i, j + 1, k]
            # For MP grid
            # b = [0, 0, 0]
        end
        # For other grids
        b = round.(Int, kpoints[k1] + dky - kpoints[k2])
        ϵ0 += epsilon(k1, k2, b, U0)
        ϵ1 += epsilon(k1, k2, b, U)

        if k == n_kz
            k2 = xyz_k[i, j, 1]
            # For MP grid
            # b = [0, 0, 1]
        else
            k2 = xyz_k[i, j, k + 1]
            # For MP grid
            # b = [0, 0, 0]
        end
        # For other grids
        b = round.(Int, kpoints[k1] + dkz - kpoints[k2])

        ϵ0 += epsilon(k1, k2, b, U0)
        ϵ1 += epsilon(k1, k2, b, U)
    end

    ϵ0 = sqrt(ϵ0) / model.n_kpts
    ϵ1 = sqrt(ϵ1) / model.n_kpts

    println("initial error = ", round(ϵ0; digits=4))
    println("final error   = ", round(ϵ1; digits=4))
    println()

    return ϵ0, ϵ1
end

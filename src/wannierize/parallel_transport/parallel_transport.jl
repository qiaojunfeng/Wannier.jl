import LinearAlgebra as LA

include("contraction.jl")


"""
Assumptions: the kpoints are contained in a N1 x N2 x N3 cartesian grid,
the neighbor list must contain the six cartesian neighbors along x,y,z directions.
"""
function parallel_transport(
    model::Model{T};
    use_A::Bool = false,
    log_interp::Bool = false,
    return_obs::Bool = false,
) where {T<:Real}

    if log_interp
        println("log interpolation")
    else
        println("parallel transport")
    end

    n_kx, n_ky, n_kz = model.kgrid
    n_kpts = model.n_kpts
    n_wann = model.n_wann

    # start from 0
    tx = collect(0:n_kx-1) / n_kx
    ty = collect(0:n_ky-1) / n_ky
    tz = collect(0:n_kz-1) / n_kz

    k_xyz, xyz_k = get_kpoint_mappings(model.kpoints, model.kgrid)

    # for overlap matrices
    M = model.M
    kpb_k = model.bvectors.kpb_k

    # the new gauge
    if use_A
        A = deepcopy(model.A)
    else
        A = zeros(Complex{T}, n_wann, n_wann, n_kpts)
        for ik = 1:n_kpts
            for i = 1:n_wann
                A[i, i, ik] = 1
            end
        end
    end

    # 1. propagate along kx
    @info "Filling (kx,0,0)"
    kpts = [xyz_k[i, 1, 1] for i = 1:n_kx]
    propagate!(A, kpts, M, kpb_k)

    # compute obstruction matrix for dimension d = 1
    # In GLS2019 paper, ũ(1) = ( τ₁ ũ(0) ) Vₒ, where Vₒ is the obstruction matrix.
    # However, here we compute Vₒ = A(nkx)' * M(nkx,1) * A(0),
    # since kpoints are discretized 1...nkx, we need approximations:
    #     ũ(1) ≈ |ψ(nkx)> * A(nkx),
    #     τ₁ ũ(0) = |ψ(1+nkx)> * A(1), 
    # so our Vₒ is actually the inverse of the Vₒ in the paper.
    k1 = xyz_k[end, 1, 1]
    k2 = xyz_k[1, 1, 1]
    O1 = orthonorm_lowdin(overlap(M, kpb_k, k1, k2, A))
    @debug "Obstruction matrix =" V = O1

    d, V = LA.eigen(O1)
    logd = log.(d)
    # Hack to avoid separating eigenvalues at -1. TODO understand that
    for i = 1:n_wann
        if imag(logd[i]) < -π + 0.01
            logd[i] += 2π * im
        end
    end
    @debug "log(d) =" logd

    # and pull it back
    for i = 1:n_kx
        ik = xyz_k[i, 1, 1]
        # Since our Vₒ is the inverse of the Vₒ in the paper,
        # we don't need a minus sign here.
        Oₖ = V * LA.diagm(0 => exp.(tx[i] * logd)) * V'
        A[:, :, ik] *= Oₖ
    end

    # 2. propagate along ky
    @info "Filling (kx,ky,0)"
    for ik = 1:n_kx
        kpts = [xyz_k[ik, j, 1] for j = 1:n_ky]
        propagate!(A, kpts, M, kpb_k)
    end

    # corner obstruction
    k1 = xyz_k[1, end, 1]
    k2 = xyz_k[1, 1, 1]
    O2 = orthonorm_lowdin(overlap(M, kpb_k, k1, k2, A))

    d, V = LA.eigen(O2)
    logd = log.(d)
    # Hack to avoid separating eigenvalues at -1. TODO understand that
    for i = 1:n_wann
        if imag(logd[i]) < -π + 0.01
            logd[i] += 2π * im
        end
    end

    # pull it back
    for i = 1:n_kx
        for j = 1:n_ky
            ik = xyz_k[i, j, 1]
            # no need a minus sign here
            Oₖ = V * LA.diagm(0 => exp.(ty[j] * logd)) * V'
            A[:, :, ik] *= Oₖ
        end
    end

    # pull back the line obstruction, at ky = 1 along kx = 0 -> 1
    Oxy = zeros(Complex{T}, n_wann, n_wann, n_kx)
    detO3 = zeros(Complex{T}, n_kx)
    # eigs = zeros(Complex{T}, n_wann, n_kx)
    for i = 1:n_kx
        k1 = xyz_k[i, n_ky, 1]
        k2 = xyz_k[i, 1, 1]
        Mᵏᵇ = overlap(M, kpb_k, k1, k2, A)
        Oxy[:, :, i] = orthonorm_lowdin(Mᵏᵇ)
        detO3[i] = LA.det(Oxy[:, :, i])
    end

    # find a continuous log of the determinant
    logD = imag(log.(detO3))
    for i = 2:n_kx
        kmin = argmin([abs(logD[i] + 2π * k - logD[i-1]) for k = -1:1])
        logD[i] += (kmin - 2) * 2π
    end
    for i = 1:n_kx
        Oxy[:, :, i] = exp(-im * logD[i] / n_wann) * Oxy[:, :, i]
        # eigs[:, i] = LA.eigvals(Oxy[:, :, i])
    end

    # Interpolate the line obstruction
    Uxy = zeros(Complex{T}, n_wann, n_wann, n_kx, n_ky)

    if !log_interp
        Uxy = matrix_transport(Oxy, ty)
    end

    for i = 1:n_kx
        O = Oxy[:, :, i]
        d, V = LA.eigen(O)

        if log_interp
            for j = 1:n_ky
                Uxy[:, :, i, j] = powm(O, ty[j])
            end
        end

        for j = 1:n_ky
            ik = xyz_k[i, j, 1]

            A[:, :, ik] *= exp(im * logD[i] * ty[j] / 2)
            A[:, :, ik] *= Uxy[:, :, i, j]
        end
    end

    # Propagate along the third dimension
    @info "Filling (k1,k2,k3)"
    for i = 1:n_kx, j = 1:n_ky
        kpts = [xyz_k[i, j, k] for k = 1:n_kz]
        propagate!(A, kpts, M, kpb_k)
    end

    # Fix corner
    k1 = xyz_k[1, 1, n_kz]
    k2 = xyz_k[1, 1, 1]
    O4 = orthonorm_lowdin(overlap(M, kpb_k, k1, k2, A))
    d, V = LA.eigen(O4)
    logd = log.(d)

    # Hack to avoid separating eigenvalues at -1. TODO understand that
    for i = 1:n_wann
        if imag(logd[i]) < -π + 0.01
            logd[i] += 2π * im
        end
    end

    for k = 1:n_kz
        # fixer = powm(O4, t3[k])
        W = V * LA.diagm(0 => exp.(tz[k] * logd)) * V'

        for i = 1:n_kx, j = 1:n_ky
            ik = xyz_k[i, j, k]
            A[:, :, ik] *= W
        end
    end

    # Fix first edge, in x-z plane, at kz = 1 along kx = 0 -> 1
    Oxz = zeros(Complex{T}, n_wann, n_wann, n_kx)
    Uxz = zeros(Complex{T}, n_wann, n_wann, n_kx, n_kz)

    for i = 1:n_kx
        k1 = xyz_k[i, 1, n_kz]
        k2 = xyz_k[i, 1, 1]
        Oxz[:, :, i] = orthonorm_lowdin(overlap(M, kpb_k, k1, k2, A))
    end

    if !log_interp
        Uxz = matrix_transport(Oxz, tz)
    end

    for i = 1:n_kx
        for k = 1:n_kz

            if log_interp
                W = powm(Oxz[:, :, i], tz[k])

                for j = 1:n_ky
                    ik = xyz_k[i, j, k]
                    A[:, :, ik] *= W
                end
            else
                for j = 1:n_ky
                    ik = xyz_k[i, j, k]
                    A[:, :, ik] *= Uxz[:, :, i, k]
                end
            end
        end
    end

    # Fix second edge, in y-z plane, at kz = 1 along ky = 0 -> 1
    Oyz = zeros(Complex{T}, n_wann, n_wann, n_ky)
    Uyz = zeros(Complex{T}, n_wann, n_wann, n_ky, n_kz)

    for j = 1:n_ky
        k1 = xyz_k[1, j, n_kz]
        k2 = xyz_k[1, j, 1]
        Oyz[:, :, j] = orthonorm_lowdin(overlap(M, kpb_k, k1, k2, A))
    end

    if !log_interp
        Uyz = matrix_transport(Oyz, tz)
    end

    for j = 1:n_ky
        for k = 1:n_kz

            if log_interp
                W = powm(Oyz[:, :, j], tz[k])

                for i = 1:n_kx
                    ik = xyz_k[i, j, k]
                    A[:, :, ik] *= W
                end
            else
                for i = 1:n_kx
                    ik = xyz_k[i, j, k]
                    A[:, :, ik] *= Uyz[:, :, j, k]
                end
            end
        end
    end

    # Fix whole surface
    for i = 1:n_kx, j = 1:n_ky
        k1 = xyz_k[i, j, n_kz]
        k2 = xyz_k[i, j, 1]
        O = orthonorm_lowdin(overlap(M, kpb_k, k1, k2, A))

        for k = 1:n_kz
            ik = xyz_k[i, j, k]
            A[:, :, ik] *= powm(O, tz[k])
        end
    end

    compute_error(model, A)

    if return_obs
        obs = Obstruction(Oxy, Oxz, Oyz, Uxy, Uxz, Uyz)
        return A, obs
    end

    A
end


function compute_error(model::Model{T}, A::Array{Complex{T},3}) where {T<:Real}
    # initial error
    ϵ0 = 0.0
    # final error
    ϵ1 = 0.0

    n_kx, n_ky, n_kz = model.kgrid
    kpb_k = model.bvectors.kpb_k

    k_xyz, xyz_k = get_kpoint_mappings(model.kpoints, model.kgrid)

    M = model.M
    A0 = model.A

    epsilon(i, j, B) = LA.norm(orthonorm_lowdin(overlap(M, kpb_k, i, j, B)) - LA.I)^2

    for i = 1:n_kx, j = 1:n_ky, k = 1:n_kz
        k1 = xyz_k[i, j, k]
        k2 = xyz_k[i%n_kx+1, j, k]
        ϵ0 += epsilon(k1, k2, A0)
        ϵ1 += epsilon(k1, k2, A)

        k2 = xyz_k[i, j%n_ky+1, k]
        ϵ0 += epsilon(k1, k2, A0)
        ϵ1 += epsilon(k1, k2, A)

        k2 = xyz_k[i, j, k%n_kz+1]
        ϵ0 += epsilon(k1, k2, A0)
        ϵ1 += epsilon(k1, k2, A)
    end

    ϵ0 = sqrt(ϵ0) / model.n_kpts
    ϵ1 = sqrt(ϵ1) / model.n_kpts


    println("initial error = ", round(ϵ0; digits=4))
    println("final error   = ", round(ϵ1; digits=4))
    println()

    ϵ0, ϵ1
end

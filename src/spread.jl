import LinearAlgebra as LA


@doc raw"""
compute the MV energy

From MV: Omega = sum_n <r2>n - |<r>n|^2
<r>n = -1/N sum_k,b wb b Im ln(Mnn,kb)
<r2>n = 1/N sum_k,b wb [(1 - |Mnn,kb|^2) + (Im ln(Mnn,kb))^2]
"""
struct Spread{T<:Real}
    # Total, Ω = ΩI + Ω̃
    Ω::T

    # gauge-invarient part
    ΩI::T

    # off-diagonal part
    ΩOD::T

    # diagonal part, not implemented yet
    ΩD::T

    # Ω̃ = ΩOD + ΩD
    Ω̃::T

    # Ω of each WF, n_wann
    ω::Vector{T}

    # WF center, cartesian coordinates, 3 * n_wann
    r::Matrix{T}

    # frozen_weight::T
    # fix_centers :: Array{Float64,2} #3 x nwannier
end


@views function omega(
    M::Array{Complex{FT},4},
    A::Array{Complex{FT},3},
    bvectors::BVectors{FT},
    only_r2::Bool = false,
) where {FT<:Real}

    n_bands, n_wann, n_kpts = size(A)
    n_bvecs = size(M, 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    # # keep in case we want to do this later on
    # μ::FT = 0.0
    # n_froz = 0
    # # frozen weight
    # w_froz::FT = 0.0

    r = zeros(FT, 3, n_wann)
    r² = zeros(FT, n_wann)

    ΩI::FT = 0.0
    ΩOD::FT = 0.0
    ΩD::FT = 0.0

    b = zeros(FT, 3)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    MAᵏᵇ = zeros(Complex{FT}, n_bands, n_wann)

    for ik = 1:n_kpts

        # w_froz -= μ * sum(abs2, A[1:n_froz, :, ik])

        for ib = 1:n_bvecs
            ikpb = kpb_k[ib, ik]

            MAᵏᵇ .= overlap(M, kpb_k, ik, ikpb) * A[:, :, ikpb]
            # compute-intensive, but should be true
            # @assert overlap(M, kpb_k, ik, ikpb)' ≈ overlap(M, kpb_k, ikpb, ik)
            Nᵏᵇ .= A[:, :, ik]' * MAᵏᵇ
            b .= recip_lattice * (kpoints[:, ikpb] + kpb_b[:, ib, ik] - kpoints[:, ik])

            w_ib = wb[ib]

            ΩI += w_ib * (n_wann - sum(abs2, Nᵏᵇ))
            ΩOD += w_ib * sum(abs2, Nᵏᵇ .- LA.diagm(0 => LA.diag(Nᵏᵇ)))

            for n = 1:n_wann
                if !only_r2
                    r[:, n] -= w_ib * imaglog(Nᵏᵇ[n, n]) * b
                end

                r²[n] += w_ib * (1 - abs(Nᵏᵇ[n, n])^2 + imaglog(Nᵏᵇ[n, n])^2)
                # r²[n] += w_ib * 2*(1 - real(Nᵏᵇ[n,n]))
            end
        end
    end

    r /= n_kpts
    r² /= n_kpts
    ΩI /= n_kpts
    ΩOD /= n_kpts
    ΩD /= n_kpts
    # w_froz /= n_kpts

    # @debug "Spreads" r r²' ΩI ΩOD ΩD

    # Ω of each WF
    ω = r² - dropdims(sum(abs.(r) .^ 2, dims = 1), dims = 1)
    # total Ω
    Ω = sum(ω)
    # Ω += w_froz
    Ω̃ = Ω - ΩI

    Spread(Ω, ΩI, ΩOD, ΩD, Ω̃, ω, r)
    # Spread(Ω, ΩI, ΩOD, ΩD, Ω̃, ω, r, w_froz)
end


"""
dΩ/dU, n_bands * n_wann * n_kpts

r: WF centers, cartesian coordinates, 3 * n_wann
"""
@views function omega_grad(
    M::Array{Complex{FT},4},
    A::Array{Complex{FT},3},
    bvectors::BVectors{FT},
    r::Matrix{FT},
    only_r2::Bool = false,
) where {FT<:Real}

    n_bands, n_wann, n_kpts = size(A)
    n_bvecs = size(M, 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    # # keep in case we want to do this later on
    # μ::FT = 0.0
    # n_froz = 0
    # # frozen weight
    # w_froz::FT = 0.0

    G = zeros(Complex{FT}, n_bands, n_wann, n_kpts)

    R = zeros(Complex{FT}, n_bands, n_wann)
    T = zeros(Complex{FT}, n_bands, n_wann)

    b = zeros(FT, 3)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    MAᵏᵇ = zeros(Complex{FT}, n_bands, n_wann)

    for ik = 1:n_kpts

        # w_froz -= μ * sum(abs2, A[1:n_froz, :, ik])
        # G[1:n_froz, :, ik] = -2 * μ * A[1:n_froz, :, ik]

        for ib = 1:n_bvecs
            ikpb = kpb_k[ib, ik]

            MAᵏᵇ .= overlap(M, kpb_k, ik, ikpb) * A[:, :, ikpb]
            Nᵏᵇ .= A[:, :, ik]' * MAᵏᵇ
            b .= recip_lattice * (kpoints[:, ikpb] + kpb_b[:, ib, ik] - kpoints[:, ik])

            w_ib = wb[ib]

            # MV way
            # fA(B) = (B - B') / 2
            # fS(B) = (B + B') / (2 * im)
            # q = imaglog.(LA.diag(Nᵏᵇ)) + r' * b
            # for m = 1:n_wann, n = 1:n_wann
            #     R[m, n] = Nᵏᵇ[m, n] * conj(Nᵏᵇ[n, n])
            #     T[m, n] = Nᵏᵇ[m, n] / Nᵏᵇ[n, n] * q[n]
            # end
            # G[:, :, ik] += 4 * w_ib * (fA(R) .- fS(T))

            q = imaglog.(LA.diag(Nᵏᵇ))
            if !only_r2
                q += r' * b
            end

            for n = 1:n_wann
                # error if division by zero. Should not happen if the initial gauge is not too bad
                if abs(Nᵏᵇ[n, n]) < 1e-10
                    display(Nᵏᵇ)
                    println()
                    error("Nᵏᵇ too small! $ik -> $ikpb")
                end

                Tfac = -im * q[n] / Nᵏᵇ[n, n]

                for m = 1:n_bands
                    R[m, n] = -MAᵏᵇ[m, n] * conj(Nᵏᵇ[n, n])
                    # T[m, n] = -im * MAᵏᵇ[m, n] / (Nᵏᵇ[n, n]) * q[n]
                    T[m, n] = Tfac * MAᵏᵇ[m, n]
                end
            end

            G[:, :, ik] .+= 4 * w_ib .* (R .+ T)
        end
    end

    G /= n_kpts

    G
end


@views function omega_grad(
    M::Array{Complex{FT},4},
    A::Array{Complex{FT},3},
    bvectors::BVectors{FT},
    only_r2::Bool = false,
) where {FT<:Real}
    r = omega(M, A, bvectors, only_r2).r
    omega_grad(M, A, bvectors, r, only_r2)
end


"""
local part of the contribution to r^2
"""
function omega_loc(
    M::Array{Complex{FT},4},
    A::Array{Complex{FT},3},
    bvectors::BVectors{FT},
) where {FT<:Real}
    n_bands, n_wann, n_kpts = size(A)
    n_bvecs = size(M, 3)

    kpb_k = bvectors.kpb_k
    wb = bvectors.weights

    loc = zeros(FT, n_kpts)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    MAᵏᵇ = zeros(Complex{FT}, n_bands, n_wann)

    for ik = 1:n_kpts
        for ib = 1:n_bvecs
            ikpb = kpb_k[ib, ik]

            MAᵏᵇ .= overlap(M, kpb_k, ik, ikpb) * A[:, :, ikpb]
            Nᵏᵇ .= A[:, :, ik]' * MAᵏᵇ

            for n = 1:n_wann
                loc[ik] += wb[ib] * (1 - abs(Nᵏᵇ[n, n])^2 + imaglog(Nᵏᵇ[n, n])^2)
            end
        end
    end

    loc
end


"""
WF postion operator matrix
"""
@views function position(
    M::Array{Complex{FT},4},
    A::Array{Complex{FT},3},
    bvectors::BVectors{FT},
) where {FT<:Real}

    n_bands, n_wann, n_kpts = size(A)
    n_bvecs = size(M, 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    # along x, y, z directions
    R = zeros(Complex{FT}, n_wann, n_wann, 3)

    b = zeros(FT, 3)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    MAᵏᵇ = zeros(Complex{FT}, n_bands, n_wann)

    for ik = 1:n_kpts
        for ib = 1:n_bvecs

            ikpb = kpb_k[ib, ik]

            MAᵏᵇ .= overlap(M, kpb_k, ik, ikpb) * A[:, :, ikpb]
            Nᵏᵇ .= A[:, :, ik]' * MAᵏᵇ
            b .= recip_lattice * (kpoints[:, ikpb] + kpb_b[:, ib, ik] - kpoints[:, ik])

            w_ib = wb[ib]

            for m = 1:n_wann
                for n = 1:n_wann
                    R[m, n, :] += w_ib * Nᵏᵇ[m, n] * b

                    if m == n
                        R[m, n, :] -= w_ib * b
                    end
                end
            end

        end
    end

    R /= -im * n_kpts

    R
end


"""
Berry connection at each kpoint
"""
@views function berry_connection(
    M::Array{Complex{FT},4},
    A::Array{Complex{FT},3},
    bvectors::BVectors{FT},
) where {FT<:Real}

    n_bands, n_wann, n_kpts = size(A)
    n_bvecs = size(M, 3)

    kpb_k = bvectors.kpb_k
    kpb_b = bvectors.kpb_b
    wb = bvectors.weights
    recip_lattice = bvectors.recip_lattice
    kpoints = bvectors.kpoints

    # along x, y, z directions
    A = zeros(Complex{FT}, n_wann, n_wann, 3, n_kpts)

    b = zeros(FT, 3)

    Nᵏᵇ = zeros(Complex{FT}, n_wann, n_wann)
    MAᵏᵇ = zeros(Complex{FT}, n_bands, n_wann)

    for ik = 1:n_kpts
        for ib = 1:n_bvecs

            ikpb = kpb_k[ib, ik]

            MAᵏᵇ .= overlap(M, kpb_k, ik, ikpb) * A[:, :, ikpb]
            Nᵏᵇ .= A[:, :, ik]' * MAᵏᵇ
            b .= recip_lattice * (kpoints[:, ikpb] + kpb_b[:, ib, ik] - kpoints[:, ik])

            w_ib = wb[ib]

            for m = 1:n_wann
                for n = 1:n_wann
                    A[m, n, :, ik] += w_ib * Nᵏᵇ[m, n] * b

                    if m == n
                        A[m, n, :, ik] -= w_ib * b
                    end
                end
            end

        end
    end

    A *= im

    A
end

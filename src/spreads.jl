module Spreads

# using StaticArrays
import LinearAlgebra as LA
import ..Utilities: overlap

# computes the MV energy
# From MV: Omega = sum_n <r2>n - |<r>n|^2
# <r>n = -1/N sum_k,b wb b Im ln(Mnn,kb)
# <r2>n = 1/N sum_k,b wb [(1 - |Mnn,kb|^2) + (Im ln(Mnn,kb))^2]
struct Omega_res
    Ωtot::Float64
    ΩI::Float64
    ΩOD::Float64
    ΩD::Float64 # not implemented yet
    Ωtilde::Float64
    frozen_weight::Float64
    # num_wann
    spreads::Array{Float64,1}
    # 3 * num_wann
    centers::Array{Float64,2}
    # num_bands * num_bands * num_kpts
    gradient::Array{ComplexF64,3}
    # fix_centers :: Array{Float64,2} #3 x nwannier
end

imaglog(z) = atan(imag(z), real(z))

@views function omega(params, A, compute_grad=false, only_r2=false)
    mu = 0.
    nfrozen = 0 # keep in case we want to do this later on
    if compute_grad
        if !only_r2
            fix_center = false # true
            if fix_center
                fix_centers = [ # bohr
                2.55000   2.55000   2.55000; 
                2.55000   2.55000   2.55000; 
                2.55000   2.55000   2.55000;
                2.55000   2.55000   2.55000;
                0.00000   0.00000   0.00000;
                0.00000   0.00000   0.00000;
                0.00000   0.00000   0.00000;
                0.00000   0.00000   0.00000;
                ]'
                centers = fix_centers
            else
                centers = omega(params, A, false).centers
            end
        end
    end
    grad = zeros(ComplexF64, params.num_bands, params.num_wann, params.num_kpts)
    r = zeros(3, params.num_wann)
    r2 = zeros(params.num_wann)
    ΩI = 0.
    ΩOD = 0.
    ΩD = 0.
    frozen_weight = 0.
    R = zeros(ComplexF64, params.num_bands, params.num_wann)
    T = zeros(ComplexF64, params.num_bands, params.num_wann)
    b = zeros(3)
    Mkb = zeros(ComplexF64, params.num_wann, params.num_wann)
    Okb = zeros(ComplexF64, params.num_bands, params.num_wann)
    for ik = 1:params.num_kpts
        frozen_weight -= mu * sum(abs2, A[1:nfrozen,:,ik])
        if compute_grad
            grad[1:nfrozen,:,ik] = -2 * mu * A[1:nfrozen,:,ik]
        end
        
        for ib = 1:params.num_bvecs
            ikpb = params.kpbs[ib, ik]
            # Mkb = overlap_A([i,j,k],[i_n,j_n,k_n],p)
            Okb .= overlap(params, ik, ikpb) * A[:,:,ikpb]
            # @assert overlap(params,ik,ikpb)' ≈ overlap(params,ikpb,ik) #compute-intensive, but should be true
            Mkb .= A[:,:,ik]' * Okb
            b .= params.recip_cell * (params.kpts[:,ikpb] .+ params.kpbs_disp[:,ib,ik] .- params.kpts[:,ik])

            ib_weight = params.kpbs_weight[ib,ik]
            if compute_grad
                # #MV way
                # fA(B) = (B-B')/2
                # fS(B) = (B+B')/(2*im)
                # q = imaglog.(diag(Mkb)) + centers'*b
                # for m=1:p.nwannier,n=1:p.nwannier
                #     R[m,n] = Mkb[m,n]*conj(Mkb[n,n])
                #     T[m,n] = Mkb[m,n]/Mkb[n,n]*q[n]
                # end
                # grad[:,:,i,j,k] += 4*p.wb*(fA(R) .- fS(T))


                q = imaglog.(LA.diag(Mkb))
                if !only_r2
                    q += centers' * b
                end
                for n = 1:params.num_wann
                    # error if division by zero. Should not happen if the initial gauge is not too bad
                    if abs(Mkb[n,n]) < 1e-10
                        println("Mkbnn too small! $ik -> $ikpb")
                        display(Mkb)
                        error()
                    end
                    @assert abs(Mkb[n,n]) > 1e-10
                    Tfac = -im * q[n] / Mkb[n,n]
                    for m = 1:params.num_bands
                        R[m,n] = -Okb[m,n] * conj(Mkb[n,n])
                        # T[m,n] = -im*Okb[m,n]/(Mkb[n,n])*q[n]
                        T[m,n] = Tfac * Okb[m,n]
                    end
                end
                grad[:,:,ik] .+= 4 * ib_weight .* (R .+ T)
            end

            ΩI += ib_weight * (params.num_wann - sum(abs2, Mkb))
            ΩOD += ib_weight * sum(abs2, Mkb .- LA.diagm(0 => LA.diag(Mkb)))
            for n = 1:params.num_wann
                if !only_r2
                    r[:,n] -= ib_weight* imaglog(Mkb[n,n]) * b
                end
                r2[n] += ib_weight * (1 - abs(Mkb[n,n])^2 + imaglog(Mkb[n,n])^2)
                # r2[n] += p.wb*2*(1 - real(Mkb[n,n]))
            end
        end
    end
    r /= params.num_kpts
    r2 /= params.num_kpts
    ΩI /= params.num_kpts
    ΩOD /= params.num_kpts
    ΩD /= params.num_kpts
    frozen_weight /= params.num_kpts
    grad /= params.num_kpts

    # @debug "Spreads" r r2' ΩI ΩOD ΩD

    spreads = r2 - dropdims(sum(abs.(r).^2, dims=1), dims=1)
    Ωtot = sum(spreads) + frozen_weight
    Ωtilde = Ωtot - ΩI

    return Omega_res(Ωtot, ΩI, ΩOD, ΩD, Ωtilde, frozen_weight, spreads, r, grad)
end

# local part of the contribution to r^2
function omega_loc(p, A)
    r = zeros(3, p.nwannier)
    r2 = zeros(p.nwannier)
    loc = zeros(Float64, p.N1, p.N2, p.N3)
    for i in 1:p.N1,j in 1:p.N2,k in 1:p.N3
        K = p.ijk_to_K[i,j,k]
        for ib = 1:p.nntot
            neighbor = p.neighbors[K,ib]
            i_n, j_n, k_n = p.K_to_ijk[neighbor,:]
            # Mkb = overlap_A([i,j,k],[i_n,j_n,k_n],p)
            Okb = overlap(K, neighbor, p) * A[:,:,i_n,j_n,k_n]
            Mkb = A[:,:,i,j,k]' * Okb
            for n = 1:p.nwannier
                loc[i,j,k] += p.wb * (1 - abs(Mkb[n,n])^2 + imaglog(Mkb[n,n])^2)
            end
        end
    end
    return loc
end

end

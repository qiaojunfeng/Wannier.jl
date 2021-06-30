using StaticArrays

# computes the MV energy
# From MV: Omega = sum_n <r2>n - |<r>n|^2
# <r>n = -1/N sum_k,b wb b Im ln(Mnn,kb)
# <r2>n = 1/N sum_k,b wb [(1 - |Mnn,kb|^2) + (Im ln(Mnn,kb))^2]
struct Omega_res
    Ωtot :: Float64
    ΩI :: Float64
    ΩOD :: Float64
    ΩD :: Float64 #not implemented yet
    Ωtilde :: Float64
    frozen_weight :: Float64
    spreads :: Array{Float64,1} #nband
    centers :: Array{Float64,2} #3 x nband
    gradient :: Array{ComplexF64, 5}#nband x nband x n1 x n2 x n3
    # fix_centers :: Array{Float64,2} #3 x nwannier
end

imaglog(z) = atan(imag(z), real(z))

@views function omega(p,A,compute_grad=false,only_r2=false)
    mu=0.
    nfrozen=0 #keep in case we want to do this later on
    if compute_grad
        if !only_r2
            fix_center = true
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
            centers = omega(p,A,false).centers
            end
        end
    end
    grad = zeros(ComplexF64,p.nband,p.nwannier,p.N1,p.N2,p.N3)
    r = zeros(3,p.nwannier)
    r2 = zeros(p.nwannier)
    ΩI = 0.
    ΩOD = 0.
    ΩD = 0.
    frozen_weight = 0.
    R = zeros(ComplexF64,p.nband,p.nwannier)
    T = zeros(ComplexF64,p.nband,p.nwannier)
    b = zeros(3)
    Mkb = zeros(ComplexF64,p.nwannier,p.nwannier)
    Okb = zeros(ComplexF64,p.nband,p.nwannier)
    for i=1:p.N1,j=1:p.N2,k=1:p.N3
        K = p.ijk_to_K[i,j,k]

        frozen_weight -= mu*sum(abs2, A[1:nfrozen,:,i,j,k])
        if compute_grad
            grad[1:nfrozen,:,i,j,k] = -2*mu*A[1:nfrozen,:,i,j,k]
        end
        
        for ib = 1:p.nntot
            neighbor = p.neighbors[K,ib]
            i_n,j_n,k_n = p.K_to_ijk[neighbor,:]
            # Mkb = overlap_A([i,j,k],[i_n,j_n,k_n],p)
            Okb .= overlap(K,neighbor,p)*A[:,:,i_n,j_n,k_n]
            # @assert overlap(K,neighbor,p)' ≈ overlap(neighbor,K,p) #compute-intensive, but should be true
            Mkb .= A[:,:,i,j,k]'*Okb
            b .= p.recp_unit_cell*([p.t1[i_n];p.t2[j_n];p.t3[k_n]] .+ p.displ_vecs[:,K,ib] .- [p.t1[i];p.t2[j];p.t3[k]])

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


                q = imaglog.(diag(Mkb))
                if !only_r2
                    q += centers'*b
                end
                for n=1:p.nwannier
                    if abs(Mkb[n,n]) < 1e-10 # error if division by zero. Should not happen if the initial gauge is not too bad
                        println("Mkbnn too large! $j $k -> $(j_n) $(k_n)")
                        display(Mkb)
                        error()
                    end
                    @assert abs(Mkb[n,n]) > 1e-10
                    Tfac = -im*q[n]/Mkb[n,n]
                    for m=1:p.nband
                        R[m,n] = -Okb[m,n]*conj(Mkb[n,n])
                        # T[m,n] = -im*Okb[m,n]/(Mkb[n,n])*q[n]
                        T[m,n] = Tfac*Okb[m,n]
                    end
                end
                grad[:,:,i,j,k] .+= 4*p.wb.*(R.+T)
            end

            ΩI += p.wb*(p.nwannier - sum(abs2,Mkb))
            ΩOD += p.wb*sum(abs2,Mkb .- diagm(0=>diag(Mkb)))
            for n=1:p.nwannier
                if !only_r2
                    r[:,n] -= p.wb*imaglog(Mkb[n,n])*b
                end
                r2[n] += p.wb*(1-abs(Mkb[n,n])^2 + imaglog(Mkb[n,n])^2)
                # r2[n] += p.wb*2*(1 - real(Mkb[n,n]))
            end
        end
    end
    Ntot = (p.N1*p.N2*p.N3)
    r /= Ntot
    r2 /= Ntot
    ΩI /= Ntot
    ΩOD /= Ntot
    ΩD /= Ntot
    frozen_weight /= Ntot
    (grad /= Ntot)

    spreads = (r2 - dropdims(sum(abs.(r).^2,dims=1),dims=1))
    Ωtot = sum(spreads) + frozen_weight
    Ωtilde = Ωtot - ΩI
    return Omega_res(Ωtot,ΩI,ΩOD,ΩD,Ωtilde,frozen_weight,spreads,r,grad)
end

# local part of the contribution to r^2
function omega_loc(p,A)
    r = zeros(3,p.nwannier)
    r2 = zeros(p.nwannier)
    loc = zeros(Float64,p.N1,p.N2,p.N3)
    for i=1:p.N1,j=1:p.N2,k=1:p.N3
        K = p.ijk_to_K[i,j,k]
        for ib = 1:p.nntot
            neighbor = p.neighbors[K,ib]
            i_n,j_n,k_n = p.K_to_ijk[neighbor,:]
            # Mkb = overlap_A([i,j,k],[i_n,j_n,k_n],p)
            Okb = overlap(K,neighbor,p)*A[:,:,i_n,j_n,k_n]
            Mkb = A[:,:,i,j,k]'*Okb
            for n=1:p.nwannier
                loc[i,j,k] += p.wb*(1-abs(Mkb[n,n])^2 + imaglog(Mkb[n,n])^2)
            end
        end
    end
    return loc
end

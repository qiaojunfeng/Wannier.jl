# include("../../wannierize_utils.jl")

#=
We look at the free Hamiltonian H = -Delta in the unit cell 2pi x [-.5,.5]^d, with BZ [0,1]^d
For every k, we therefore solve the Hamiltonian H_k (-i nabla + k)^2 on 2pi [-.5,.5]^d with periodic BC
This is discretized in the basis of K/G vectors e_K(r) = exp(iK.r), with K integers
In this basis, H_KK' = delta_KK' |K|^2
Therefore unk(r) = e^(i K_n r), where K_n corresponds to the n-th element of the K ordered by increasing value
if perm = sortperm(eig_k), then perm gives the mapping from n to K
Ki is limited to -Li:Li. By Fourier duality this gives a real-space grid of the unit cell


i,j,k are indices in the BZ (from 1 to N1,N2,N3), kind is the flat version of (i,j,k) (from 1 to N1*N2*N3).
The BZ loops are always on [0,1]^d, which is shifted by "corner" when computing λ, so that the actual BZ is corner + [0,1]^d
=#

using SpecialFunctions
using LinearAlgebra
using Dates

nwannier = 2
# k-grid
N1 = 1
N2 = 1
N3 = 20

# G-grid
L1 = 0
L2 = 0
L3 = 2

dim = 3 - count(x->x==1, [N1,N2,N3])

#SCDM parameters
mu = 0.0;
sigma = 2.0;

do_shift = false
corner = [0.,0.,0.]
# corner = [0.,0.,0.] # do not use directly in the code, only the functions below

# all loops in the code refer to a k point grid like (0:N-1)/N. This transfers to the actual kpt grid
function k_red_to_real(k)::Vector{Float64} #uses non-constant globals, so use explicit type
    if do_shift
        shift = [1/N1/2, 1/N2/2, 1/N3/2]
        shift[1:3-dim] = 0 #don't shift singleton dimensions
        k+shift+corner
    else
        k+corner
    end
end
function k_real_to_red(k)::Vector{Float64}
    if do_shift
        shift = [1/N1/2, 1/N2/2, 1/N3/2]
        shift[1:3-dim] = 0 #don't shift singleton dimensions
        k-shift-corner
    else
        k-corner
    end
end

k = randn(3); @assert k ≈ k_red_to_real(k_real_to_red(k)) #sanity check

λ(k,Ks) = dropdims(sum(abs2, (Ks .+ k_red_to_real(k)), dims=1),dims=1) #eigenvalues at k in [0,1]^d


nband = (2L1+1)*(2L2+1)*(2L3+1) #size of the Hamiltonian
Ks = zeros(Int,3,nband)
let ind = 1
    for l1=-L1:L1
        for l2=-L2:L2
            for l3=-L3:L3
                Ks[:,ind] = [l1, l2, l3]
                ind += 1
            end
        end
    end
end
function write_files(N1, N2, N3, L1, L2, L3, nwannier)
    nband = (2L1+1)*(2L2+1)*(2L3+1) #size of the Hamiltonian
    Ks = zeros(Int,3,nband) #list of K vectors
    Ntot = N1*N2*N3
    Ns = [N1,N2,N3]
    R = zeros(ComplexF64,nband,3)
    Rfrac = zeros(ComplexF64,nband,3)
    ind = 1
    for l1=-L1:L1
        for l2=-L2:L2
            for l3=-L3:L3
                Ks[:,ind] = [l1, l2, l3]
                R[ind,:]  = 2*pi*[l1/(2*L1+1), l2/(2*L2+1), l3/(2*L3+1)]
                Rfrac[ind,:] = [(l1+L1)/(2*L1+1), (l2+L2)/(2*L2+1), (l3+L3)/(2*L3+1)]
                ind += 1
            end
        end
    end

    win = open("free.win","w")
    write(win,"""
num_bands = $nband
num_wann = $nwannier

begin unit_cell_cart
bohr
$(2pi) 0 0
0 $(2pi) 0
0 0 $(2pi)
end unit_cell_cart

mp_grid = $N1 $N2 $N3

begin kpoints
""")
    kind = 1
    kind_to_ijk = zeros(Int64,Ntot,3) #kind runs from 1:Ntot, i,j,k from 1:N1,N2,N3
    ijk_to_kind = zeros(Int64,N1,N2,N3)
    kpt_order = [1,2,3] #potentially allow different ordering
    for i=0:N1-1
        for j=0:N2-1
            for k=0:N3-1
                write(win, "$(i/N1) $(j/N2) $(k/N3)\n")

                kind += 1
                ijk = [i,j,k]
                big_offset = Ns[kpt_order[1]]*Ns[kpt_order[2]]
                med_offset = Ns[kpt_order[2]]
                kind = ijk[kpt_order[1]]*big_offset + ijk[kpt_order[2]]*med_offset + ijk[kpt_order[3]]+1
                ijk_to_kind[i+1,j+1,k+1] = kind
                kind_to_ijk[kind,:] = [i+1 j+1 k+1]

            end
        end
    end

    write(win,"end kpoints")
    close(win)

    neighbors =
        [1 0 0;
         0 1 0;
         0 0 1;
         -1 0 0;
         0 -1 0;
         0 0 -1]

    # overlap matrix Mmnkb between k and k+b, with k+b = displ_vec + kpb
    # <un(k)|um(kpb + displ_vec)> = <un(k)|e^-(displ_vec)x um(kpb)>
    # = delta_KK', but mn are sorted versions of K K'
    function get_ovl(k, kpb, displ_vec,Ks)
        eig_k = λ(k,Ks)
        eig_kpb = λ(kpb,Ks)

        sort_k = sortperm(eig_k)
        sort_kpb = sortperm(eig_kpb)
        A = zeros(Int64, size(Ks,2),size(Ks,2))
        for m=1:nband
            for n=1:nband
                Km = sort_k[m]
                @assert eig_k[Km] == sort(eig_k)[m]
                Kpbn = sort_kpb[n]
                A[m,n] = Ks[:,Km] == (Ks[:,Kpbn] .- displ_vec)
            end
        end
        return A
    end

    eig = open("free.eig","w")
    eig_allk = zeros(ComplexF64,nband,N1*N2*N3)
    ind = 1
    for i=0:N1-1
        for j=0:N2-1
            for k=0:N3-1
                ijk = [i,j,k]
                kind = ijk_to_kind[i+1,j+1,k+1]
                eig_k = sort(λ(ijk./Ns,Ks))
                eig_allk[:,ind] = eig_k
                ind += 1
                for iband=1:nband
                    write(eig,"$iband $kind $(eig_k[iband])\n")
                end
            end
        end
    end
    close(eig)

    mmn = open("free.mmn","w")
    write(mmn, "Created by free.jl $L1 $L2 $L3\n")
    nneighbors = count(x -> x != 1, Ns)*2
    write(mmn, "$nband $Ntot $nneighbors\n")
    for i=0:N1-1
        for j=0:N2-1
            for k=0:N3-1
                ijk = [i,j,k]
                for ib = 1:6
                    b = neighbors[ib,:]
                    ijkpb_recen = mod.(ijk+b,[N1; N2; N3]) #bring back to BZ
                    kind = ijk_to_kind[(ijk.+1)...]
                    kindb = ijk_to_kind[(ijkpb_recen.+1)...]
                    if kind == kindb
                        continue #ignore singleton dimensions
                    end
                    displ_vec = div.((ijk+b) - ijkpb_recen,Ns)
                    write(mmn,"$kind $kindb $(displ_vec[1]) $(displ_vec[2]) $(displ_vec[3])\n")
                    ovl = get_ovl(ijk./Ns,(ijkpb_recen)./Ns,displ_vec,Ks)
                    for n=1:nband
                        for m=1:nband
                            write(mmn,string(ovl[m,n]), " 0\n")
                        end
                    end
                end
            end
        end
    end
    close(mmn)

    # Unk and SCDM stuff
    # same for all kpoints
    Amn = zeros(ComplexF64,N1,N2,N3,nband,nwannier)

    gamma = [0.0,0.0,0.0]
    eig_k = λ(k_real_to_red(gamma),Ks) # λ is computed in reduced coordinates
    Unk = exp.(1im*R*Ks)
    occ = (1/2)*erfc.((real(eig_k).-mu)./sigma)

    Qtemp, Rtemp, piv = qr(diagm(0=>occ)*Unk',Val(true))

    cols = piv[1:nwannier]
    min_S = Inf
    for i=0:N1-1
        for j=0:N2-1
            for k=0:N3-1
                kpt = [i,j,k]./Ns

                phase = exp.(1im*R[cols,:]*k_red_to_real(kpt))

                eig_k = λ(kpt,Ks)
                n_to_K = sortperm(eig_k)
                K_to_n = invperm(n_to_K)
                eig_k = eig_k[n_to_K]
                # Note that the order of Ks is kpt dependent
                Unk = zeros(ComplexF64,nband,nband)
                for m = 1 : nband
                    Km = n_to_K[m]
                    Unk[:,m] = exp.(1im*R*(Ks[:,Km]))
                end
                occ = (1/2)*erfc.((real(eig_k).-mu)./sigma)

                Y = diagm(0=>phase)*Unk[cols,:]*diagm(0=>occ)
                Y = Y'

                U,S,V = svd(Y)
                Amn[i+1,j+1,k+1,:,:] = U*V'
                min_S = min(min_S, minimum(S))
            end
        end
    end
    println("Min singular value: $min_S")

    for file in ("free.amn", "free.SCDM.amn")
        out = open(file,"w")
        write(out, "Created by free.jl ", string(Dates.now()), "\n")
        write(out, "$(nband) $(N1*N2*N3) $(nwannier)\n")
        for i=1:N1
            for j=1:N2
                for k=1:N3
                    kind = ijk_to_kind[i,j,k]
                    for n=1:nwannier
                        for m = 1:nband
                            coeff = Amn[i,j,k,m,n]
                            write(out, "$m $n $kind $(real(coeff)) $(imag(coeff))\n")
                        end
                    end
                end
            end
        end
        close(out)
    end
end

write_files(N1,N2,N3,L1,L2,L3,nwannier)

include("../../run_optim.jl")

loc = omega_loc(p,A)
# writedlm("data",loc)
println(maximum(loc))
println("Max loc = ", maximum(loc))
include("plot_free_wannier.jl")

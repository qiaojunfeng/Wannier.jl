export symmetrize, fbz2ibz, unfold_amn, find, rotate, get_symop

seedname = "/scratch/gcistaro/Calculations/GaAs_SAWF/GaAs"
#

function symmetrize(Op, isym, Rn)
    println("Reading isym file")

    #get R for the factor in D
    n_wann = size(isym.gk[1].h[1])[1]
    nk_ibz = size(isym.kpoints)[1]
    
    
    println("nwann=", n_wann, "  nk_ibz", nk_ibz)

    Op_sym = [zeros(ComplexF64, n_wann, n_wann) for _ in 1:nk_ibz]
    

    for ik in 1:nk_ibz
        for ih in 1:size(isym.gk[ik].h)[1]
            # D[isym,m,n] = <g_m| S^-1 |g_n>
            # Op_h = h[ik]*Op*D[h^-1][k]
            ih_index = isym.gk[ik].isym[ih]
            phase = [exp(-im*2.0*π*dot(isym.kpoints[ik],Rn[ih_index][iwann][:])) for iwann in 1:n_wann]
            Op_h = apply_sym_Uk(Op[ik], phase, isym.D[ih_index], isym.symops[ih_index].t_rev, isym.gk[ik].h[ih])
            #if ik == 2
            #   println("isym=", ih_index-1, "    t_rev=", isym.symops[ih_index].t_rev)
            #   #println("rn", Rn[ih_index][:][:])
            #   #println("D")
            #    #for iwann in 1:4
            #    #    println(isym.D[ih_index][iwann,:])
            #    #end 
            #    #println("Op")
            #    #for iwann in 1:4
            #    #    println(Op[ik][iwann,:])
            #    #end 
            #    #println("d")
            #    #for iwann in 1:4
            #    #    println(isym.gk[ik].h[ih][iwann,:])
            #    #end 
            #    #println("tot")
            #    dpy = zeros(ComplexF64, n_wann, n_wann)
            #    open("/scratch/gcistaro/codes/symWannier/d"*string(ih_index-1)*".txt") do io
            #        for iwann in 1:n_wann
            #            line = split(readline(io))
            #            for jwann in 1:n_wann
            #                dpy[iwann,jwann] = parse.(Float64, line[jwann]) +im*parse.(Float64, line[n_wann+jwann])
            #            end
            #        end
            #    end
            #    #if(norm(dpy-Op_h)>1.e-06)
            #    if(norm(dpy-isym.gk[ik].h[ih])>1.e-06)
            #    for iwann in 1:4
            #            println(isym.gk[ik].h[ih][iwann, :])
            #            println(dpy[iwann,:])
            #        end
            #    end
            #end
            
            Op_sym[ik] = Op_sym[ik] + Op_h
        end
        Op_sym[ik][:,:] ./= size(isym.gk[ik].h)[1]
        #println("Op_sym:::")
        #for iwann in 1:4
        #    println(Op_sym[ik][iwann,:])
        #end
    end

    return Op_sym
end

function apply_sym_Uk(Op, phase, D, t_rev, d=nothing)
    Op_h = Op*(phase.*D)
    if t_rev == 1
        Op_h = conj.(Op_h)
    end
    if d != nothing
        Op_h = d*Op_h
    end
    return Op_h 
end


function unfold_amn(A, f2i, Rn, isym)
    num_kpoints = size(f2i)[1]
    n_wann      = size(A[1])[1]
    Afull = [zeros(ComplexF64, size(A[1])) for _ in 1:num_kpoints]
    for ik in 1:num_kpoints
        ik_ibz      = f2i[ik][1][1]
        sym_index   = f2i[ik][1][2]
        phase = [exp(-im*2.0*π*dot(isym.kpoints[ik_ibz],Rn[sym_index][iwann][:])) for iwann in 1:n_wann]
        Afull[ik] = apply_sym_Uk(A[ik_ibz], phase, isym.D[sym_index],isym.symops[sym_index].t_rev)
    end
    return Afull
end

function find_pos(symops,proj,D)
    n_wann = size(proj)[1]
    n_symm = size(symops)[1]

    #we find where D sends the wannier center
    R = [[zeros(3) for _ in 1: n_wann] for _ in 1:n_symm]
    for isym in 1:n_symm
        for iwann in 1:n_wann
            #calculate rotated center s*center
            Sw = proj[iwann].center'*symops[isym].s - symops[isym].ft'
            #evaluate which of the SAWF we are jumping into
            iRw = findall(>(1.e-07), abs.(D[isym][:, iwann])   )         
            #println("abs.(D[isym][:, iwann]) ", abs.(D[isym][:, iwann]))
            #println(iRw)
            @assert all(x -> x == proj[iRw[1]].center, [proj[iRw[iwann_]].center for iwann_ in 1:size(iRw)[1] ])
            R[isym][iwann][:] = Sw' - proj[iRw[1]].center 
            #println(Sw, proj[iRw[1]].center, iRw[1], R[isym][iwann][:])
        end
    end
    return R
end


function rotate(kvector, symops)
    sk = symops.s*kvector
    if symops.t_rev == 1
        sk = -sk
    end 
    return sk
end

function find(kvector, kset)
    i = findall(<(1.e-05), [norm(kvector -k -round.(kvector -k)) for k in kset])
    return i[1]
end

function find_strict(kvector, kset)
    i = findall(<(1.e-05), [norm(kvector -k ) for k in kset])
    return i[1]
end

function fbz2ibz(k_fbz, k_ibz, symops)
    i = []
    for kpt in k_fbz
        push!(i,[])
    end

    for ik in 1:size(k_ibz)[1]
        for is in 1:size(symops)[1]
            ifbz = find(rotate(k_ibz[ik], symops[is]), k_fbz)
            push!(i[ifbz], [ik,is])
        end 
    end 
    return i
end


function get_symop(symops, symlist)
    # initialize everything to identity
    s0    = Matrix(1.0I, 3, 3)
    ft0   = zeros(3,1)
    t_rev = 0 

    for info in symlist
        isym  = info[1]
        inv   = info[2]
        ft1   = zeros(3,1)

        if inv == -1
            #inv(rs-t) = (r+t)s^-1 = rs^-1 +ts^-1
            s1  = symops[symops[isym].invs].s
            ft1 = ((symops[isym].ft')*s1)'
        else
            #rs-t
            s1  = symops[isym].s
            ft1 = symops[isym].ft
        end

        #now we merge operation 0 and 1
        #r' = r*s0 - t0
        #r''= r'*s1- t1 = (r*s0 - t0)*s1 - t1 = r*s0*s1 - t0*s1 - t1
        s0 = s0*s1
        ft0 = (((ft0')*s1)' + ft1)
        t_rev = mod( t_rev + symops[isym].t_rev, 2 )
    end

    #find operation in symops
    for isym in 1:size(symops)[1]
        symop = symops[isym]
        tDelta = ft0 - symop.ft
        if isapprox( s0, symop.s, atol = 1.e-06) && isapprox( tDelta, round.(tDelta), atol = 1.e-06 ) && ( t_rev == symop.t_rev )
            return isym, tDelta
        end
    end
end
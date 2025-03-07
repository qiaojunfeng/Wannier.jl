export write_full_data, to_crys, unfold_mmn


function to_crys(nnkp, bvectors)
    bvectors_cryst = [zeros(3) for _ in 1:size(bvectors)[1]]
    for ib in 1:size(bvectors)[1]
        bvectors_cryst[ib] = inv(nnkp.recip_lattice)*bvectors[ib]
    end
    return bvectors_cryst
end 



function write_full_data(seedname)
    nnkp = read_nnkp(seedname*".nnkp")
    Kstencil = read_nnkp_compute_bweights(seedname*".nnkp")
    isym = read_w90_isym(seedname*".isym")
    rescale_repmat(isym.gk) #not sure why this is needed (to check)
    f2i = fbz2ibz(nnkp.kpoints, isym.kpoints, isym.symops)
    println(f2i)

    #Eig
    ieig = read_eig(seedname*".ieig")
    eig  = unfold_eig(ieig, f2i)

    #Amn
    A  = read_amn(seedname*".iamn")
    #fill vector for the factor e^{-ikR_n'} appearing in Koretsune (9) 
    Rn = find_pos(isym.symops,  
                nnkp.projections,
                isym.D
    )
    println(Rn)
    Asymm = symmetrize(A, isym, Rn)
    Afull = unfold_amn(Asymm, f2i, Rn, isym)

    #Mmn
    M, kpb_k, kpb_G = read_mmn(seedname*".immn")
    bvectors_cryst = to_crys(nnkp, Kstencil.bvectors)
    Mfull, kpb_kfull, kpb_Gfull = unfold_mmn(M, kpb_k, kpb_G, nnkp, isym, f2i, bvectors_cryst)


    #check Eig with symwannier
    eigpy = read_eig(seedname*".eig")
    for ik in 1:size(f2i)[1]
        i = findall( >(1.e-10), abs.(eigpy[ik]-eig[ik]))
        println("ik = ", ik , "    " , i)
    end 



    #check Amn with symwannier
    Apy = read_amn(seedname*".amn")
    for ik in 1:size(f2i)[1]
        i = findall( >(1.e-10), abs.(Afull[ik]-Apy[ik]))
        println("ik = ", ik , "    " , i)
    end 



    #check Mmn with symwannier
    Mpy, kpb_k, kpb_G = read_mmn(seedname*".mmn")
    for ik in 1:size(f2i)[1]
        for ib in 1:8
            i = findall( >(1.e-10), abs.(Mfull[ik][ib]-Mpy[ik][ib]))
            println("ik = ", ik , "    " ,i)
        end
    end 

end



function unfold_mmn(M, kpb_k, kpb_G, nnkp, isym, f2i, bvectors_cryst)
#    pydata = [[zeros(Int32, 6) for ib in 1:8] for ik in 1:8]
#    pydata[ 1 ][ 1 ] = [ 1 ,  2 ,  2 ,  1 ,  1 ,  2 ,  1 ]
#    pydata[ 1 ][ 2 ] = [ 2 ,  3 ,  2 ,  1 ,  3 ,  3 ,  3 ]
#    pydata[ 1 ][ 3 ] = [ 3 ,  5 ,  2 ,  1 ,  4 ,  5 ,  4 ]
#    pydata[ 1 ][ 4 ] = [ 4 ,  8 ,  2 ,  1 ,  2 ,  8 ,  2 ]
#    pydata[ 1 ][ 5 ] = [ 5 ,  2 ,  2 ,  1 ,  1 ,  2 ,  1 ]
#    pydata[ 1 ][ 6 ] = [ 6 ,  3 ,  2 ,  1 ,  3 ,  3 ,  3 ]
#    pydata[ 1 ][ 7 ] = [ 7 ,  5 ,  2 ,  1 ,  4 ,  5 ,  4 ]
#    pydata[ 1 ][ 8 ] = [ 8 ,  8 ,  2 ,  1 ,  2 ,  8 ,  2 ]
#    pydata[ 2 ][ 1 ] = [ 1 ,  1 ,  1 ,  1 ,  1 ,  1 ,  1 ]
#    pydata[ 2 ][ 2 ] = [ 2 ,  4 ,  3 ,  1 ,  1 ,  4 ,  1 ]
#    pydata[ 2 ][ 3 ] = [ 3 ,  6 ,  3 ,  1 ,  9 ,  6 ,  9 ]
#    pydata[ 2 ][ 4 ] = [ 4 ,  7 ,  3 ,  1 ,  5 ,  7 ,  5 ]
#    pydata[ 2 ][ 5 ] = [ 5 ,  1 ,  1 ,  1 ,  1 ,  1 ,  1 ]
#    pydata[ 2 ][ 6 ] = [ 6 ,  4 ,  3 ,  1 ,  1 ,  4 ,  1 ]
#    pydata[ 2 ][ 7 ] = [ 7 ,  6 ,  3 ,  1 ,  9 ,  6 ,  9 ]
#    pydata[ 2 ][ 8 ] = [ 8 ,  7 ,  3 ,  1 ,  5 ,  7 ,  5 ]
#    pydata[ 3 ][ 1 ] = [ 2 ,  4 ,  3 ,  3 ,  1 ,  4 ,  1 ]
#    pydata[ 3 ][ 2 ] = [ 1 ,  1 ,  1 ,  3 ,  1 ,  1 ,  1 ]
#    pydata[ 3 ][ 3 ] = [ 8 ,  7 ,  3 ,  3 ,  5 ,  7 ,  5 ]
#    pydata[ 3 ][ 4 ] = [ 7 ,  6 ,  3 ,  3 ,  9 ,  6 ,  9 ]
#    pydata[ 3 ][ 5 ] = [ 6 ,  4 ,  3 ,  3 ,  1 ,  4 ,  1 ]
#    pydata[ 3 ][ 6 ] = [ 5 ,  1 ,  1 ,  3 ,  1 ,  1 ,  1 ]
#    pydata[ 3 ][ 7 ] = [ 4 ,  7 ,  3 ,  3 ,  5 ,  7 ,  5 ]
#    pydata[ 3 ][ 8 ] = [ 3 ,  6 ,  3 ,  3 ,  9 ,  6 ,  9 ]
#    pydata[ 4 ][ 1 ] = [ 1 ,  3 ,  2 ,  1 ,  3 ,  3 ,  3 ]
#    pydata[ 4 ][ 2 ] = [ 2 ,  2 ,  2 ,  1 ,  1 ,  2 ,  1 ]
#    pydata[ 4 ][ 3 ] = [ 3 ,  8 ,  2 ,  1 ,  2 ,  8 ,  2 ]
#    pydata[ 4 ][ 4 ] = [ 4 ,  5 ,  2 ,  1 ,  4 ,  5 ,  4 ]
#    pydata[ 4 ][ 5 ] = [ 5 ,  3 ,  2 ,  1 ,  3 ,  3 ,  3 ]
#    pydata[ 4 ][ 6 ] = [ 6 ,  2 ,  2 ,  1 ,  1 ,  2 ,  1 ]
#    pydata[ 4 ][ 7 ] = [ 7 ,  8 ,  2 ,  1 ,  2 ,  8 ,  2 ]
#    pydata[ 4 ][ 8 ] = [ 8 ,  5 ,  2 ,  1 ,  4 ,  5 ,  4 ]
#    pydata[ 5 ][ 1 ] = [ 3 ,  6 ,  3 ,  4 ,  9 ,  6 ,  9 ]
#    pydata[ 5 ][ 2 ] = [ 8 ,  7 ,  3 ,  4 ,  5 ,  7 ,  5 ]
#    pydata[ 5 ][ 3 ] = [ 1 ,  1 ,  1 ,  4 ,  1 ,  1 ,  1 ]
#    pydata[ 5 ][ 4 ] = [ 6 ,  4 ,  3 ,  4 ,  1 ,  4 ,  1 ]
#    pydata[ 5 ][ 5 ] = [ 7 ,  6 ,  3 ,  4 ,  9 ,  6 ,  9 ]
#    pydata[ 5 ][ 6 ] = [ 4 ,  7 ,  3 ,  4 ,  5 ,  7 ,  5 ]
#    pydata[ 5 ][ 7 ] = [ 5 ,  1 ,  1 ,  4 ,  1 ,  1 ,  1 ]
#    pydata[ 5 ][ 8 ] = [ 2 ,  4 ,  3 ,  4 ,  1 ,  4 ,  1 ]
#    pydata[ 6 ][ 1 ] = [ 3 ,  8 ,  2 ,  9 ,  2 ,  5 ,  4 ]
#    pydata[ 6 ][ 2 ] = [ 2 ,  2 ,  2 ,  9 ,  1 ,  8 ,  2 ]
#    pydata[ 6 ][ 3 ] = [ 8 ,  5 ,  2 ,  9 ,  4 ,  2 ,  1 ]
#    pydata[ 6 ][ 4 ] = [ 5 ,  3 ,  2 ,  9 ,  3 ,  3 ,  3 ]
#    pydata[ 6 ][ 5 ] = [ 7 ,  8 ,  2 ,  9 ,  2 ,  5 ,  4 ]
#    pydata[ 6 ][ 6 ] = [ 6 ,  2 ,  2 ,  9 ,  1 ,  8 ,  2 ]
#    pydata[ 6 ][ 7 ] = [ 4 ,  5 ,  2 ,  9 ,  4 ,  2 ,  1 ]
#    pydata[ 6 ][ 8 ] = [ 1 ,  3 ,  2 ,  9 ,  3 ,  3 ,  3 ]
#    pydata[ 7 ][ 1 ] = [ 8 ,  5 ,  2 ,  5 ,  4 ,  8 ,  2 ]
#    pydata[ 7 ][ 2 ] = [ 2 ,  2 ,  2 ,  5 ,  1 ,  5 ,  4 ]
#    pydata[ 7 ][ 3 ] = [ 1 ,  3 ,  2 ,  5 ,  3 ,  3 ,  3 ]
#    pydata[ 7 ][ 4 ] = [ 7 ,  8 ,  2 ,  5 ,  2 ,  2 ,  1 ]
#    pydata[ 7 ][ 5 ] = [ 4 ,  5 ,  2 ,  5 ,  4 ,  8 ,  2 ]
#    pydata[ 7 ][ 6 ] = [ 6 ,  2 ,  2 ,  5 ,  1 ,  5 ,  4 ]
#    pydata[ 7 ][ 7 ] = [ 5 ,  3 ,  2 ,  5 ,  3 ,  3 ,  3 ]
#    pydata[ 7 ][ 8 ] = [ 3 ,  8 ,  2 ,  5 ,  2 ,  2 ,  1 ]
#    pydata[ 8 ][ 1 ] = [ 8 ,  7 ,  3 ,  2 ,  5 ,  7 ,  5 ]
#    pydata[ 8 ][ 2 ] = [ 3 ,  6 ,  3 ,  2 ,  9 ,  6 ,  9 ]
#    pydata[ 8 ][ 3 ] = [ 2 ,  4 ,  3 ,  2 ,  1 ,  4 ,  1 ]
#    pydata[ 8 ][ 4 ] = [ 5 ,  1 ,  1 ,  2 ,  1 ,  1 ,  1 ]
#    pydata[ 8 ][ 5 ] = [ 4 ,  7 ,  3 ,  2 ,  5 ,  7 ,  5 ]
#    pydata[ 8 ][ 6 ] = [ 7 ,  6 ,  3 ,  2 ,  9 ,  6 ,  9 ]
#    pydata[ 8 ][ 7 ] = [ 6 ,  4 ,  3 ,  2 ,  1 ,  4 ,  1 ]
#    pydata[ 8 ][ 8 ] = [ 1 ,  1 ,  1 ,  2 ,  1 ,  1 ,  1 ]
#
#sympy = [[zeros(5) for ib in 1:8] for ik in 1:8]
#
#sympy[ 1 ][ 1 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 1 ][ 2 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 1 ][ 3 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 1 ][ 4 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 1 ][ 5 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 1 ][ 6 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 1 ][ 7 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 1 ][ 8 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 2 ][ 1 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 2 ][ 2 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 2 ][ 3 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 2 ][ 4 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 2 ][ 5 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 2 ][ 6 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 2 ][ 7 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 2 ][ 8 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 3 ][ 1 ] =  [ 3 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 3 ][ 2 ] =  [ 3 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 3 ][ 3 ] =  [ 4 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 3 ][ 4 ] =  [ 2 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 3 ][ 5 ] =  [ 3 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 3 ][ 6 ] =  [ 3 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 3 ][ 7 ] =  [ 4 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 3 ][ 8 ] =  [ 2 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 4 ][ 1 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 4 ][ 2 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 4 ][ 3 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 4 ][ 4 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 4 ][ 5 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 4 ][ 6 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 4 ][ 7 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 4 ][ 8 ] =  [ 1 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 5 ][ 1 ] =  [ 3 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 5 ][ 2 ] =  [ 2 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 5 ][ 3 ] =  [ 4 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 5 ][ 4 ] =  [ 4 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 5 ][ 5 ] =  [ 3 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 5 ][ 6 ] =  [ 2 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 5 ][ 7 ] =  [ 4 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 5 ][ 8 ] =  [ 4 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 6 ][ 1 ] =  [ 8 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 6 ][ 2 ] =  [ 8 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 6 ][ 3 ] =  [ 8 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 6 ][ 4 ] =  [ 8 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 6 ][ 5 ] =  [ 8 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 6 ][ 6 ] =  [ 8 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 6 ][ 7 ] =  [ 8 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 6 ][ 8 ] =  [ 8 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 7 ][ 1 ] =  [ 10 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 7 ][ 2 ] =  [ 10 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 7 ][ 3 ] =  [ 10 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 7 ][ 4 ] =  [ 10 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 7 ][ 5 ] =  [ 10 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 7 ][ 6 ] =  [ 10 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 7 ][ 7 ] =  [ 10 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 7 ][ 8 ] =  [ 10 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 8 ][ 1 ] =  [ 3 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 8 ][ 2 ] =  [ 4 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 8 ][ 3 ] =  [ 2 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 8 ][ 4 ] =  [ 2 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 8 ][ 5 ] =  [ 3 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 8 ][ 6 ] =  [ 4 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 8 ][ 7 ] =  [ 2 , 1 , 0.0 , 0.0 , 0.0 ]
#sympy[ 8 ][ 8 ] =  [ 2 , 1 , 0.0 , 0.0 , 0.0 ]


    Mfull = [[zeros(ComplexF64, size(M[1][1])) for _ in 1:size(M[1])[1]] for _ in 1:size(nnkp.kpoints)[1]]
    kpb_kfull = [[0 for _ in 1:size(M[1])[1]] for _ in 1:size(nnkp.kpoints)[1]]
    kpb_Gfull = [[zeros(Float64, 3) for _ in 1:size(M[1])[1]] for _ in 1:size(nnkp.kpoints)[1]]

    println(size(M)[1], " ", size(M[1])[1])
    for ikf in 1:size(nnkp.kpoints)[1]
        kf       = nnkp.kpoints[ikf]
        iki      = f2i[ikf][1][1]
        isym_ikf = f2i[ikf][1][2]
        ki       = isym.kpoints[iki]
        
        for ibf in 1:size(bvectors_cryst)[1]
            # getting b point as (k+b) + G - k 
            #ikbf = nnkp.kpb_k[ikf][ibf] #index of k+b
            #G = nnkp.kpb_G[ikf][ibf]
            bf  = bvectors_cryst[ibf]#nnkp.kpoints[ikbf] + G - nnkp.kpoints[ikf] # unrotated b vector
            ikbf = find(kf + bf, nnkp.kpoints)#nnkp.kpb_k[ikf][ibf] #index of k+b
            #bi = h^{-1} bf
            bi = rotate(bf, isym.symops[isym.symops[isym_ikf].invs]) 
            ibi = find_strict(bi, bvectors_cryst)
            ikbi= find( ki + bi, nnkp.kpoints)

            kpb_kfull[ikf][ibf] = ikbf
            kpb_Gfull[ikf][ibf] = nnkp.kpoints[ikf] + bf - nnkp.kpoints[ikbf]  
            
            ikbi_ibz  = f2i[ikbi][1][1]            
            isym_ikbi = f2i[ikbi][1][2] #symmetry to get ki+bi in IBZ
            isym_ikbf = f2i[ikbf][1][2] #symmetry to get kf+bf in IBZ

            #if norm(pydata[ikf][ibf] - [ibi, ikbi, ikbi_ibz, isym_ikf, isym_ikbi, ikbf, isym_ikbf]) > 1.e-07
            #    println(ikf, " ", ibf)
            #    println(ibi, " ", ikbi, " " ,ikbi_ibz," ",isym_ikf,  " ", isym_ikbi, " " , ikbf, " " ,isym_ikbf)#isym_tot, " ", tDelta)
            #    println(pydata[ikf][ibf])
            #end




            #get equivalent operation of g^{-1}[isym_ikbi]*g^{-1}[isym_ikf]g[isym_ikbf] 
            isym_h, tDelta = get_symop(isym.symops, [[isym_ikbi, -1], [isym_ikf, -1], [isym_ikbf, +1]])

            #if norm(sympy[ikf][ibf] - [isym_h, 1., tDelta[1], tDelta[2], tDelta[3]]) > 1.e-07
            #    println(ikf, " ", ibf)
            #    println(isym_h, " ", 1.,  " ", tDelta[1], " ",tDelta[2], " ",tDelta[3])#isym_tot, " ", tDelta)
            #    println(sympy[ikf][ibf])
            #end

            factor = 1#changes with spinors

            index_isymh = findall(<(1.e-06), abs.(isym.gk[ikbi_ibz].isym.-isym_h))[1]
            if isym.symops[isym_ikbi].t_rev == 1
                Mfull[ikf][ibf] = M[iki][ibi]*conj.(isym.gk[ikbi_ibz].h[index_isymh])*factor
            else
                Mfull[ikf][ibf] = M[iki][ibi]*isym.gk[ikbi_ibz].h[index_isymh]*factor
            end

            if isym.symops[isym_ikf].t_rev == 1
                Mfull[ikf][ibf] = conj.(Mfull[ikf][ibf])
            end

            kbi_ibz = isym.kpoints[ikbi_ibz]
            Mfull[ikf][ibf] = Mfull[ikf][ibf]*exp(im*2.0*π*(dot(bi, isym.symops[isym_ikf].ft) + dot(kbi_ibz, tDelta)))
        end
    end
    return (; Mfull, kpb_kfull, kpb_Gfull)
end


function unfold_eig(ieig, f2i)
    nkpoints_full = size(f2i)[1]
    num_bands     = size(ieig[1])[1]
    eig = [zeros(Float64, num_bands) for _ in 1:nkpoints_full ]

    for ik_fbz in 1:nkpoints_full
        ik_ibz = f2i[ik_fbz][1][1]
        eig[ik_fbz] = ieig[ik_ibz]
    end

    return eig
end
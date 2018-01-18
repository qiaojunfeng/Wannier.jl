type WannierParameters
	N1::Int
	N2::Int
	N3::Int
	t1::Array{Float64,1}
	t2::Array{Float64,1}
	t3::Array{Float64,1}
	M::Array{Complex128,4} # nband x nband x Ntot x nntot
	A::Array{Complex128,5} # nband x nwannier x N1 x N2 x N3
	eigs::Array{Float64,2} # Ntot x nband
	neighbors::Array{Int,2} # Ntot x nntot
        kpt_order::Array{Int,1} # permutation such that the order of the kpoints in the file has perm[3] as fastest index
	ijk_to_K :: Array{Int,3}
	K_to_ijk :: Array{Int,2}
	nntot::Int
	nband::Int
        nwannier :: Int
	filename::String
	map::Bool
	logMethod::Bool
        recp_unit_cell :: Array{Float64,2} #3 x 3
        wb :: Float64 #must be inputed manually for now, and be the same for all b vectors
        displ_vecs :: Array{Float64,3} # integer array of displacements between Brillouin zones, 3 x Ntot x nntot
end

# read mmn, win and optionally amn
function read_system(filename, read_amn=true,read_eig=true)
    println("Reading $filename.win")
    win = open("$filename.win")
    nwannier = -1
    nband = -1
    N1,N2,N3 = -1,-1,-1
    unit_cell = zeros(3,3)
    fill!(unit_cell,NaN)
    Ns = (-1,-1,-1)
    kpt_order = [1,2,3]
    while !eof(win)
        line = readline(win)
	line = lowercase(line) # handle case insensitive win files (relic of Fortran)
        line = replace(line, "="," ")
        line = replace(line, ":"," ")
        if(startswith(line,"!"))
            continue
        elseif contains(line,"mp_grid")
            N1,N2,N3 = Base.map(x -> parse(Int,x), split(line)[2:4])
            Ns = [N1,N2,N3]
        elseif contains(line,"num_bands")
            nband = parse(Int, split(line)[2])
        elseif contains(line,"num_wann")
            nwannier = parse(Int, split(line)[2])
        elseif contains(line,"begin unit_cell_cart")
            readline(win) # skip "bohr"
            for i = 1:3
                unit_cell[:,i] = Base.map(x -> parse(Float64,x), split(readline(win))) #weird convention...
            end
        elseif contains(line,"begin kpoints")
            # try to guess the kpoint format
            read_arr() = Base.map(x -> parse(Float64,x), split(readline(win)))
            @assert read_arr() == [0.,0.,0.] #first line
            #second line: detect fastest changing index -> kpt_order[3]
            arr = read_arr()
            nonzero_index = find(x-> x != 0., arr)
            @assert length(nonzero_index) == 1
            @assert isapprox(sum(arr), 1/Ns[nonzero_index[1]],rtol=1e-2)
            kpt_order[3] = nonzero_index[1]

            #read rest of the first part to make sure it's correct
            for i=2:Ns[nonzero_index[1]]-1
                arr = read_arr()
                nonzero_index = find(x-> x != 0., arr)
                @assert length(nonzero_index) == 1
                @assert isapprox(sum(arr), i/Ns[nonzero_index[1]],rtol=1e-2)
            end

            if count(x -> x != 1, Ns) >= 2
                # get the second fastest changing index
                arr = read_arr()
                nonzero_index = find(x-> x != 0., arr)
                @assert length(nonzero_index) == 1
                @assert isapprox(sum(arr), 1/Ns[nonzero_index[1]],rtol=1e-2)
                kpt_order[2] = nonzero_index[1]
            else
                # special 1D case
                kpt_order[2] = kpt_order[3] == 1 ? 2 : 1
            end

            # get last fastest changing index
            kpt_order[1] = 6 - sum(kpt_order[2:3])
        end
    end
    @assert N1 != -1
    @assert N2 != -1
    @assert N3 != -1
    @assert nwannier != -1
    if nband == -1
        nband = nwannier
    end
    @assert all(x->!isnan(x), unit_cell)

    # TODO implement that logic (a point has two neighbors that are the same)
    @assert N1 != 2
    @assert N2 != 2
    @assert N3 != 2

    recp_unit_cell = 2pi * inv(unit_cell)'
    
    println("Reading $filename.mmn")
    mmn = open("$filename.mmn")
    readline(mmn) # skip header
    line = readline(mmn)
    nband2,Ntot,nntot = Base.map(x -> parse(Int64,x),split(line))
    @assert nband2 == nband
    @assert Ntot == N1*N2*N3

    M = zeros(Complex128,nband,nband,Ntot,nntot) # overlap matrix, (K) representation
    neighbors = zeros(Int64,Ntot,nntot) # for each point, list of neighbors, (K) representation
    displ_vecs = zeros(Int64,3,Ntot,nntot)
    while !eof(mmn)
        for nneighbor = 1:nntot
            line = readline(mmn)
            arr = split(line)
            K = parse(Int64, arr[1])
            Kpb = parse(Int64, arr[2])
            displ_vecs[:,K,nneighbor] = [parse(Int64, arr[3]),parse(Int64, arr[4]),parse(Int64, arr[5])]
            neighbors[K,nneighbor,:] = Kpb
            for mn = 0:nband^2-1
                m,n = mod(mn,nband)+1, div(mn,nband)+1
                line = readline(mmn)
                arr = split(line)
                ol = parse(Float64, arr[1]) + im*parse(Float64, arr[2])
                M[m,n,K,nneighbor] = ol
		@assert !isnan(ol)
            end
        end
    end


    t1 = collect(0:N1-1)/N1
    t2 = collect(0:N2-1)/N2
    t3 = collect(0:N3-1)/N3

    # We switch between a big index K=1:N^3 and three indices i,j,k = 1:N using these arrays
    K_to_ijk = zeros(Int64,Ntot,3)
    ijk_to_K = zeros(Int64,N1,N2,N3)
    Ns = [N1,N2,N3]
    for i=0:N1-1
        for j=0:N2-1
            for k=0:N3-1
                ijk = [i,j,k]
                big_offset = Ns[kpt_order[1]]*Ns[kpt_order[2]]
                med_offset = Ns[kpt_order[2]]
                K = ijk[kpt_order[1]]*big_offset + ijk[kpt_order[2]]*med_offset + ijk[kpt_order[3]]+1
                ijk_to_K[i+1,j+1,k+1] = K
                K_to_ijk[K,:] = [i+1 j+1 k+1]
            end
        end
    end

    
    #try to guess wb from (B1) of MV
    sumb = zeros(3,3)
    @assert ijk_to_K[1,1,1] == 1
    for ib = 1:nntot
        neighbor = neighbors[1,ib]
        i_n,j_n,k_n = K_to_ijk[neighbor,:]
	b = recp_unit_cell*([t1[i_n];t2[j_n];t3[k_n]] .+ displ_vecs[:,1,ib] .- [t1[1];t2[1];t3[1]])
        sumb += b*b'
    end
    @assert sumb / norm(sumb) ≈ diagm([Int(N1>1), Int(N2>1), Int(N3>1)])
    wb = 1/norm(sumb)
    println("Finite difference condition satisfied, wb=$wb")

    A = zeros(Complex128,nband,nwannier,N1,N2,N3)
    if read_amn
        println("Reading $filename.amn")
        amn = open("$filename.amn")
        readline(amn)
        arr = split(readline(amn))
        @assert parse(Int64, arr[1]) == nband
        @assert parse(Int64, arr[2]) == N1*N2*N3
        @assert parse(Int64, arr[3]) == nwannier
        lines = readlines(amn)
        for line in lines
            arr = split(line)
            m = parse(Int64, arr[1])
            n = parse(Int64, arr[2])
            kpt = parse(Int64, arr[3])
            Aijkmn = parse(Float64, arr[4]) + im*parse(Float64, arr[5])
            A[m,n,K_to_ijk[kpt,1], K_to_ijk[kpt,2], K_to_ijk[kpt,3]] = Aijkmn
        end

        # normalization
        for i=1:N1,j=1:N2,k=1:N3
            A[:,:,i,j,k] = normalize_matrix(A[:,:,i,j,k])
        end
    else
        for n = 1:nwannier
            A[n,n,:,:,:] = 1
        end
    end

    eigs = zeros(Ntot,nband)
    if read_eig
        println("Reading $filename.eig")
        eigfile = open("$filename.eig")
        while(!eof(eigfile))
            arr = split(readline(eigfile))
            iband = parse(Int,arr[1])
            iK = parse(Int,arr[2])
            val = parse(Float64,arr[3])
            eigs[iK,iband] = val
        end
    end

    map = true
    logMethod = true

    println("$nband bands, $nwannier Wannier, $N1 x $N2 x $N3 grid, $nntot neighbors")

    return WannierParameters(N1,N2,N3,t1,t2,t3,M,A,eigs,neighbors,kpt_order,ijk_to_K,K_to_ijk,nntot,nband,nwannier,filename,map,logMethod,recp_unit_cell,wb,displ_vecs)
end
    

# Computes overlap between two neighboring K points
function overlap(K1,K2,p)
    #fix useful if N3 = 1 (e.g. for 2D models)
    if(K1 == K2)
    	return eye(Complex128,p.nband)
    end
    for nneighbor=1:p.nntot
        if p.neighbors[K1,nneighbor] == K2
            return view(p.M, :, :, K1, nneighbor)
        end
    end
    error("No neighbors found, K1 = $(K1), K2 = $(K2)")
    return eye(Complex128,p.nband)
end
# Computes the overlap between K1 and K2, rotated by A1 and A2 respectively
function overlap(v1,v2,A1,A2,p)
    return (A1')*overlap(v1, v2, p)*A2
end
# Computes the overlap between K1 and K2, using the rotation contained in the array A
function overlap_A(v1,v2,p)
    i1,j1,k1 = v1
    i2,j2,k2 = v2
    K1 = p.ijk_to_K[i1,j1,k1]
    K2 = p.ijk_to_K[i2,j2,k2]
    A1 = p.A[:,:,i1,j1,k1]
    A2 = p.A[:,:,i2,j2,k2]
    return overlap(K1,K2,A1,A2,p)
end

# Power of a unitary (or at least, normal) matrix A
function powm(A,p)
    d,V = eig(A)
    V = normalize_matrix(V)
    accuracy = vecnorm(V*Diagonal(d)*V'-A)
    @assert(accuracy < 1e-10)
    return V*Diagonal(d.^p)*V'
end

# Normalize a matrix A to be unitary. If X is a matrix with orthogonal columns and A a non-singular matrix, then Löwdin-orthogonalizing X*A is equivalent to computing X*normalize_matrix(A)
function normalize_matrix(A)
    U,S,V = svd(A)
    accuracy = vecnorm(U*Diagonal(S)*V'-A)
    @assert(accuracy < 1e-10)
    return U*V'
end

function normalize_matrix_chol(A)
    ovl = A'A
    return A/chol(ovl)
end

# normalize and freeze a block of a matrix
# A = [Uf* Ur*]
# A*A = I, Uf* Uf = I
# UfUf* + UrUr* = I
# From this, obtain Uf*Ur = 0
# Strategy: first orthogonalize Uf, then project Uf out of Ur, then orthogonalize the range of Ur
function normalize_and_freeze(A,frozen,not_frozen)
    # orthogonalize Uf
    Uf = A[frozen,:]'
    U,S,V = svd(Uf)
    Uf = U*V'
    # Uf = normalize_matrix_chol(Uf)

    # project Uf out of Ur
    Ur = A[not_frozen,:]'
    Ur -= Uf*Uf'*Ur

    # # alternative method, maybe more stable but slower
    # ovl = Ur'Ur
    # S, U = eig(Hermitian(ovl))
    # S = real(S)
    # @assert !any(x -> 1e-11 <= x <= 1e-9, S)
    # @assert count(x -> x > 1e-10, S) == size(A,2) - nfrozen
    # Sm12 = map(x-> x < 1e-10 ? 0. : 1/sqrt(x), S)
    # Ur = Ur*(U*diagm(Sm12)*U')

    # renormalize the range of Ur
    U,S,V = svd(Ur)
    eps = 1e-10
    @assert !any(x -> 1e-11 <= x <= 1e-9, S)
    @assert count(x -> x > 1e-10, S) == size(A,2) - count(frozen)
    S[S .> eps] = 1
    S[S .< eps] = 0
    Ur = U*Diagonal(S)*V'

    A[not_frozen,:] = Ur'

    
    B = vcat(Uf',Ur')
    B[frozen,:] = Uf'
    B[not_frozen,:] = Ur'
    @assert B'B ≈ eye(size(A,2)) rtol=1e-12
    @assert B[frozen,:]*B[frozen,:]' ≈ eye(count(frozen)) rtol=1e-12
    @assert vecnorm(Uf'*Ur)<1e-10
    return B
end


# Propagate A0, defined at the first kpt, to the given list of kpts.
# Those must be neighbors, and only the first kpoint is assumed to have been rotated
function propagate(A0, kpts,p)
    N = length(kpts)
    n,m = size(A0)
    @assert(n==m,"Non square matrix given as argument in propagate")
    As = zeros(Complex128,n,n,N)
    As[:,:,1] = A0
    for i0=2:N
        As[:,:,i0] = normalize_matrix(overlap(kpts[i0],kpts[i0-1],p)*As[:,:,i0-1])
    end
    return As
end

# Plot obstructions
function plot_surface_obstructions(p, suffix="")
    if p.N3 != 1
        phases = zeros(p.N1,p.N2,nwannier)
        for i=1:p.N1,j=1:p.N2
            Obs = normalize_matrix(overlap_A([i,j,p.N3],[i,j,1],p)) #rotation at image point
            phases[i,j,:] = sort(imag(log(eig(Obs)[1])))
        end
        figure()
        xx = [p.t1[i] for i=1:p.N1,j=1:p.N2]
        yy = [p.t2[j] for i=1:p.N1,j=1:p.N2]
        for n=1:nwannier
            plot_surface(xx,yy,phases[:,:,n],rstride=1,cstride=1)
        end
        savefig("wannierize$filename.$suffix.pdf")
        close()
    end
end
 
function write_amn(p,A,filename)
    ## Output amn file
    out = open("$filename.amn","w")
    write(out, "Created by wannierize.jl ", string(now()), "\n")
    write(out, "$(p.nband) $(p.N1*p.N2*p.N3) $(p.nwannier)\n")
    for K=1:p.N1*p.N2*p.N3
        for n=1:p.nwannier
            for m = 1:p.nband
                coeff = A[m,n,p.K_to_ijk[K,1], p.K_to_ijk[K,2], p.K_to_ijk[K,3]]
                write(out, "$m $n $K $(real(coeff)) $(imag(coeff))\n")
            end
        end
    end
    close(out)
end

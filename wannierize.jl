module wannierize

const DO_PLOT = false # plot obstruction phases

if(DO_PLOT)
    using PyPlot
end

## Assumptions: the kpoints are contained in a NxNxN cartesian grid, the neighbor list must contain the six cartesian neighbors

## J=nbbands, NxNxN grid, filename.mmn must contain the overlaps
function make_wannier(J,N,nntot,filename)
    t = collect(0:N-1)/N

    # We switch between a big index K=1:N^3 and three indices i,j,k = 1:N using these arrays
    K_to_ijk = zeros(Int64,N^3,3)
    ijk_to_K = zeros(Int64,N,N,N)
    for i=0:N-1
        for j=0:N-1
            for k=0:N-1
                K = i*(N^2)+j*N+k+1
                ijk_to_K[i+1,j+1,k+1] = K
                K_to_ijk[K,:] = [i+1 j+1 k+1]
            end
        end
    end

    A = zeros(Complex128,N,N,N,J,J) # unitary rotation matrix at each k-point, (i,j,k) representation
    M = zeros(Complex128,N^3,nntot,J,J) # overlap matrix, (K) representation
    neighbors = zeros(Int64,N^3,nntot) # for each point, list of neighbors, (K) representation


    # ## In case we want to read an amn file at some point
    # amn = open("$filename.amn")
    # readline(amn)
    # readline(amn)
    # lines = readlines(amn)
    # for line in lines
    #     # print(line)
    #     arr = split(line)
    #     m = parse(Int64, arr[1])
    #     n = parse(Int64, arr[2])
    #     kpt = parse(Int64, arr[3])
    #     Aijkmn = parse(Float64, arr[4]) + im*parse(Float64, arr[5])
    #     A[K_to_ijk[kpt,1], K_to_ijk[kpt,2], K_to_ijk[kpt,3],m,n] = Aijkmn
    # end

    println("Reading $filename.mmn")
    mmn = open("$filename.mmn")
    readline(mmn) # skip headers
    readline(mmn)
    while !eof(mmn)
        for nneighbor = 1:nntot
            line = readline(mmn)
            arr = split(line)
            K = parse(Int64, arr[1])
            Kpb = parse(Int64, arr[2])
            neighbors[K,nneighbor,:] = Kpb
            for mn = 0:J^2-1
                line = readline(mmn)
                arr = split(line)
                ol = parse(Float64, arr[1]) + im*parse(Float64, arr[2])
                m,n = mod(mn,J)+1, div(mn,J)+1
                M[K,nneighbor,m,n] = ol
            end
        end
    end

    fill!(A,NaN) #protection: A must be filled by the algorithm

    # Computes overlap between two neighboring K points
    function overlap(K1,K2)
        for nneighbor=1:nntot
            if neighbors[K1,nneighbor] == K2
                return M[K1,nneighbor,:,:]
            end
        end
        error("No neighbors found")
    end
    # Computes the overlap between K1 and K2, rotated by A1 and A2 respectively
    function overlap(K1,K2,A1,A2)
        return (A1')*overlap(K1, K2)*A2
    end
    # Computes the overlap between K1 and K2, using the rotation contained in the array A
    function overlap_A(K1,K2)
        i1,j1,k1 = K1
        i2,j2,k2 = K2
        return overlap(ijk_to_K[i1,j1,k1],ijk_to_K[i2,j2,k2],A[i1,j1,k1,:,:], A[i2,j2,k2,:,:])
    end

    # Power of a unitary (or at least, normal) matrix A
    function powm(A,p)
        d,V = eig(A)
        return V*diagm(d.^p)*V'
    end

    # Normalize a matrix A to be unitary. If X is a matrix with orthogonal columns and A a non-singular matrix, then Löwdin-orthogonalizing X*A is equivalent to computing X*normalize(A)
    function normalize(A)
        U,S,V = svd(A)
        return U*V'
    end

    # Propagate A0, defined at the first kpt, to the given list of kpts.
    # Those must be neighbors, and only the first kpoint is assumed to have been rotated
    function propagate(A0, kpts)
        N = length(kpts)
        As = zeros(Complex128,N,J,J)
        As[1,:,:] = A0
        for i=2:N
            As[i,:,:] = normalize(overlap(kpts[i],kpts[i-1])*As[i-1,:,:])
            # println("Before/After")
            # println(norm(overlap(kpts[i],kpts[i-1])*As[i-1,:,:] - eye(Complex128,J)))
            # println(norm(As[i,:,:]'*overlap(kpts[i],kpts[i-1])*As[i-1,:,:] - eye(Complex128,J)))
        end
        return As
    end

    function ensure_no_wraparound(A)
        if any(abs(imag(log(eig(A)[1]))) .> 3*pi/4)
            error("Unusually large log! Is your system a topological insulator, or your mesh too coarse? If not, please contact antoine.levitt@inria.fr")
        end
    end

    println("Filling (k,0,0)")
    A[:,1,1,:,:] = propagate(eye(J), [ijk_to_K[i,1,1] for i=1:N])

    # for i=2:N
    #     println(norm(overlap_A([i-1,1,1],[i,1,1]) - eye(Complex128,J)))
    # end

    # compute obstruction matrix
    Obs = normalize(overlap_A([N,1,1],[1,1,1]))
    # and pull it back
    for i=1:N
        A[i,1,1,:,:] = A[i,1,1,:,:]*powm(Obs,t[i])
    end

    # # test
    # for i=2:N
    #     println(norm(overlap_A([i-1,1,1],[i,1,1]) - eye(Complex128,J)))
    # end
    # println(norm(overlap_A([N,1,1],[1,1,1]) - eye(Complex128,J)))


    println("Filling (k1,k2,0)")
    for i=1:N
        A[i,:,1,:,:] = propagate(A[i,1,1,:,:], [ijk_to_K[i,j,1] for j=1:N])
    end

    # corner obstruction
    Obs = normalize(overlap_A([1,N,1],[1,1,1]))
    # pull it back
    for i=1:N
        for j=1:N
            A[i,j,1,:,:] = A[i,j,1,:,:]*powm(Obs,t[j])
        end
    end

    # Pull back the line obstruction
    phases = zeros(N,J)
    for i=1:N
        Obs = normalize(overlap_A([i,N,1],[i,1,1])) #rotation at image point
        ensure_no_wraparound(Obs)
        for j=1:N
            A[i,j,1,:,:] = A[i,j,1,:,:]*powm(Obs,t[j])
        end
        phases[i,:] = imag(log(eig(Obs)[1]))
    end

    if DO_PLOT
        figure()
        plot(phases,"x")
    end


    # phases = zeros(N,J)
    # for i=1:N
    #     Obs = normalize(overlap_A([i,N,1],[i,1,1])) #rotation at image point
    #     phases[i,:] = imag(log(eig(Obs)[1]))
    # end
    # figure()
    # plot(phases,"x")

    # omegaright = zeros(N,N)
    # omegaup = zeros(N,N)
    # for i=1:N,j=1:N
    #     right = i==N ? 1 : i+1
    #     up = j==N ? 1 : j+1
    #     omegaright[i,j] = norm(overlap_A([right,j,1],[i,j,1]) - eye(Complex128,J))
    #     omegaup[i,j] = norm(overlap_A([i,up,1],[i,j,1]) - eye(Complex128,J))
    # end
    # # figure()
    # matshow(omegaright)
    # colorbar()
    # matshow(omegaup)
    # colorbar()

        
    # Plot obstructions
    function plot_surface_obstructions()
        if DO_PLOT
            phases = zeros(N,N,J)
            for i=1:N,j=1:N
                Obs = normalize(overlap_A([i,j,N],[i,j,1])) #rotation at image point
                phases[i,j,:] = sort(imag(log(eig(Obs)[1])))
            end
            figure()
            xx = [t[i] for i=1:N,j=1:N]
            yy = [t[j] for i=1:N,j=1:N]
            for n=1:J
                plot_surface(xx,yy,phases[:,:,n],rstride=1,cstride=1)
            end
        end
    end

    println("Filling (k1,k2,k3)")
    for i=1:N,j=1:N
        A[i,j,:,:,:] = propagate(A[i,j,1,:,:], [ijk_to_K[i,j,k] for k=1:N])
    end

    plot_surface_obstructions()        

    # Fix corner
    Obs = normalize(overlap_A([1,1,N],[1,1,1]))
    for k=1:N
        fixer = powm(Obs,t[k])
        for i=1:N,j=1:N
            A[i,j,k,:,:] = A[i,j,k,:,:]*fixer
        end
    end
    
    # Fix first edge
    for i=1:N
        Obs = normalize(overlap_A([i,1,N], [i,1,1]))
        ensure_no_wraparound(Obs)
        for k=1:N
            fixer = powm(Obs, t[k])
            for j=1:N
                A[i,j,k,:,:] = A[i,j,k,:,:]*fixer
            end
        end
    end
    # Fix second edge
    for j=1:N
        Obs = normalize(overlap_A([1,j,N], [1,j,1]))
        ensure_no_wraparound(Obs)
        for k=1:N
            fixer = powm(Obs, t[k])
            for i=1:N
                A[i,j,k,:,:] = A[i,j,k,:,:]*fixer
            end
        end
    end
    
    plot_surface_obstructions()
    
    # Fix whole surface
    for i=1:N,j=1:N
        Obs = normalize(overlap_A([i,j,N],[i,j,1]))
        ensure_no_wraparound(Obs)
        for k=1:N
            A[i,j,k,:,:] = A[i,j,k,:,:]*powm(Obs,t[k])
        end
    end
        
    plot_surface_obstructions()

    ## Output amn file
    out = open("$filename.amn","w")
    write(out, "Created by wannierize.jl", string(now()),"\n")
    write(out, "$J $(N^3) $J\n")
    for K=1:N^3
        for n=1:J,m=1:J
            coeff = A[K_to_ijk[K,1], K_to_ijk[K,2], K_to_ijk[K,3],m,n]
            write(out, "$m $n $K $(real(coeff)) $(imag(coeff))\n")
        end
    end
    close(out)
end

function read_parameters(filename)
    # Get parameters from mmn file
    mmn = open("$filename.mmn")
    readline(mmn)
    line = readline(mmn)
    nin,N3in,nntotin = split(line)
    nin,N3in,nntotin = parse(Int64,nin), parse(Int64,N3in), parse(Int64,nntotin)
    J = nin
    N = Int64(round(N3in^(1/3)))
    nntot = nntotin
    @assert N^3 == N3in # grid must be N^3
    close(mmn)

    return J,N,nntot
    
    # amn = open("$filename.amn")
    # readline(amn) # skip comment
    # line = readline(amn) # read params
    # nin,N3in,nin2 = split(line)
    # nin,N3in,nin2 = parse(Int64,nin), parse(Int64,N3in), parse(Int64,nin2)
    # @assert nin == J
    # @assert N^3 == N3in
    # @assert nin2 == J
    # close(amn)
end
end

if(length(ARGS) >= 1)
    filename = ARGS[1]
else
    filename = "silicon"
end

J,N,nntot = wannierize.read_parameters(filename)
println("$J bands, $N grid points, $nntot neighbors")
wannierize.make_wannier(J,N,nntot,filename)

using LinearAlgebra
using Dates



# Computes overlap between two neighboring K points
function overlap(K1,K2,p)
    #fix useful if N3 = 1 (e.g. for 2D models)
    if(K1 == K2)
        return Matrix((1.0+0.0im)I,p.nband,p.nband)
    end
    for nneighbor=1:p.nntot
        if p.neighbors[K1,nneighbor] == K2
            return view(p.M, :, :, K1, nneighbor)
        end
    end
    error("No neighbors found, K1 = $(K1), K2 = $(K2)")
    return Matrix((1.0+0.0im)I,p.nband,p.nband)
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
# Actualize the overlap matrix (if A is changed, the overlap should be changed accordingly)
function actualize_M!(p)
    M = p.M
    for K1 = 1:(p.N1*p.N2*p.N3)
	for nneighbor = 1:p.nntot
	    i1,j1,k1 = p.K_to_ijk[K1,:]
	    i2,j2,k2 = p.K_to_ijk[p.neighbors[K1,nneighbor],:]
	    A1 = normalize_matrix(p.A[:,:,i1,j1,k1])
	    A2 = normalize_matrix(p.A[:,:,i2,j2,k2])
	    M[:,:,K1,nneighbor] = A1'*A2
	end
    end
    return M
end
# Power of a unitary (or at least, normal) matrix A
function powm(A,p)
    #Workaround, eigen incompatible with lazy adjoint.
    dV = eigen(Matrix(A))
    d = dV.values
    V = dV.vectors
    V = normalize_matrix(V)
    accuracy = norm(V*Diagonal(d)*V'-A)
    #@assert(accuracy < 1e-10)
    return V*Diagonal(d.^p)*V'
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
    S[S .> eps] .= 1
    S[S .< eps] .= 0
    Ur = U*Diagonal(S)*V'

    A[not_frozen,:] = Ur'


    B = vcat(Uf',Ur')
    B[frozen,:] .= Uf'
    B[not_frozen,:] .= Ur'
    @assert isapprox(B'B, I, rtol=1e-12)
    @assert isapprox(B[frozen,:]*B[frozen,:]', I, rtol=1e-12)
    @assert norm(Uf'*Ur)<1e-10
    return B
end


# Propagate A0, defined at the first kpt, to the given list of kpts.
# Those must be neighbors, and only the first kpoint is assumed to have been rotated
function propagate(A0, kpts,p)
    N = length(kpts)
    n,m = size(A0)
    @assert(n==m,"Non square matrix given as argument in propagate")
    As = zeros(ComplexF64,n,n,N)
    As[:,:,1] = A0
    for i0=2:N
        As[:,:,i0] = normalize_matrix(overlap(kpts[i0],kpts[i0-1],p)*As[:,:,i0-1])
    end
    return As
end

# Plot obstructions
function plot_surface_obstructions(p, suffix="")
    if p.N3 != 1
        phases = zeros(p.N1,p.N2,p.nwannier)
        for i=1:p.N1,j=1:p.N2
            Obs = normalize_matrix(overlap_A([i,j,p.N3],[i,j,1],p)) #rotation at image point
            phases[i,j,:] = sort(imag(log.(eigvals(Obs))))
        end
        figure()
        xx = [p.t1[i] for i=1:p.N1,j=1:p.N2]
        yy = [p.t2[j] for i=1:p.N1,j=1:p.N2]
        for n=1:p.nwannier
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

function interpfft(A,newSize)
    Y = fft(A)
    (N1,N2,N3)  = size(A)
    ##Assuming N1,N2,N3 are even
    Y_padded = zeros(ComplexF64,newSize,newSize,newSize)
    Y_padded[div(newSize-N1,2)+1:div(newSize+N1,2),div(newSize-N2,2)+1:div(newSize+N2,2),div(newSize-N3,2)+1:div(newSize+N3,2)] = fftshift(Y)
    Y_padded = ifftshift(Y_padded)
    Ainterp = (newSize^3/(N1*N2*N3))*real(ifft(Y_padded))
    return Ainterp
end

function interpFourier(A)
    Y = fftshift(fft(A))
    half_size = floor.(size(A)./2)
    function interpA(x)
	poly = 0
	for k = 1:size(A)[3], j = 1:size(A)[2], i = 1:size(A)[1]
	    poly += Y[i,j,k]*exp(im*2*pi*dot([(i-1);(j-1);(k-1)].-half_size,x))
	end
	return real(poly/length(A))
    end
    return interpA
end

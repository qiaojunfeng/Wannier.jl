import Optim
using Random

include("MV.jl")
include("wannierize_utils.jl")

filename = "free"
read_amn = true #read $file.amn as input
read_eig = true #read the eig file (can be set to false when not disentangling)
do_write_amn = true #write $file.optimize.amn at the end
nfrozen = 1 #will freeze if either n <= nfrozen or eigenvalue in frozen window
frozen_window_low = -Inf
frozen_window_high = -Inf
ftol = 1e-20 #tolerance on spread
gtol = 1e-4 #tolerance on gradient
maxiter = 3000 #maximum optimization iterations
m = 100 #history size of BFGS

# expert/experimental features
do_normalize_phase = false # perform a global rotation by a phase factor at the end
do_randomize_gauge = false #randomize initial gauge
cluster_size = 1e-6 #will also freeze additional eigenvalues if the freezing cuts a cluster. Set to 0 to disable
only_r2 = false #only minimize sum_n <r^2>_n, not sum_n <r^2>_n - <r>_n^2

Random.seed!(0)

p = read_system(filename,read_amn,read_eig)

# local number of frozen bands
function local_nfrozen(p, nfrozen, i, j, k, win_low, win_high)
    nf = 0
    @assert win_low == -Inf #for now
    K = p.ijk_to_K[i,j,k]
    nf = count(a -> a < win_high, p.eigs[K,:])
    ret = max(nf,nfrozen)
    @assert ret <= p.nwannier

    # if cluster of eigenvalues, take them all
    @assert issorted(p.eigs[K,:])
    while p.eigs[K,ret+1] < p.eigs[k,ret]+cluster_size
        ret = ret + 1
    end

    return ret
end

function local_frozen_sets(p, nfrozen, i, j, k, win_low, win_high)
    nf = 0
    K = p.ijk_to_K[i,j,k]
    
    l_frozen = (p.eigs[K,:] .>= win_low).&(p.eigs[K,:] .<= win_high)
    l_not_frozen = .!l_frozen

    # if cluster of eigenvalues and nfrozen > 0, take them all
    if nfrozen > 0 
        @assert issorted(p.eigs[K,:])
        while nfrozen < p.nwannier && p.eigs[K,nfrozen+1] < p.eigs[k,nfrozen]+cluster_size
            nfrozen = nfrozen + 1
        end

        l_frozen[1:nfrozen] .= true
        l_not_frozen = .!l_frozen
    end

    @assert count(l_frozen) <= p.nwannier
    return l_frozen, l_not_frozen
end
# There are three formats: A, (X,Y) and XY stored contiguously in memory
# A is the format used in the rest of the code, XY is the format used in the optimizer, (X,Y) is intermediate
# A is nb x nw
# X is nw x nw, Y is nb x nw and has first block equal to [I 0]

function XY_to_A(p,X,Y)
    A = zeros(ComplexF64,p.nband,p.nwannier,p.N1,p.N2,p.N3)
    for i=1:p.N1,j=1:p.N2,k=1:p.N3
        l_frozen, l_not_frozen = local_frozen_sets(p, nfrozen, i, j, k, frozen_window_low, frozen_window_high)
        lnf = count(l_frozen)
        @assert Y[:,:,i,j,k]'Y[:,:,i,j,k] ≈ I
        @assert X[:,:,i,j,k]'X[:,:,i,j,k] ≈ I
        @assert Y[l_frozen,1:lnf,i,j,k] ≈ I
        @assert norm(Y[l_not_frozen,1:lnf,i,j,k]) ≈ 0
        @assert norm(Y[l_frozen,lnf+1:end,i,j,k]) ≈ 0
        A[:,:,i,j,k] = Y[:,:,i,j,k]*X[:,:,i,j,k]
        # @assert normalize_and_freeze(A[:,:,i,j,k],lnf) ≈ A[:,:,i,j,k] rtol=1e-4
    end
    A
end

function A_to_XY(p,A)
    X = zeros(ComplexF64,p.nwannier,p.nwannier,p.N1,p.N2,p.N3)
    Y = zeros(ComplexF64,p.nband,p.nwannier,p.N1,p.N2,p.N3)
    for i=1:p.N1,j=1:p.N2,k=1:p.N3
        l_frozen, l_not_frozen = local_frozen_sets(p, nfrozen, i, j, k, frozen_window_low, frozen_window_high)
        lnf = count(l_frozen)
        Afrozen = normalize_and_freeze(A[:,:,i,j,k], l_frozen, l_not_frozen)
        Af = Afrozen[l_frozen,:]
        Ar = Afrozen[l_not_frozen,:]

        #determine Y
        if lnf != p.nwannier
            proj = Ar*Ar'
            proj = Hermitian((proj+proj')/2)
            D,V = eigen(proj) #sorted by increasing eigenvalue
        end
        Y[l_frozen,1:lnf,i,j,k] = Matrix(I,lnf,lnf)
        if lnf != p.nwannier
            Y[l_not_frozen,lnf+1:end,i,j,k] = V[:,end-p.nwannier+lnf+1:end]
        end
        
        #determine X
        Xleft, S, Xright = svd(Y[:,:,i,j,k]'*Afrozen)
        X[:,:,i,j,k] = Xleft*Xright'

        @assert Y[:,:,i,j,k]'Y[:,:,i,j,k] ≈ I
        @assert X[:,:,i,j,k]'X[:,:,i,j,k] ≈ I
        @assert Y[l_frozen,1:lnf,i,j,k] ≈ I
        @assert norm(Y[l_not_frozen,1:lnf,i,j,k]) ≈ 0
        @assert norm(Y[l_frozen,lnf+1:end,i,j,k]) ≈ 0
        @assert Y[:,:,i,j,k]*X[:,:,i,j,k] ≈ Afrozen
    end
    X,Y
end

function XY_to_XY(p,XY) #XY to (X,Y)
    X = zeros(ComplexF64,p.nwannier,p.nwannier,p.N1,p.N2,p.N3)
    Y = zeros(ComplexF64,p.nband,p.nwannier,p.N1,p.N2,p.N3)
    for i=1:p.N1,j=1:p.N2,k=1:p.N3
        XYk = XY[:,i,j,k]
        X[:,:,i,j,k] = reshape(XYk[1:p.nwannier*p.nwannier], (p.nwannier, p.nwannier))
        Y[:,:,i,j,k] = reshape(XYk[p.nwannier*p.nwannier+1:end], (p.nband, p.nwannier))
    end
    X,Y
end
function obj(p,X,Y)
    A = XY_to_A(p,X,Y)

    res = omega(p,A,true,only_r2)
    func = res.Ωtot
    grad = res.gradient

    gradX = zero(X)
    gradY = zero(Y)
    for i=1:p.N1,j=1:p.N2,k=1:p.N3
        l_frozen, l_not_frozen = local_frozen_sets(p, nfrozen, i, j, k, frozen_window_low, frozen_window_high)
        lnf = count(l_frozen)
        gradX[:,:,i,j,k] = Y[:,:,i,j,k]'*grad[:,:,i,j,k]
        gradY[:,:,i,j,k] = grad[:,:,i,j,k]*X[:,:,i,j,k]'

        gradY[l_frozen,:,i,j,k] .= 0
        gradY[:,1:lnf,i,j,k] .= 0

        # to compute the projected gradient: redundant, taken care of by the optimizer
        # function proj_stiefel(G,X)
        #     G .- X*((X'G .+ G'X)./2)
        # end
        # gradX[:,:,i,j,k] = proj_stiefel(gradX[:,:,i,j,k],X[:,:,i,j,k])
        # gradY[:,:,i,j,k] = proj_stiefel(gradY[:,:,i,j,k],Y[:,:,i,j,k])
    end
    return func, gradX,gradY, res
end

function minimize(p,A)
    # initial X,Y
    X0,Y0 = A_to_XY(p,A)
    
    if do_randomize_gauge
        if read_amn
            error("don't set do_randomize_gauge and read_amn")
        end
        X0 = zeros(ComplexF64,p.nwannier,p.nwannier,p.N1,p.N2,p.N3)
        Y0 = zeros(ComplexF64,p.nband,p.nwannier,p.N1,p.N2,p.N3)
        for i=1:p.N1,j=1:p.N2,k=1:p.N3
            l_frozen, l_not_frozen = local_frozen_sets(p, nfrozen, i, j, k, frozen_window_low, frozen_window_high)
            lnf = count(l_frozen)
            X0[:,:,i,j,k] = normalize_matrix(randn(p.nwannier,p.nwannier) + im*randn(p.nwannier,p.nwannier))
            Y0[l_frozen,1:lnf,i,j,k] = eye(lnf)
            Y0[l_not_frozen, lnf+1:p.nwannier,i,j,k] = normalize_matrix(randn(p.nband-lnf,p.nwannier-lnf) + im*randn(p.nband-lnf,p.nwannier-lnf))
        end
    end

    M = p.nwannier*p.nwannier+p.nband*p.nwannier
    XY0 = zeros(ComplexF64, M, p.N1, p.N2, p.N3)
    for i=1:p.N1,j=1:p.N2,k=1:p.N3
        XY0[:,i,j,k] = vcat(vec(X0[:,:,i,j,k]),vec(Y0[:,:,i,j,k]))
    end

    # We have three formats:
    # (X,Y): Ntot x nw x nw, Ntot x nb x nw
    # A: ntot x nb x nw
    # XY: (nw*nw + nb*nw) x Ntot
    function fg!(G, XY)
        @assert size(G) == size(XY)
        X,Y = XY_to_XY(p,XY)

        f, gradX, gradY, res = obj(p,X,Y)

        for i=1:p.N1,j=1:p.N2,k=1:p.N3
            G[:,i,j,k] = vcat(vec(gradX[:,:,i,j,k]),vec(gradY[:,:,i,j,k]))
        end
        return f
    end

    function f(XY)
        return fg!(similar(XY),XY)
    end
    function g!(g, XY)
        fg!(g,XY)
        return g
    end

    # need QR orthogonalization rather than SVD to preserve the sparsity structure of Y
    XYkManif = Optim.ProductManifold(Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (p.nwannier, p.nwannier), (p.nband, p.nwannier))
    XYManif = Optim.PowerManifold(XYkManif, (M,), (p.N1,p.N2,p.N3))

    # stepsize_mult = 1
    # step = 0.5/(4*8*p.wb)*(p.N1*p.N2*p.N3)*stepsize_mult
    # ls = LineSearches.Static(step)
    ls = Optim.HagerZhang()
    # ls = LineSearches.BackTracking()
    # meth = Optim.GradientDescent
    # meth = Optim.ConjugateGradient
    meth = Optim.LBFGS
    res = Optim.optimize(f,g!,XY0, meth(manifold=XYManif, linesearch=ls,m=m), Optim.Options(show_trace=true,iterations=maxiter,f_tol=ftol, g_tol=gtol, allow_f_increases=true))
    display(res)
    XYmin = Optim.minimizer(res)
    
    Xmin,Ymin = XY_to_XY(p, XYmin)
    Amin = XY_to_A(p,Xmin,Ymin)
end



if read_amn
    A0 = copy(p.A)
else
    A0 = randn(size(p.A)) + im*randn(size(p.A))
end

for i=1:p.N1,j=1:p.N2,k=1:p.N3
    l_frozen, l_not_frozen = local_frozen_sets(p, nfrozen, i, j, k, frozen_window_low, frozen_window_high)
    A0[:,:,i,j,k] = normalize_and_freeze(A0[:,:,i,j,k],l_frozen,l_not_frozen)
end

A = minimize(p,A0)

# fix global phase
if do_normalize_phase
    for i=1:nwannier
        imax = indmax(abs.(A[:,i,1,1,1]))
        @assert abs(A[imax,i,1,1,1]) > 1e-2
        A[:,i,:,:,:] *= conj(A[imax,i,1,1,1]/abs(A[imax,i,1,1,1]))
    end
end

if do_write_amn
    write_amn(p,A,"$(p.filename).optimized")
end

module Disentangle

import LinearAlgebra as LA
import Optim
# using ..Parameters: InputParams
using ..Spreads: omega

function get_frozen_bands_k(params, ik)# , num_frozen, dis_froz_min, dis_froz_max)
    # number of frozen bands
    # FIX: test only using num_frozen 0, w/o window,projectability
    num_frozen = 0
    # select frozen bands based on energy window
    freeze_energy_window = true
    # dis_froz_min < energy < dis_froz_max bands are frozen
    dis_froz_min = -Inf
    dis_froz_max = 12 #6.5
    # select frozen bands based on projectability
    freeze_projectability = false
    # < threshold bands are excluded
    dis_proj_min = 1 / 100
    # >= threshold bands are frozen
    dis_proj_max = 98 / 100
    # will also freeze additional eigenvalues if the freezing cuts a cluster. Set to 0 to disable
    cluster_size = 1e-6
    
    l_frozen = BitVector()

    if freeze_energy_window
        l_frozen = (params.eig[:, ik] .>= dis_froz_min) .& (params.eig[:, ik] .<= dis_froz_max)
    end

    # if cluster of eigenvalues and nfrozen > 0, take them all
    if num_frozen > 0
        while num_frozen < params.num_wann
            if params.eig[num_frozen + 1, ik] < params.eig[num_frozen, ik] + cluster_size
                num_frozen = num_frozen + 1
            else
                break
            end
        end

        l_frozen[1:num_frozen] .= true
    end

    if freeze_projectability
        # size = nbands * nwannier
        amn_k = params.amn[:,:,ik]
        proj = dropdims(real(sum(amn_k .* conj(amn_k), dims=2)), dims=2)
        @debug "projectability @ $ik" proj

        if any(.!l_frozen[proj .>= dis_proj_max])
            @debug "l_frozen before dis_proj" l_frozen
            l_frozen[proj .>= dis_proj_max] .= true
            @debug "l_frozen after dis_proj" l_frozen
            @debug "proj .>= dis_proj_max" proj .>= dis_proj_max
        end
    end

    @assert count(l_frozen) <= params.num_wann
    return l_frozen
end

"""
normalize and freeze a block of a matrix

Block form:  A = vcat(Uf, Ur)
Semiunitary: A' * A = I
             Uf' * Uf + Ur' * Ur = I
Frozen:      Uf * Uf' = I
Also:        Uf * Ur' = 0
Strategy: first orthogonalize Uf, then project Uf out of Ur, then orthogonalize the range of Ur
"""
function orthonormalize_and_freeze(A::Matrix{ComplexF64}, frozen::BitVector, not_frozen::BitVector)
    # orthogonalize Uf
    Uf = A[frozen,:]
    U, S, V = LA.svd(Uf)
    Uf = U * V'
    # Uf = normalize_matrix_chol(Uf)

    # project Uf out of Ur
    Ur = A[not_frozen,:]
    Ur -= Ur * Uf' * Uf

    # # alternative method, maybe more stable but slower
    # ovl = Ur'Ur
    # S, U = eig(Hermitian(ovl))
    # S = real(S)
    # @assert !any(x -> 1e-11 <= x <= 1e-9, S)
    # @assert count(x -> x > 1e-10, S) == size(A,2) - nfrozen
    # Sm12 = map(x-> x < 1e-10 ? 0. : 1/sqrt(x), S)
    # Ur = Ur*(U*diagm(Sm12)*U')

    # renormalize the range of Ur
    U, S, V = LA.svd(Ur)
    eps = 1e-10
    @assert !any(x -> 1e-11 <= x <= 1e-9, S)
    @assert count(x -> x > 1e-10, S) == size(A, 2) - count(frozen)
    S[S .> eps] .= 1
    S[S .< eps] .= 0
    Ur = U * LA.Diagonal(S) * V'

    A[not_frozen,:] = Ur

    B = vcat(Uf, Ur)
    B[frozen,:] .= Uf
    B[not_frozen,:] .= Ur
    @assert isapprox(B'B, LA.I, rtol=1e-12)
    @assert isapprox(B[frozen,:] * B[frozen,:]', LA.I, rtol=1e-12)
    @assert LA.norm(Uf * Ur') < 1e-10
    
    return B
end

# There are three formats: A, (X,Y) and XY stored contiguously in memory
# A is the format used in the rest of the code, XY is the format used in the optimizer, (X,Y) is intermediate
# A is num_bands * num_wann
# X is num_wann * num_wann, Y is num_bands * num_wann and has first block equal to [I 0]
function XY_to_A(params, X, Y)
    A = zeros(ComplexF64, params.num_bands, params.num_wann, params.num_kpts)
    for ik = 1:params.num_kpts
        # check
        l_frozen = get_frozen_bands_k(params, ik)
        l_non_frozen = .!l_frozen
        num_frozen = count(l_frozen)
        @assert Y[:,:,ik]' * Y[:,:,ik] ≈ LA.I
        @assert X[:,:,ik]' * X[:,:,ik] ≈ LA.I
        @assert Y[l_frozen,1:num_frozen,ik] ≈ LA.I
        @assert LA.norm(Y[l_non_frozen,1:num_frozen,ik]) ≈ 0
        @assert LA.norm(Y[l_frozen,num_frozen + 1:end,ik]) ≈ 0

        A[:,:,ik] = Y[:,:,ik] * X[:,:,ik]
        # @assert normalize_and_freeze(A[:,:,i,j,k],lnf) ≈ A[:,:,i,j,k] rtol=1e-4
    end
    return A
end

function A_to_XY(params, A)
    X = zeros(ComplexF64, params.num_wann, params.num_wann, params.num_kpts)
    Y = zeros(ComplexF64, params.num_bands, params.num_wann, params.num_kpts)
    for ik = 1:params.num_kpts
        l_frozen = get_frozen_bands_k(params, ik)
        l_non_frozen = .!l_frozen
        num_frozen = count(l_frozen)
        Afrozen = orthonormalize_and_freeze(A[:,:,ik], l_frozen, l_non_frozen)
        Af = Afrozen[l_frozen,:]
        Ar = Afrozen[l_non_frozen,:]

        # determine Y
        if num_frozen != params.num_wann
            proj = Ar * Ar'
            proj = LA.Hermitian((proj + proj') / 2)
            D, V = LA.eigen(proj) # sorted by increasing eigenvalue
        end
        Y[l_frozen,1:num_frozen,ik] = Matrix(LA.I, num_frozen, num_frozen)
        if num_frozen != params.num_wann
            Y[l_non_frozen,num_frozen + 1:end,ik] = V[:,end - params.num_wann + num_frozen + 1:end]
        end
        
        # determine X
        Xleft, S, Xright = LA.svd(Y[:,:,ik]' * Afrozen)
        X[:,:,ik] = Xleft * Xright'

        @assert Y[:,:,ik]' * Y[:,:,ik] ≈ LA.I
        @assert X[:,:,ik]' * X[:,:,ik] ≈ LA.I
        @assert Y[l_frozen,1:num_frozen,ik] ≈ LA.I
        @assert LA.norm(Y[l_non_frozen,1:num_frozen,ik]) ≈ 0
        @assert LA.norm(Y[l_frozen,num_frozen + 1:end,ik]) ≈ 0
        @assert Y[:,:,ik] * X[:,:,ik] ≈ Afrozen
    end
    return X, Y
end

function XY_to_XY(params, XY) # XY to (X,Y)
    X = zeros(ComplexF64, params.num_wann, params.num_wann, params.num_kpts)
    Y = zeros(ComplexF64, params.num_bands, params.num_wann, params.num_kpts)
    for ik = 1:params.num_kpts
        XYk = XY[:,ik]
        X[:,:,ik] = reshape(XYk[1:params.num_wann^2], (params.num_wann, params.num_wann))
        Y[:,:,ik] = reshape(XYk[params.num_wann^2 + 1:end], (params.num_bands, params.num_wann))
    end
    return X, Y
end

function obj(params, X, Y)
    A = XY_to_A(params, X, Y)

    #
    only_r2 = false
    res = omega(params, A, true, only_r2)
    func = res.Ωtot
    grad = res.gradient

    gradX = zero(X)
    gradY = zero(Y)
    for ik = 1:params.num_wann
        l_frozen = get_frozen_bands_k(params, ik)
        l_non_frozen = .!l_frozen
        num_frozen = count(l_frozen)
        gradX[:,:,ik] = Y[:,:,ik]' * grad[:,:,ik]
        gradY[:,:,ik] = grad[:,:,ik] * X[:,:,ik]'

        gradY[l_frozen,:,ik] .= 0
        gradY[:,1:num_frozen,ik] .= 0

        # to compute the projected gradient: redundant, taken care of by the optimizer
        # function proj_stiefel(G,X)
        #     G .- X*((X'G .+ G'X)./2)
        # end
        # gradX[:,:,i,j,k] = proj_stiefel(gradX[:,:,i,j,k],X[:,:,i,j,k])
        # gradY[:,:,i,j,k] = proj_stiefel(gradY[:,:,i,j,k],Y[:,:,i,j,k])
    end
    return func, gradX, gradY, res
end

function minimize(params, A)
    # tolerance on spread
    ftol = 1e-20
    # tolerance on gradient
    gtol = 1e-4
    # maximum optimization iterations
    maxiter = 200 # 3000
    # history size of BFGS
    m = 100

    # initial X,Y
    X0, Y0 = A_to_XY(params, A)
    
    do_randomize_gauge = false
    if do_randomize_gauge
        X0 = zeros(ComplexF64, params.num_wann, params.num_wann, params.num_kpts)
        Y0 = zeros(ComplexF64, params.num_bands, params.num_wann, params.num_kpts)
        for ik = 1:params.num_kpts
            l_frozen = get_frozen_bands_k(params, ik)
            l_non_frozen = .!l_frozen
            num_frozen = count(l_frozen)
            X0[:,:,ik] = normalize_matrix(
                randn(params.num_wann, params.num_wann) + im * randn(params.num_wann, params.num_wann))
            Y0[l_frozen,1:num_frozen,ik] = eye(num_frozen)
            Y0[l_non_frozen, num_frozen + 1:params.num_wann,ik] = normalize_matrix(
                randn(params.num_bands - num_frozen, params.num_wann - num_frozen) +
                im * randn(params.num_bands - num_frozen, params.num_wann - num_frozen))
        end
    end

    M = params.num_wann^2 + params.num_bands * params.num_wann
    XY0 = zeros(ComplexF64, M, params.num_kpts)
    for ik = 1:params.num_kpts
        XY0[:,ik] = vcat(vec(X0[:,:,ik]), vec(Y0[:,:,ik]))
    end

    # We have three formats:
    # (X,Y): num_wann * num_wann * num_kpts, num_bands * num_wann * num_kpts
    # A: num_bands * num_wann * num_kpts
    # XY: (num_wann * num_wann + num_bands * num_wann) * num_kpts
    function fg!(G, XY)
        @assert size(G) == size(XY)
        X, Y = XY_to_XY(params, XY)

        f, gradX, gradY, res = obj(params, X, Y)

        for ik = 1:params.num_kpts
            G[:,ik] = vcat(vec(gradX[:,:,ik]), vec(gradY[:,:,ik]))
        end
        return f
    end

    function f(XY)
        return fg!(similar(XY), XY)
    end
    function g!(g, XY)
        fg!(g, XY)
        return g
    end

    tmp = omega(params, A, false, false)
    @info "Initial centers & spreads" tmp.centers tmp.spreads' sum(tmp.spreads)

    # need QR orthogonalization rather than SVD to preserve the sparsity structure of Y
    XYkManif = Optim.ProductManifold(Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), 
        (params.num_wann, params.num_wann), (params.num_bands, params.num_wann))
    XYManif = Optim.PowerManifold(XYkManif, (M,), (params.num_kpts,))

    # stepsize_mult = 1
    # step = 0.5/(4*8*p.wb)*(p.N1*p.N2*p.N3)*stepsize_mult
    # ls = LineSearches.Static(step)
    ls = Optim.HagerZhang()
    # ls = LineSearches.BackTracking()
    # meth = Optim.GradientDescent
    # meth = Optim.ConjugateGradient
    meth = Optim.LBFGS
    res = Optim.optimize(f, g!, XY0, meth(manifold=XYManif, linesearch=ls, m=m), 
        Optim.Options(show_trace=true, iterations=maxiter, f_tol=ftol, g_tol=gtol, allow_f_increases=true))
    display(res)
    XYmin = Optim.minimizer(res)
    
    Xmin, Ymin = XY_to_XY(params, XYmin)
    Amin = XY_to_A(params, Xmin, Ymin)

    tmp = omega(params, Amin, false, false)
    @info "Final centers & spreads" tmp.centers tmp.spreads' sum(tmp.spreads)

    return Amin
end

end

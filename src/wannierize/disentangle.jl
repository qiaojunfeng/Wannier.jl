import LinearAlgebra as LA
import Optim
include("../param.jl")
include("../spread.jl")


function set_frozen_bands!(data::CoreData, params::InputParams)

    data.frozen = falses(data.num_bands, data.num_kpts)

    for ik = 1:data.num_kpts
        frozen_k = falses(data.num_bands)
        num_frozen_k = 0

        if params.dis_froz_win
            frozen_k = (data.eig[:, ik] .>= params.dis_froz_win_min) .&
                       (data.eig[:, ik] .<= params.dis_froz_win_max)
            # @debug "ik = $ik, l_frozen_k = $(l_frozen_k)"
            num_frozen_k = count(frozen_k)
        end

        if params.dis_froz_proj
            # size = nbands * nwannier
            amn_k = data.amn[:, :, ik]
            proj = dropdims(real(sum(amn_k .* conj(amn_k), dims=2)), dims=2)
            # @debug "projectability @ $ik" proj

            if any(.!frozen_k[proj.>=params.dis_froz_proj_max])
                # @debug "l_frozen before dis_proj" l_frozen'
                frozen_k[proj.>=params.dis_froz_proj_max] .= true
                # @debug "l_frozen after dis_proj" l_frozen'
                # @debug "proj .>= dis_proj_max" (proj .>= dis_proj_max)'
                num_frozen_k = count(frozen_k)
            end
        end

        # if cluster of eigenvalues and nfrozen > 0, take them all
        if params.dis_froz_degen && params.dis_froz_degen_thres > 0 && num_frozen_k > 0
            while num_frozen_k < data.num_wann
                if data.eig[num_frozen_k+1, ik] < data.eig[num_frozen_k, ik] + params.dis_froz_degen_thres
                    num_frozen_k += 1
                    frozen_k[num_frozen_k] .= true
                else
                    break
                end
            end
        end

        if params.dis_froz_num > 0
            frozen_k[1:params.dis_froz_num] .= true
        end

        @assert count(frozen_k) <= data.num_wann
        # @debug "num_frozen_k = $num_frozen_k, ik = $ik" frozen_k'
        data.frozen[:, ik] = frozen_k
    end
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
function orthonormalize_and_freeze(A::Matrix{ComplexF64}, frozen::BitVector, non_frozen::BitVector)
    # Make sure Uf can fully represent frozen bands.
    # Uf = <ψ|g>, where |ψ> is Bloch wfcs, |g> is the guiding functions (GF).
    # We do a Lowdin orthonormalization on Uf' = <g|ψ> so (Uf')' * Uf' = I, i.e. orthonormal,
    # this also means Uf * Uf' = I, i.e. <ψ|g><g|ψ> = I
    # =>  |g> can fully represent frozen bands without loss.
    Uf = A[frozen, :]
    U, S, V = LA.svd(Uf')
    Uf = V * U'
    # Uf = normalize_matrix_chol(Uf')'

    # Remove Uf out of Ur, i.e. do not destroy frozen space
    # The projector (i.e. density matrix) of the frozen states represented on the GF basis is 
    #     |ψf><ψf| = |g><g|ψf><ψf|g><g| = |g> Uf' * Uf <g|
    # The projector of the non-frozen states on GF basis is
    #     |ψr><ψr| = |g><g|ψr><ψr|g><g| = |g> Ur' * Ur <g|
    # The space |ψr><ψr| should not destroy the frozen space |ψf><ψf|.
    # To achieve this, we remove |ψf><ψf| components out of |ψr><ψr|
    #     |ψr> -= |ψf><ψf|ψr>
    # =>  |g><g|ψr> -= |g><g|ψf><ψf|g><g|ψr>
    # =>  Ur' -= Uf' * Uf * Ur'
    Ur = A[non_frozen, :]
    Ur -= Ur * Uf' * Uf

    # alternative method, maybe more stable but slower
    # ovl = Ur'Ur
    # S, U = eig(Hermitian(ovl))
    # S = real(S)
    # @assert !any(x -> 1e-11 <= x <= 1e-9, S)
    # @assert count(x -> x > 1e-10, S) == size(A,2) - nfrozen
    # Sm12 = map(x-> x < 1e-10 ? 0. : 1/sqrt(x), S)
    # Ur = Ur*(U*diagm(Sm12)*U')

    # Renormalize the range of Ur
    # The remaining Wannier function (WF) |wr> = |ψ> Ur,
    # after removal of Uf, we need to renormalize so the |wr> are orthonormal.
    # I = <wr|wr> = Ur' <ψ|ψ> Ur  =>  Ur' * Ur = I
    # Use Lowdin normalization but needs to limit the number of independent vectors.
    U, S, V = LA.svd(Ur)
    eps = 1e-10
    @assert !any(x -> 1e-11 <= x <= 1e-9, S)
    @assert count(x -> x > eps, S) == size(A, 2) - count(frozen)
    S[S.>eps] .= 1
    S[S.<eps] .= 0
    Ur = U * LA.Diagonal(S) * V'

    # FIX: remove this?
    # A[non_frozen,:] = Ur

    B = vcat(Uf, Ur)
    B[frozen, :] .= Uf
    B[non_frozen, :] .= Ur

    # Semiunitary
    @assert isapprox(B' * B, LA.I, rtol=1e-12)
    # Frozen
    @assert isapprox(B[frozen, :] * B[frozen, :]', LA.I, rtol=1e-12)
    # Independent
    @assert LA.norm(Uf * Ur') < 1e-10

    return B
end

function max_projectability(A::Matrix{ComplexF64})
    proj = A * A'
    proj_ortho = LA.I - proj
    U, S, V = LA.svd(proj_ortho)
    n = abs(size(A, 1) - size(A, 2))
    @assert count(S .> 1e-5) >= n
    R = hcat(A, V[:, 1:n])
    @assert R * R' ≈ LA.I
    return R' * A, R

    U, S, V = LA.svd(A)
    return V' * A, V
    proj = A * A'
    # A * A' = V D V'  =>  (V'A) (A'V) = D
    D, V = LA.eigen(proj)
    # sort eigenvalues in descending order
    V = V[:, sortperm(D, rev=true)]
    A_new = V' * A
    @debug "projectability" real(LA.diag(proj)') real(LA.diag(A_new * A_new')')
    return A_new, V
end

function maxproj_froz(A::Matrix{ComplexF64}, froz::BitVector)
    @assert length(froz) == size(A, 1)
    m, n = size(A)
    D = zeros(ComplexF64, n + count(froz), m)
    D[1:n, :] = transpose(A)
    c = n + 1
    for i = 1:length(froz)
        if froz[i]
            D[c, i] = 1
            c += 1
        end
    end
    U, S, V = LA.svd(D)

    nbasis = count(S .> 1e-5)
    V_new = V[:, 1:nbasis]
    proj = V_new * V_new'
    proj_ortho = LA.I - proj
    U, S, V = LA.svd(proj_ortho)
    R = hcat(V_new, V[:, 1:m-nbasis])
    @assert R * R' ≈ LA.I

    return R' * A, R
end

# There are three formats: A, (X,Y) and XY stored contiguously in memory
# A is the format used in the rest of the code, XY is the format used in the optimizer, (X,Y) is intermediate
# A is num_bands * num_wann
# X is num_wann * num_wann, Y is num_bands * num_wann and has first block equal to [I 0]
function XY_to_A(data, X, Y)
    A = zeros(ComplexF64, data.num_bands, data.num_wann, data.num_kpts)
    for ik = 1:data.num_kpts
        # check
        l_frozen = data.frozen[:, ik]
        l_non_frozen = .!l_frozen
        num_frozen = count(l_frozen)
        @assert Y[:, :, ik]' * Y[:, :, ik] ≈ LA.I
        @assert X[:, :, ik]' * X[:, :, ik] ≈ LA.I
        @assert Y[l_frozen, 1:num_frozen, ik] ≈ LA.I
        @assert LA.norm(Y[l_non_frozen, 1:num_frozen, ik]) ≈ 0
        @assert LA.norm(Y[l_frozen, num_frozen+1:end, ik]) ≈ 0

        A[:, :, ik] = Y[:, :, ik] * X[:, :, ik]
        # @assert normalize_and_freeze(A[:,:,i,j,k],lnf) ≈ A[:,:,i,j,k] rtol=1e-4
    end
    return A
end

function A_to_XY(data::CoreData, A)
    X = zeros(ComplexF64, data.num_wann, data.num_wann, data.num_kpts)
    Y = zeros(ComplexF64, data.num_bands, data.num_wann, data.num_kpts)
    for ik = 1:data.num_kpts
        l_frozen = data.frozen[:, ik]
        l_non_frozen = .!l_frozen
        num_frozen = count(l_frozen)
        Afrozen = orthonormalize_and_freeze(A[:, :, ik], l_frozen, l_non_frozen)
        Af = Afrozen[l_frozen, :]
        Ar = Afrozen[l_non_frozen, :]

        # determine Y
        if num_frozen != data.num_wann
            proj = Ar * Ar'
            proj = LA.Hermitian((proj + proj') / 2)
            D, V = LA.eigen(proj) # sorted by increasing eigenvalue
        end
        Y[l_frozen, 1:num_frozen, ik] = Matrix(LA.I, num_frozen, num_frozen)
        if num_frozen != data.num_wann
            Y[l_non_frozen, num_frozen+1:end, ik] = V[:, end-data.num_wann+num_frozen+1:end]
        end

        # determine X
        Xleft, S, Xright = LA.svd(Y[:, :, ik]' * Afrozen)
        X[:, :, ik] = Xleft * Xright'

        @assert Y[:, :, ik]' * Y[:, :, ik] ≈ LA.I
        @assert X[:, :, ik]' * X[:, :, ik] ≈ LA.I
        @assert Y[l_frozen, 1:num_frozen, ik] ≈ LA.I
        @assert LA.norm(Y[l_non_frozen, 1:num_frozen, ik]) ≈ 0
        @assert LA.norm(Y[l_frozen, num_frozen+1:end, ik]) ≈ 0
        @assert Y[:, :, ik] * X[:, :, ik] ≈ Afrozen
    end
    return X, Y
end

function XY_to_XY(data, XY) # XY to (X,Y)
    X = zeros(ComplexF64, data.num_wann, data.num_wann, data.num_kpts)
    Y = zeros(ComplexF64, data.num_bands, data.num_wann, data.num_kpts)
    for ik = 1:data.num_kpts
        XYk = XY[:, ik]
        X[:, :, ik] = reshape(XYk[1:data.num_wann^2], (data.num_wann, data.num_wann))
        Y[:, :, ik] = reshape(XYk[data.num_wann^2+1:end], (data.num_bands, data.num_wann))
    end
    return X, Y
end

function obj(data::CoreData, params::InputParams, X, Y)
    A = XY_to_A(data, X, Y)

    #
    only_r2 = false
    res = omega(data, params, A, true, only_r2)
    func = res.Ωtot
    grad = res.gradient

    gradX = zero(X)
    gradY = zero(Y)
    for ik = 1:data.num_kpts
        l_frozen = data.frozen[:, ik]
        l_non_frozen = .!l_frozen
        num_frozen = count(l_frozen)
        gradX[:, :, ik] = Y[:, :, ik]' * grad[:, :, ik]
        gradY[:, :, ik] = grad[:, :, ik] * X[:, :, ik]'

        gradY[l_frozen, :, ik] .= 0
        gradY[:, 1:num_frozen, ik] .= 0

        # to compute the projected gradient: redundant, taken care of by the optimizer
        # function proj_stiefel(G,X)
        #     G .- X*((X'G .+ G'X)./2)
        # end
        # gradX[:,:,i,j,k] = proj_stiefel(gradX[:,:,i,j,k],X[:,:,i,j,k])
        # gradY[:,:,i,j,k] = proj_stiefel(gradY[:,:,i,j,k],Y[:,:,i,j,k])
    end
    return func, gradX, gradY, res
end

function minimize(data::CoreData, params::InputParams, A)
    # initial X,Y
    X0, Y0 = A_to_XY(data, A)

    do_randomize_gauge = false
    if do_randomize_gauge
        X0 = zeros(ComplexF64, data.num_wann, data.num_wann, data.num_kpts)
        Y0 = zeros(ComplexF64, data.num_bands, data.num_wann, data.num_kpts)
        for ik = 1:data.num_kpts
            l_frozen = data.frozen[:, ik]
            l_non_frozen = .!l_frozen
            num_frozen = count(l_frozen)
            X0[:, :, ik] = normalize_matrix(
                randn(data.num_wann, data.num_wann) + im * randn(data.num_wann, data.num_wann))
            Y0[l_frozen, 1:num_frozen, ik] = eye(num_frozen)
            Y0[l_non_frozen, num_frozen+1:data.num_wann, ik] = normalize_matrix(
                randn(data.num_bands - num_frozen, data.num_wann - num_frozen) +
                im * randn(data.num_bands - num_frozen, data.num_wann - num_frozen))
        end
    end

    M = data.num_wann^2 + data.num_bands * data.num_wann
    XY0 = zeros(ComplexF64, M, data.num_kpts)
    for ik = 1:data.num_kpts
        XY0[:, ik] = vcat(vec(X0[:, :, ik]), vec(Y0[:, :, ik]))
    end

    # We have three formats:
    # (X,Y): num_wann * num_wann * num_kpts, num_bands * num_wann * num_kpts
    # A: num_bands * num_wann * num_kpts
    # XY: (num_wann * num_wann + num_bands * num_wann) * num_kpts
    function fg!(G, XY)
        @assert size(G) == size(XY)
        X, Y = XY_to_XY(data, XY)

        f, gradX, gradY, res = obj(data, params, X, Y)

        for ik = 1:data.num_kpts
            G[:, ik] = vcat(vec(gradX[:, :, ik]), vec(gradY[:, :, ik]))
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

    tmp = omega(data, params, A, false, false)
    @info "Initial centers & spreads" tmp.centers tmp.spreads' sum(tmp.spreads)

    # need QR orthogonalization rather than SVD to preserve the sparsity structure of Y
    XYkManif = Optim.ProductManifold(Optim.Stiefel_SVD(), Optim.Stiefel_SVD(),
        (data.num_wann, data.num_wann), (data.num_bands, data.num_wann))
    XYManif = Optim.PowerManifold(XYkManif, (M,), (data.num_kpts,))

    # stepsize_mult = 1
    # step = 0.5/(4*8*p.wb)*(p.N1*p.N2*p.N3)*stepsize_mult
    # ls = LineSearches.Static(step)
    ls = Optim.HagerZhang()
    # ls = LineSearches.BackTracking()
    # meth = Optim.GradientDescent
    # meth = Optim.ConjugateGradient
    meth = Optim.LBFGS
    res = Optim.optimize(f, g!, XY0, meth(manifold=XYManif, linesearch=ls, m=params.history_size),
        Optim.Options(show_trace=true, iterations=params.max_iter,
            f_tol=params.omega_tol, g_tol=params.gradient_tol, allow_f_increases=true))
    display(res)
    XYmin = Optim.minimizer(res)

    Xmin, Ymin = XY_to_XY(data, XYmin)
    Amin = XY_to_A(data, Xmin, Ymin)

    tmp = omega(data, params, Amin, false, false)
    @info "Final centers & spreads" tmp.centers tmp.spreads' sum(tmp.spreads)

    return Amin
end

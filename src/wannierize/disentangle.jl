using LinearAlgebra
using Optim: Optim

function get_frozen_bands(
    E::AbstractMatrix{T}, dis_froz_max::T, dis_froz_min::T=-Inf
) where {T<:Real}
    n_bands, n_kpts = size(E)
    frozen_bands = falses(n_bands, n_kpts)

    # For each kpoint
    frozen_k = falses(n_bands)

    for ik in 1:n_kpts
        fill!(frozen_k, false)

        frozen_k .= (E[:, ik] .>= dis_froz_min) .& (E[:, ik] .<= dis_froz_max)
        frozen_bands[:, ik] = frozen_k
    end

    return frozen_bands
end

function set_frozen_degen!(
    frozen_bands::AbstractMatrix{Bool}, E::AbstractMatrix{T}, atol::T=1e-4
) where {T<:Real}
    atol <= 0 && error("atol must be positive")

    for ik in 1:n_kpts
        frozen_k = frozen_bands[:, ik]

        # if cluster of eigenvalues and count(frozen_k) > 0, take them all
        if degen && count(frozen_k) > 0
            ib = findlast(frozen_k)

            while ib < n_bands
                if E[ib + 1, ik] < E[ib, ik] + atol
                    ib += 1
                    frozen_k[ib] .= true
                else
                    break
                end
            end
        end

        frozen_bands[:, ik] .= frozen_k
    end

    return nothing
end

"""Make sure number of frozen bands at each kpoint <= n_wann"""
function check_frozen_bands(frozen_bands::AbstractMatrix{Bool}, n_wann::Int)
    n_bands, n_kpts = size(frozen_bands)
    n_wann > n_bands && error("n_wann > n_bands")

    for ik in 1:n_kpts
        frozen_k = frozen_bands[:, ik]

        if count(frozen_k) > n_wann
            error("Too many frozen bands")
        end
    end
end

"""
Set frozen bands according to two energy windows
"""
function set_frozen_win!(
    model::Model{T},
    dis_froz_max::T,
    dis_froz_min::T=-Inf;
    degen::Bool=false,
    degen_atol::T=1e-4,
) where {T<:Real}
    frozen_bands = get_frozen_bands(model.E, dis_froz_max, dis_froz_min)

    if degen
        set_frozen_degen!(frozen_bands, model.E, degen_atol)
    end

    check_frozen_bands(frozen_bands, model.n_wann)

    model.frozen_bands .= frozen_bands

    return nothing
end

"""
Get frozen bands according to projectability ∈ [0.0, 1.0]
"""
function get_frozen_proj(
    E::AbstractMatrix{T}, A::AbstractArray{Complex{T}}, dis_proj_max::T
) where {T<:Real}
    n_bands, n_kpts = size(E)
    frozen_bands = falses(n_bands, n_kpts)

    # For each kpoint
    frozen_k = falses(n_bands)

    for ik in 1:n_kpts
        fill!(frozen_k, false)

        # n_bands * n_wann
        Aₖ = A[:, :, ik]
        # projectability
        p = dropdims(real(sum(Aₖ .* conj(Aₖ); dims=2)); dims=2)
        # @debug "projectability" ik p

        frozen_k[p .>= dis_proj_max] .= true
        frozen_bands[:, ik] = frozen_k
    end

    return frozen_bands
end

"""
Set frozen bands according to projectability ∈ [0.0, 1.0]
"""
function set_frozen_proj!(
    model::Model{T}, dis_proj_max::T; degen::Bool=false, degen_atol::T=1e-4
) where {T<:Real}
    frozen_bands = get_frozen_bands_proj(model.E, model.A, dis_proj_max)

    if degen
        set_frozen_degen!(frozen_bands, model.E, degen_atol)
    end

    check_frozen_bands(frozen_bands, model.n_wann)

    model.frozen_bands .= frozen_bands

    return nothing
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
function orthonorm_freeze(A::Matrix{T}, frozen::BitVector) where {T<:Complex}
    n_bands, n_wann = size(A)
    non_frozen = .!frozen

    # Make sure Uf can fully represent frozen bands.
    # Uf = <ψ|g>, where |ψ> is Bloch wfcs, |g> is the guiding functions.
    # We do a Lowdin orthonormalization on Uf so that Uf * Uf' = I,
    # i.e. <ψ|g'><g'|ψ> = I -> |g'>s span the frozen |ψ>s.
    Uf = A[frozen, :]
    Uf = orthonorm_lowdin(Uf)
    # Uf = orthonorm_cholesky(Uf)

    # Remove Uf out of Ur, i.e. do not destroy frozen space
    # The projector of the frozen states represented on the |g> basis is
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
    # M = Ur' * Ur
    # S, U = eig(Hermitian(M))
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
    U, S, V = svd(Ur)
    atol = 1e-10
    @assert count(x -> x > atol, S) == n_wann - count(frozen)
    S[S .> atol] .= 1
    S[S .< atol] .= 0
    Ur = U * Diagonal(S) * V'

    B = similar(A)
    B[frozen, :] .= Uf
    B[non_frozen, :] .= Ur

    # Semiunitary
    @assert isapprox(B' * B, I; atol=atol)
    # Frozen
    @assert isapprox(B[frozen, :] * B[frozen, :]', I; atol=atol)
    # Independent
    @assert norm(Uf * Ur') < atol

    return B
end

# function max_projectability(A::Matrix{ComplexF64})
#     proj = A * A'
#     proj_ortho = I - proj
#     U, S, V = svd(proj_ortho)
#     n = abs(size(A, 1) - size(A, 2))
#     @assert count(S .> 1e-5) >= n
#     R = hcat(A, V[:, 1:n])
#     @assert R * R' ≈ I
#     return R' * A, R

#     U, S, V = svd(A)
#     return V' * A, V
#     proj = A * A'
#     # A * A' = V D V'  =>  (V'A) (A'V) = D
#     D, V = eigen(proj)
#     # sort eigenvalues in descending order
#     V = V[:, sortperm(D, rev = true)]
#     A_new = V' * A
#     @debug "projectability" real(diag(proj)') real(diag(A_new * A_new')')
#     return A_new, V
# end

# function maxproj_froz(A::Matrix{ComplexF64}, froz::BitVector)
#     @assert length(froz) == size(A, 1)
#     m, n = size(A)
#     D = zeros(ComplexF64, n + count(froz), m)
#     D[1:n, :] = transpose(A)
#     c = n + 1
#     for i = 1:length(froz)
#         if froz[i]
#             D[c, i] = 1
#             c += 1
#         end
#     end
#     U, S, V = svd(D)

#     nbasis = count(S .> 1e-5)
#     V_new = V[:, 1:nbasis]
#     proj = V_new * V_new'
#     proj_ortho = I - proj
#     U, S, V = svd(proj_ortho)
#     R = hcat(V_new, V[:, 1:m-nbasis])
#     @assert R * R' ≈ I

#     return R' * A, R
# end

"""
There are three formats: A, (X, Y), and XY stored contiguously in memory.
A is the format used in the rest of the code, XY is the format used in the optimizer,
(X, Y) is intermediate format.
A is n_bands * n_wann,
X is n_wann * n_wann,
Y is n_bands * n_wann and has first block equal to [I 0].
"""
function X_Y_to_A(X::Array{T,3}, Y::Array{T,3}) where {T<:Complex}
    n_bands, n_wann, n_kpts = size(Y)

    A = zeros(T, n_bands, n_wann, n_kpts)

    for ik in 1:n_kpts
        # check
        # idx_f = model.frozen_bands[:, ik]
        # idx_nf = .!idx_f
        # n_froz = count(idx_f)
        # @assert Y[:, :, ik]' * Y[:, :, ik] ≈ I
        # @assert X[:, :, ik]' * X[:, :, ik] ≈ I
        # @assert Y[idx_f, 1:n_froz, ik] ≈ I
        # @assert norm(Y[idx_nf, 1:n_froz, ik]) ≈ 0
        # @assert norm(Y[idx_f, n_froz+1:end, ik]) ≈ 0

        A[:, :, ik] = Y[:, :, ik] * X[:, :, ik]
        # @assert orthonorm_freeze(A[:, :, ik], frozen_bands[:, ik]) ≈ A[:, :, ik] rtol=1e-4
    end

    return A
end

"""
size(A) = n_bands * n_wann * n_kpts
size(frozen) = n_bands * n_kpts
"""
function A_to_X_Y(A::Array{T,3}, frozen::BitMatrix) where {T<:Complex}
    n_bands, n_wann, n_kpts = size(A)

    X = zeros(T, n_wann, n_wann, n_kpts)
    Y = zeros(T, n_bands, n_wann, n_kpts)

    for ik in 1:n_kpts
        idx_f = frozen[:, ik]
        idx_nf = .!idx_f
        n_froz = count(idx_f)

        Af = orthonorm_freeze(A[:, :, ik], idx_f)
        Uf = Af[idx_f, :]
        Ur = Af[idx_nf, :]

        # determine Y
        Y[idx_f, 1:n_froz, ik] = Matrix{T}(I, n_froz, n_froz)

        if n_froz != n_wann
            Pr = Ur * Ur'
            Pr = Hermitian((Pr + Pr') / 2)
            D, V = eigen(Pr) # sorted by increasing eigenvalue
            Y[idx_nf, (n_froz + 1):end, ik] = V[:, (end - n_wann + n_froz + 1):end]
        end

        # determine X
        X[:, :, ik] = orthonorm_lowdin(Y[:, :, ik]' * Af)

        @assert Y[:, :, ik]' * Y[:, :, ik] ≈ I
        @assert X[:, :, ik]' * X[:, :, ik] ≈ I
        @assert Y[idx_f, 1:n_froz, ik] ≈ I
        @assert norm(Y[idx_nf, 1:n_froz, ik]) ≈ 0
        @assert norm(Y[idx_f, (n_froz + 1):end, ik]) ≈ 0
        @assert Y[:, :, ik] * X[:, :, ik] ≈ Af
    end

    return X, Y
end

"""XY to (X, Y)"""
function XY_to_X_Y(XY::Matrix{T}, n_bands::Int, n_wann::Int) where {T<:Complex}
    n_kpts = size(XY, 2)

    X = zeros(T, n_wann, n_wann, n_kpts)
    Y = zeros(T, n_bands, n_wann, n_kpts)

    n = n_wann^2

    for ik in 1:n_kpts
        XYk = XY[:, ik]
        X[:, :, ik] = reshape(XYk[1:n], (n_wann, n_wann))
        Y[:, :, ik] = reshape(XYk[(n + 1):end], (n_bands, n_wann))
    end

    return X, Y
end

"""XY to (X, Y)"""
function X_Y_to_XY(X::Array{T,3}, Y::Array{T,3}) where {T<:Complex}
    n_bands, n_wann, n_kpts = size(Y)
    # n_wann, n_wann, n_kpts = size(X)

    n = n_wann^2
    XY = zeros(T, n + n_bands * n_wann, n_kpts)

    for ik in 1:n_kpts
        # XY[:, ik] = vcat(vec(X[:, :, ik]), vec(Y[:, :, ik]))
        XY[1:n, ik] = reshape(X[:, :, ik], n)
        XY[(n + 1):end, ik] = reshape(Y[:, :, ik], n_bands * n_wann)
    end

    return XY
end

function omega(
    bvectors::BVectors{FT},
    M::Array{Complex{FT},4},
    X::Array{Complex{FT},3},
    Y::Array{Complex{FT},3},
) where {FT<:Real}
    A = X_Y_to_A(X, Y)
    return omega(bvectors, M, A)
end

"""
size(M) = n_bands * n_bands * n_bvecs * n_kpts
size(Y) = n_wann * n_wann * n_kpts
size(Y) = n_bands * n_wann * n_kpts
size(frozen) = n_bands * n_kpts
"""
function omega_grad(
    bvectors::BVectors{FT},
    M::Array{Complex{FT},4},
    X::Array{Complex{FT},3},
    Y::Array{Complex{FT},3},
    frozen::BitMatrix,
) where {FT<:Real}
    n_kpts = size(Y, 3)

    A = X_Y_to_A(X, Y)

    G = omega_grad(bvectors, M, A)

    GX = zero(X)
    GY = zero(Y)

    for ik in 1:n_kpts
        idx_f = frozen[:, ik]
        n_froz = count(idx_f)

        GX[:, :, ik] = Y[:, :, ik]' * G[:, :, ik]
        GY[:, :, ik] = G[:, :, ik] * X[:, :, ik]'

        GY[idx_f, :, ik] .= 0
        GY[:, 1:n_froz, ik] .= 0

        # to compute the projected gradient: redundant, taken care of by the optimizer
        # function proj_stiefel(G,X)
        #     G .- X*((X'G .+ G'X)./2)
        # end
        # GX[:,:,ik] = proj_stiefel(GX[:,:,ik], X[:,:,ik])
        # GY[:,:,ik] = proj_stiefel(GY[:,:,ik], Y[:,:,ik])
    end

    return GX, GY
end

"""
Set grad of frozen bands to 0.

This is used in test.
"""
function zero_froz_grad!(G::Matrix, model::Wannier.Model)
    GX, GY = Wannier.XY_to_X_Y(G, model.n_bands, model.n_wann)
    for ik in 1:size(G, 2)
        idx_f = model.frozen_bands[:, ik]
        n_froz = count(idx_f)
        GY[idx_f, :, ik] .= 0
        GY[:, 1:n_froz, ik] .= 0
    end
    G .= Wannier.X_Y_to_XY(GX, GY)
    return nothing
end

function get_fg!_disentangle(model::Model)
    function f(XY)
        X, Y = XY_to_X_Y(XY, model.n_bands, model.n_wann)
        return omega(model.bvectors, model.M, X, Y).Ω
    end

    """size(G) == size(XY)"""
    function g!(G, XY)
        X, Y = XY_to_X_Y(XY, model.n_bands, model.n_wann)
        GX, GY = omega_grad(model.bvectors, model.M, X, Y, model.frozen_bands)

        n = model.n_wann^2

        for ik in 1:(model.n_kpts)
            # G[:, ik] .= vcat(vec(GX[:, :, ik]), vec(GY[:, :, ik]))
            G[1:n, ik] = vec(GX[:, :, ik])
            G[(n + 1):end, ik] = vec(GY[:, :, ik])
        end

        return nothing
    end

    return f, g!
end

function disentangle(
    model::Model{T};
    random_gauge::Bool=false,
    f_tol::T=1e-7,
    g_tol::T=1e-5,
    max_iter::Int=200,
    history_size::Int=20,
) where {T<:Real}
    n_bands = model.n_bands
    n_wann = model.n_wann
    n_kpts = model.n_kpts

    # initial X, Y
    if random_gauge
        X0 = zeros(Complex{T}, n_wann, n_wann, n_kpts)
        Y0 = zeros(Complex{T}, n_bands, n_wann, n_kpts)

        for ik in 1:n_kpts
            idx_f = model.frozen_bands[:, ik]
            idx_nf = .!idx_f
            n_froz = count(l_frozen)

            m = n_wann
            n = n_wann
            M = randn(T, m, n) + im * randn(T, m, n)
            X0[:, :, ik] = orthonorm_lowdin(M)

            Y0[idx_f, 1:n_froz, ik] = I
            m = n_bands - n_froz
            n = n_wann - n_froz
            N = randn(T, m, n) + im * randn(m, n)
            Y0[idx_nf, (n_froz + 1):n_wann, ik] = orthonorm_lowdin(N)
        end
    else
        X0, Y0 = A_to_X_Y(model.A, model.frozen_bands)
    end

    # compact storage
    XY0 = X_Y_to_XY(X0, Y0)

    # We have three storage formats:
    # (X, Y): n_wann * n_wann * n_kpts, n_bands * n_wann * n_kpts
    # A: n_bands * n_wann * n_kpts
    # XY: (n_wann * n_wann + n_bands * n_wann) * n_kpts
    f, g! = get_fg!_disentangle(model)

    Ωⁱ = omega(model.bvectors, model.M, model.A)
    @info "Initial spread"
    show(Ωⁱ)
    println("\n")

    # need QR orthogonalization rather than SVD to preserve the sparsity structure of Y
    XYkManif = Optim.ProductManifold(
        Optim.Stiefel_SVD(), Optim.Stiefel_SVD(), (n_wann, n_wann), (n_bands, n_wann)
    )
    XYManif = Optim.PowerManifold(XYkManif, (n_wann^2 + n_bands * n_wann,), (n_kpts,))

    # stepsize_mult = 1
    # step = 0.5/(4*8*p.wb)*(p.N1*p.N2*p.N3)*stepsize_mult
    # ls = LineSearches.Static(step)
    ls = Optim.HagerZhang()
    # ls = LineSearches.BackTracking()

    # meth = Optim.GradientDescent
    # meth = Optim.ConjugateGradient
    meth = Optim.LBFGS

    opt = Optim.optimize(
        f,
        g!,
        XY0,
        meth(; manifold=XYManif, linesearch=ls, m=history_size),
        Optim.Options(;
            show_trace=true,
            iterations=max_iter,
            f_tol=f_tol,
            g_tol=g_tol,
            allow_f_increases=true,
        ),
    )
    display(opt)

    XYmin = Optim.minimizer(opt)

    Xmin, Ymin = XY_to_X_Y(XYmin, n_bands, n_wann)
    Amin = X_Y_to_A(Xmin, Ymin)

    Ωᶠ = omega(model.bvectors, model.M, Amin)
    @info "Final spread"
    show(Ωᶠ)
    println("\n")

    return Amin
end

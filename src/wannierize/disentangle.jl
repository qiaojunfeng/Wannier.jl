using LinearAlgebra
using Optim: Optim

export disentangle

"""
    get_frozen_bands(E, dis_froz_max, dis_froz_min)

Generate a `BitMatrix` of frozen bands by checking the two frozen windows.

# Arguments
- `E`: the energy eigenvalues of the Hamiltonian
- `dis_froz_max`: the upper bound of the frozen window
- `dis_froz_min`: the lower bound of the frozen window

!!! note

    The `dis_froz_max` and `dis_froz_min` work similarly as `Wannier90`.
"""
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

"""
    set_frozen_degen!(frozen_bands, E, atol=1e-4)

Freeze bands which are degenerate with the `frozen_bands`.

In some cases, we might want to freeze the whole set of degenerated eigen vectors.

# Arguments
- `frozen_bands`: the `BitMatrix` of frozen bands
- `E`: the energy eigenvalues of the Hamiltonian
- `atol`: the tolerance of degeneracy
"""
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

"""
    check_frozen_bands(frozen_bands, n_wann)

Sanity check that the number of frozen bands at each kpoint <= `n_wann`.

# Arguments
- `frozen_bands`: the `BitMatrix` of frozen bands
- `n_wann`: the number of wannier functions
"""
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
    set_frozen_win!(model, dis_froz_max, dis_froz_min=-Inf; degen=false, degen_atol=1e-4)

Set frozen bands of the `Model` according to two energy windows.

# Arguments
- `model`: the `Model` to be set
- `dis_froz_max`: the upper bound of the frozen window
- `dis_froz_min`: the lower bound of the frozen window

# Keyword Arguments
- `degen`: whether to freeze the whole set of degenerated eigen vectors
- `degen_atol`: the tolerance of degeneracy
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
    get_frozen_proj(E, U, dis_proj_max)

Get frozen bands according to band projectability.

# Arguments
- `E`: the energy eigenvalues of the Hamiltonian
- `U`: the gauge rotation matrices
- `dis_proj_max`: the upper bound projectability.
    Bands with projectability >= `dis_proj_max` are frozen.

!!! note

    The band projectability for band ``n`` at kpoint ``\\bm{k}`` is calculated by
    ``p_{n \\bm{k}} = \\sum_{m=1}^{m=n_{wann}} | U_{nm \\bm{k}} |^2``,
    and usually each element ``p_{n \\bm{k}} \\in [0.0, 1.0]``.
    In such cases, the `dis_proj_max` is usually set to sth. like `0.9` to
    freeze high-projectability bands.
"""
function get_frozen_proj(
    E::AbstractMatrix{T}, U::AbstractArray{Complex{T}}, dis_proj_max::T
) where {T<:Real}
    n_bands, n_kpts = size(E)
    frozen_bands = falses(n_bands, n_kpts)

    # For each kpoint
    frozen_k = falses(n_bands)

    for ik in 1:n_kpts
        fill!(frozen_k, false)

        # n_bands * n_wann
        Uₖ = U[:, :, ik]
        # projectability
        p = dropdims(real(sum(Uₖ .* conj(Uₖ); dims=2)); dims=2)
        # @debug "projectability" ik p

        frozen_k[p .>= dis_proj_max] .= true
        frozen_bands[:, ik] = frozen_k
    end

    return frozen_bands
end

"""
    set_frozen_proj!(model, dis_proj_max; degen=false, degen_atol=1e-4)

Set frozen bands of the `Model` according to projectability.

# Arguments
- `model`: the `Model` to be set
- `dis_proj_max`: the upper bound projectability.
    Bands with projectability >= `dis_proj_max` are frozen.

# Keyword Arguments
- `degen`: whether to freeze the whole set of degenerated eigen vectors
- `degen_atol`: the tolerance of degeneracy
"""
function set_frozen_proj!(
    model::Model{T}, dis_proj_max::T; degen::Bool=false, degen_atol::T=1e-4
) where {T<:Real}
    frozen_bands = get_frozen_proj(model.E, model.U, dis_proj_max)

    if degen
        set_frozen_degen!(frozen_bands, model.E, degen_atol)
    end

    check_frozen_bands(frozen_bands, model.n_wann)

    model.frozen_bands .= frozen_bands

    return nothing
end

"""
    orthonorm_freeze(U, frozen)

Normalize and freeze a block of a matrix.

Conditions:
- Block form:  `U = vcat(Uf, Ur)`
- Semiunitary
  - `U' * U = I`
  - `Uf' * Uf + Ur' * Ur = I`
- Frozen:      `Uf * Uf' = I`
- Also:        `Uf * Ur' = 0`

Strategy:
1. orthogonalize `Uf`
2. project `Uf` out of `Ur`
3. orthogonalize the range of `Ur`

# Arguments
- `U`: the matrix to be orthonormalized and frozen
- `frozen`: the `BitVector` specifying which bands are frozen
"""
function orthonorm_freeze(U::Matrix{T}, frozen::BitVector) where {T<:Complex}
    n_bands, n_wann = size(U)
    non_frozen = .!frozen

    # Make sure Uf can fully represent frozen bands.
    # Uf = <ψ|g>, where |ψ> is Bloch wfcs, |g> is the guiding functions.
    # We do a Lowdin orthonormalization on Uf so that Uf * Uf' = I,
    # i.e. <ψ|g'><g'|ψ> = I -> |g'>s span the frozen |ψ>s.
    Uf = U[frozen, :]
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
    Ur = U[non_frozen, :]
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
    A, S, B = svd(Ur)
    atol = 1e-10
    @assert count(x -> x > atol, S) == n_wann - count(frozen)
    S[S .> atol] .= 1
    S[S .< atol] .= 0
    Ur = A * Diagonal(S) * B'

    V = similar(U)
    V[frozen, :] .= Uf
    V[non_frozen, :] .= Ur

    # Semiunitary
    @assert isapprox(V' * V, I; atol=atol)
    # Frozen
    @assert isapprox(V[frozen, :] * V[frozen, :]', I; atol=atol)
    # Independent
    @assert norm(Uf * Ur') < atol

    return V
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
    X_Y_to_U(X::Array{T,3}, Y::Array{T,3})

Convert the `(X, Y)` layout to the `U` layout.

There are three formats: `U`, `(X, Y)`, and `XY` stored contiguously in memory.
For each kpoint,
- `U`: `size(U) = (n_bands, n_wann)`, the format used in the rest of the code
- `(X, Y)`: `size(X) = (n_wann, n_wann)`, `size(Y) = (n_bands, n_wann)`, intermediate format
- `XY`: this is the format used in the optimizer
"""
function X_Y_to_U(X::AbstractArray{T,3}, Y::AbstractArray{T,3}) where {T<:Complex}
    n_bands, n_wann, n_kpts = size(Y)

    U = zeros(T, n_bands, n_wann, n_kpts)

    for ik in 1:n_kpts
        # sanity check, should be true but comment out to be faster
        # idx_f = model.frozen_bands[:, ik]
        # idx_nf = .!idx_f
        # n_froz = count(idx_f)
        # @assert Y[:, :, ik]' * Y[:, :, ik] ≈ I
        # @assert X[:, :, ik]' * X[:, :, ik] ≈ I
        # @assert Y[idx_f, 1:n_froz, ik] ≈ I
        # @assert norm(Y[idx_nf, 1:n_froz, ik]) ≈ 0
        # @assert norm(Y[idx_f, n_froz+1:end, ik]) ≈ 0

        U[:, :, ik] = Y[:, :, ik] * X[:, :, ik]
        # @assert orthonorm_freeze(U[:, :, ik], frozen_bands[:, ik]) ≈ U[:, :, ik] rtol=1e-4
    end

    return U
end

"""
    U_to_X_Y(U::Array{T,3}, frozen::BitMatrix) where {T<:Complex}

Convert the `U` layout to the `(X, Y)` layout.

See also [`X_Y_to_U`](@ref).

# Arguments
- `U`: `n_bands * n_wann * n_kpts`
- `frozen`: `n_bands * n_kpts`
"""
function U_to_X_Y(U::AbstractArray{T,3}, frozen::BitMatrix) where {T<:Complex}
    n_bands, n_wann, n_kpts = size(U)

    X = zeros(T, n_wann, n_wann, n_kpts)
    Y = zeros(T, n_bands, n_wann, n_kpts)

    for ik in 1:n_kpts
        idx_f = frozen[:, ik]
        idx_nf = .!idx_f
        n_froz = count(idx_f)

        Af = orthonorm_freeze(U[:, :, ik], idx_f)
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

"""
    XY_to_X_Y(XY::Matrix{T}, n_bands::Int, n_wann::Int)

Convert the `XY` layout to the `(X, Y)` layout.

See also [`X_Y_to_U`](@ref).

# Arguments
- `XY`: `n_bands * n_wann * n_kpts` contiguous array
- `n_bands`: number of bands, to be used to reshape `XY`
- `n_wann`: number of wannier functions, to be used to reshape `XY`
"""
function XY_to_X_Y(XY::AbstractMatrix{T}, n_bands::Int, n_wann::Int) where {T<:Complex}
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

"""
    X_Y_to_XY(X::Array{T,3}, Y::Array{T,3}) where {T<:Complex}

Convert the `(X, Y)` layout to the `XY` layout.

See also [`X_Y_to_U`](@ref).
"""
function X_Y_to_XY(X::AbstractArray{T,3}, Y::AbstractArray{T,3}) where {T<:Complex}
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

"""
    omega(bvectors, M, X, Y)

Compute WF spread in the `(X, Y)` layout.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `X`: `n_wann * n_wann * n_kpts` array
- `Y`: `n_bands * n_wann * n_kpts` array
"""
function omega(
    bvectors::BVectors{FT},
    M::Array{Complex{FT},4},
    X::Array{Complex{FT},3},
    Y::Array{Complex{FT},3},
) where {FT<:Real}
    U = X_Y_to_U(X, Y)
    return omega(bvectors, M, U)
end

"""
    omega(model, X, Y)

Compute WF spread in the `(X, Y)` layout.

# Arguments
- `model`: `Model`
- `X`: `n_wann * n_wann * n_kpts` array
- `Y`: `n_bands * n_wann * n_kpts` array
"""
function omega(
    model::Model{FT}, X::Array{Complex{FT},3}, Y::Array{Complex{FT},3}
) where {FT<:Real}
    return omega(model.bvectors, model.M, X, Y)
end

"""
    omega_grad(bvectors, M, X, Y, frozen)

Compute gradient of WF spread in the `(X, Y)` layout.

# Arguments
- `bvectors`: bvecoters
- `M`: `n_bands * n_bands * * n_bvecs * n_kpts` overlap array
- `X`: `n_wann * n_wann * n_kpts` array
- `Y`: `n_bands * n_wann * n_kpts` array
- `frozen`: `n_bands * n_kpts` array for frozen bands
"""
function omega_grad(
    bvectors::BVectors{FT},
    M::Array{Complex{FT},4},
    X::Array{Complex{FT},3},
    Y::Array{Complex{FT},3},
    frozen::BitMatrix,
) where {FT<:Real}
    n_kpts = size(Y, 3)

    U = X_Y_to_U(X, Y)

    G = omega_grad(bvectors, M, U)

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
    zero_froz_grad!(G, model)

Set gradient of frozen bands to 0.

This is used in test.

# Arguments
- `G`: gradient of the spread, in `XY` layout
- `model`: model
"""
function zero_froz_grad!(G::Matrix, model::Model)
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

"""
    get_fg!_disentangle(model::Model)

Return a tuple of two functions `(f, g!)` for spread and gradient, respectively.
"""
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

"""
    disentangle(model; random_gauge=false, f_tol=1e-7, g_tol=1e-5, max_iter=200, history_size=3)

Run disentangle on the `Model`.

# Arguments
- `model`: model

# Keyword arguments
- `random_gauge`: use random `U` matrices as initial guess
- `f_tol`: tolerance for spread convergence
- `g_tol`: tolerance for gradient convergence
- `max_iter`: maximum number of iterations
- `history_size`: history size of LBFGS
"""
function disentangle(
    model::Model{T};
    random_gauge::Bool=false,
    f_tol::T=1e-7,
    g_tol::T=1e-5,
    max_iter::Int=200,
    history_size::Int=3,
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
        X0, Y0 = U_to_X_Y(model.U, model.frozen_bands)
    end

    # compact storage
    XY0 = X_Y_to_XY(X0, Y0)

    # We have three storage formats:
    # (X, Y): n_wann * n_wann * n_kpts, n_bands * n_wann * n_kpts
    # U: n_bands * n_wann * n_kpts
    # XY: (n_wann * n_wann + n_bands * n_wann) * n_kpts
    f, g! = get_fg!_disentangle(model)

    Ωⁱ = omega(model)
    @info "Initial spread"
    show(Ωⁱ)
    println("\n")

    Ωⁱ = omega(model, X0, Y0)
    @info "Initial spread (with states freezed)"
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
    Umin = X_Y_to_U(Xmin, Ymin)

    Ωᶠ = omega(model, Umin)
    @info "Final spread"
    show(Ωᶠ)
    println("\n")

    return Umin
end

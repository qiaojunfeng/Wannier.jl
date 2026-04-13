"""
    X_Y_to_U(X::Array{T,3}, Y::Array{T,3})

Convert the `(X, Y)` layout to the `U` layout.

There are three formats: `U`, `(X, Y)`, and `XY` stored contiguously in memory.
For each kpoint,
- `U`: `size(U) = (n_bands, n_wann)`, the format used in the rest of the code
- `(X, Y)`: `size(X) = (n_wann, n_wann)`, `size(Y) = (n_bands, n_wann)`, intermediate format
- `XY`: this is the format used in the optimizer
"""
function X_Y_to_U(
        X::AbstractVector{<:AbstractMatrix{T}}, Y::AbstractVector{<:AbstractMatrix{T}}
    ) where {T <: Complex}
    n_bands, n_wann = size(Y[1])
    n_kpts = length(Y)

    U = [zeros(T, n_bands, n_wann) for i in 1:n_kpts]
    return X_Y_to_U!(U, X, Y)
end

function X_Y_to_U!(U::AbstractVector, X::AbstractVector, Y::AbstractVector)
    @inbounds for (u, y, x) in zip(U, Y, X)
        mul!(u, y, x)
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
function U_to_X_Y(
        U::AbstractVector{<:AbstractMatrix{T}}, frozen::Vector{BitVector}
    ) where {T <: Complex}
    nkpts = length(U)
    nbands, nwann = size(U[1])

    X = [zeros(T, nwann, nwann) for i in 1:nkpts]
    Y = [zeros(T, nbands, nwann) for i in 1:nkpts]

    @inbounds for ik in 1:nkpts
        idx_f = frozen[ik]
        idx_nf = .!idx_f
        n_froz = count(idx_f)

        Af = orthonorm_freeze(U[ik], idx_f)
        Uf = Af[idx_f, :]
        Ur = Af[idx_nf, :]

        # determine Y
        Y[ik][idx_f, 1:n_froz] .= Matrix{T}(I, n_froz, n_froz)

        if n_froz != nwann
            Pr = Ur * Ur'
            Pr = Hermitian((Pr + Pr') / 2)
            D, V = eigen(Pr) # sorted by increasing eigenvalue
            Y[ik][idx_nf, (n_froz + 1):end] .= V[:, (end - nwann + n_froz + 1):end]
        end

        # determine X
        X[ik] = orthonorm_lowdin(Y[ik]' * Af)
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
function XY_to_X_Y(XY::AbstractMatrix{T}, nbands::Int, nwann::Int) where {T <: Complex}
    nkpts = size(XY, 2)

    X = [zeros(T, nwann, nwann) for i in 1:nkpts]
    Y = [zeros(T, nbands, nwann) for i in 1:nkpts]
    return XY_to_X_Y!(X, Y, XY)
end

function XY_to_X_Y!(X::AbstractVector, Y::AbstractVector, XY::AbstractMatrix)
    n_wann2 = size(X[1], 1)^2
    @inbounds for (ik, (x, y)) in enumerate(zip(X, Y))
        for i in eachindex(x)
            x[i] = XY[i, ik]
        end
        for i in eachindex(y)
            y[i] = XY[n_wann2 + i, ik]
        end
    end
    return X, Y
end

"""
    X_Y_to_XY(X::Array{T,3}, Y::Array{T,3}) where {T<:Complex}

Convert the `(X, Y)` layout to the `XY` layout.

See also [`X_Y_to_U`](@ref).
"""
function X_Y_to_XY(
        X::AbstractVector{<:AbstractMatrix{T}}, Y::AbstractVector{<:AbstractMatrix{T}}
    ) where {T <: Complex}
    nkpts = length(Y)
    nbands, nwann = size(Y[1])
    n = nwann^2
    XY = zeros(T, n + nbands * nwann, nkpts)
    return X_Y_to_XY!(XY, X, Y)
end

function X_Y_to_XY!(XY::AbstractMatrix, X::AbstractVector, Y::AbstractVector)
    n = length(X[1])
    @inbounds for (ik, (x, y)) in enumerate(zip(X, Y))
        for i in eachindex(x)
            XY[i, ik] = x[i]
        end

        for i in eachindex(y)
            XY[n + i, ik] = y[i]
        end
    end
    return XY
end

@doc raw"""
    GU_to_GX_GY(G, X, Y, frozen)

Compute dΩ/dX and dΩ/dY from dΩ/dU.

Acutally they are the conjugate gradients, e.g., ``\frac{d \Omega}{d U^*}``.

# Arguments
- `G`: `n_bands * n_wann * n_kpts` array for gradient dΩ/dU
- `X`: `n_wann * n_wann * n_kpts` array for X
- `Y`: `n_bands * n_wann * n_kpts` array for Y
- `frozen`: `n_bands * n_kpts` BitMatrix for frozen bands
"""
function GU_to_GX_GY(
        G::Array{T, 3}, X::Vector{Matrix{T}}, Y::Vector{Matrix{T}}, frozen::Vector
    ) where {T}
    n_kpts = length(X)
    GX = [zeros(T, size(X[1])) for i in 1:n_kpts]
    GY = [zeros(T, size(Y[1])) for i in 1:n_kpts]

    @inbounds for ik in 1:n_kpts
        idx_f = frozen[ik]
        n_froz = count(idx_f)

        mul!(GX[ik], Y[ik]', G[:, :, ik])
        mul!(GY[ik], G[:, :, ik], X[ik]')

        GY[ik][idx_f, :] .= 0
        GY[ik][:, 1:n_froz] .= 0
    end

    return GX, GY
end

function GU_to_GX_GY(
        G::Vector, X::Vector{Matrix{T}}, Y::Vector{Matrix{T}}, frozen::Vector
    ) where {T}
    n_kpts = length(X)
    GX = [zeros(T, size(X[1])) for i in 1:n_kpts]
    GY = [zeros(T, size(Y[1])) for i in 1:n_kpts]

    @inbounds for ik in 1:n_kpts
        idx_f = frozen[ik]
        n_froz = count(idx_f)

        mul!(GX[ik], Y[ik]', G[ik])
        mul!(GY[ik], G[ik], X[ik]')

        GY[ik][idx_f, :] .= 0
        GY[ik][:, 1:n_froz] .= 0
    end

    return GX, GY
end

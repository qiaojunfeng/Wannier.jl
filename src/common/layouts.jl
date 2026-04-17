"""
    X_Y_to_U(X, Y)

Convert the `(X, Y)` layout to the `U` layout.

There are three formats: `U`, `(X, Y)`, and `XY` stored contiguously in memory.
For each kpoint,
- `U`: `size(U) = (n_bands, n_wann, n_kpts)`, the format used in the rest of the code
- `(X, Y)`: `size(X) = (n_wann, n_wann, n_kpts)`, `size(Y) = (n_bands, n_wann, n_kpts)`, intermediate format
- `XY`: this is the format used in the optimizer
"""
function X_Y_to_U(X::AbstractArray{T, 3}, Y::AbstractArray{T, 3}) where {T <: Complex}
    n_bands = size(Y, 1)
    n_wann = size(Y, 2)
    n_kpts = size(Y, 3)

    U = zeros(T, n_bands, n_wann, n_kpts)
    return X_Y_to_U!(U, X, Y)
end

function X_Y_to_U!(U::AbstractArray{T, 3}, X::AbstractArray{T, 3}, Y::AbstractArray{T, 3}) where {T}
    @inbounds for ik in axes(U, 3)
        mul!(view(U, :, :, ik), view(Y, :, :, ik), view(X, :, :, ik))
    end
    return U
end

"""
    U_to_X_Y(U, frozen)

Convert the `U` layout to the `(X, Y)` layout.

See also [`X_Y_to_U`](@ref).

# Arguments
- `U`: `n_bands × n_wann × n_kpts`
- `frozen`: `n_bands × n_kpts` BitMatrix
"""
function U_to_X_Y(U::AbstractArray{T, 3}, frozen::AbstractMatrix{Bool}) where {T <: Complex}
    nkpts = size(U, 3)
    nbands, nwann = size(U, 1), size(U, 2)

    X = zeros(T, nwann, nwann, nkpts)
    Y = zeros(T, nbands, nwann, nkpts)

    @inbounds for ik in 1:nkpts
        idx_f = view(frozen, :, ik)
        idx_nf = .!idx_f
        n_froz = count(idx_f)

        Af = orthonorm_freeze(view(U, :, :, ik), idx_f)
        Uf = Af[idx_f, :]
        Ur = Af[idx_nf, :]

        # determine Y
        Y[idx_f, 1:n_froz, ik] .= Matrix{T}(I, n_froz, n_froz)

        if n_froz != nwann
            Pr = Ur * Ur'
            Pr = Hermitian((Pr + Pr') / 2)
            D, V = eigen(Pr) # sorted by increasing eigenvalue
            Y[idx_nf, (n_froz + 1):end, ik] .= V[:, (end - nwann + n_froz + 1):end]
        end

        # determine X
        X[:, :, ik] .= orthonorm_lowdin(view(Y, :, :, ik)' * Af)
    end

    return X, Y
end

"""
    XY_to_X_Y(XY, n_bands, n_wann)

Convert the `XY` layout to the `(X, Y)` layout.

See also [`X_Y_to_U`](@ref).

# Arguments
- `XY`: `n_bands * n_wann * n_kpts` contiguous array
- `n_bands`: number of bands, to be used to reshape `XY`
- `n_wann`: number of wannier functions, to be used to reshape `XY`
"""
function XY_to_X_Y(XY::AbstractMatrix{T}, nbands::Int, nwann::Int) where {T <: Complex}
    nkpts = size(XY, 2)

    X = zeros(T, nwann, nwann, nkpts)
    Y = zeros(T, nbands, nwann, nkpts)
    return XY_to_X_Y!(X, Y, XY)
end

function XY_to_X_Y!(X::AbstractArray{T, 3}, Y::AbstractArray{T, 3}, XY::AbstractMatrix) where {T}
    n_wann2 = size(X, 1)^2
    @inbounds for ik in axes(X, 3)
        xk = view(X, :, :, ik)
        yk = view(Y, :, :, ik)
        for i in eachindex(xk)
            xk[i] = XY[i, ik]
        end
        for i in eachindex(yk)
            yk[i] = XY[n_wann2 + i, ik]
        end
    end
    return X, Y
end

"""
    X_Y_to_XY(X, Y)

Convert the `(X, Y)` layout to the `XY` layout.

See also [`X_Y_to_U`](@ref).
"""
function X_Y_to_XY(X::AbstractArray{T, 3}, Y::AbstractArray{T, 3}) where {T <: Complex}
    nkpts = size(Y, 3)
    nbands, nwann = size(Y, 1), size(Y, 2)
    n = nwann^2
    XY = zeros(T, n + nbands * nwann, nkpts)
    return X_Y_to_XY!(XY, X, Y)
end

function X_Y_to_XY!(XY::AbstractMatrix, X::AbstractArray{T, 3}, Y::AbstractArray{T, 3}) where {T}
    @inbounds for ik in axes(X, 3)
        xk = view(X, :, :, ik)
        yk = view(Y, :, :, ik)
        n = length(xk)
        for i in eachindex(xk)
            XY[i, ik] = xk[i]
        end
        for i in eachindex(yk)
            XY[n + i, ik] = yk[i]
        end
    end
    return XY
end

@doc raw"""
    GU_to_GX_GY(G, X, Y, frozen)

Compute dΩ/dX and dΩ/dY from dΩ/dU.

Acutally they are the conjugate gradients, e.g., ``\frac{d \Omega}{d U^*}``.

# Arguments
- `G`: `n_bands × n_wann × n_kpts` array for gradient dΩ/dU
- `X`: `n_wann × n_wann × n_kpts` array for X
- `Y`: `n_bands × n_wann × n_kpts` array for Y
- `frozen`: `n_bands × n_kpts` BitMatrix for frozen bands
"""
function GU_to_GX_GY(
        G::AbstractArray{T, 3},
        X::AbstractArray{T, 3},
        Y::AbstractArray{T, 3},
        frozen::AbstractMatrix{Bool},
    ) where {T}
    n_kpts = size(X, 3)
    GX = zeros(T, size(X)...)
    GY = zeros(T, size(Y)...)

    @inbounds for ik in 1:n_kpts
        idx_f = view(frozen, :, ik)
        n_froz = count(idx_f)

        mul!(view(GX, :, :, ik), view(Y, :, :, ik)', view(G, :, :, ik))
        mul!(view(GY, :, :, ik), view(G, :, :, ik), view(X, :, :, ik)')

        view(GY, idx_f, :, ik) .= 0
        view(GY, :, 1:n_froz, ik) .= 0
    end

    return GX, GY
end

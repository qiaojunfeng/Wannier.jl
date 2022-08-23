
struct RGrid{T<:Real,XT<:AbstractArray3,YT<:AbstractArray3,ZT<:AbstractArray3}
    # spanning vectors, 3 * 3, each column is a spanning vector
    basis::Mat3{T}

    # usually these are LazyGrids
    # x fractional coordinate, nx * ny * nz
    X::XT
    # y fractional coordinate, nx * ny * nz
    Y::YT
    # z fractional coordinate, nx * ny * nz
    Z::ZT
end

function RGrid(basis::AbstractMatrix, X, Y, Z)
    size(X) == size(Y) == size(Z) || error("X, Y, Z must have the same size")
    return RGrid(Mat3(basis), X, Y, Z)
end

@doc """
Get origin of the RGrid.
"""
function origin(rgrid::RGrid)
    O = [rgrid.X[1, 1, 1], rgrid.Y[1, 1, 1], rgrid.Z[1, 1, 1]]
    origin = rgrid.basis * O
    return origin
end

@doc """
Get the span vectors of the RGrid.
"""
function span_vectors(rgrid::RGrid)
    O = [rgrid.X[1, 1, 1], rgrid.Y[1, 1, 1], rgrid.Z[1, 1, 1]]
    v1 = [rgrid.X[end, 1, 1], rgrid.Y[end, 1, 1], rgrid.Z[end, 1, 1]] - O
    v2 = [rgrid.X[1, end, 1], rgrid.Y[1, end, 1], rgrid.Z[1, end, 1]] - O
    v3 = [rgrid.X[1, 1, end], rgrid.Y[1, 1, end], rgrid.Z[1, 1, end]] - O
    # to cartesian
    v1 = rgrid.basis * v1
    v2 = rgrid.basis * v2
    v3 = rgrid.basis * v3
    # each column is a vector
    spanvec = hcat(v1, v2, v3)
    return spanvec
end

@doc """Return X, Y, Z in cartesian coordinates"""
function cartesianize_xyz(rgrid::RGrid)
    XYZ = vcat(reshape(rgrid.X, :)', reshape(rgrid.Y, :)', reshape(rgrid.Z, :)')
    XYZᶜ = rgrid.basis * XYZ
    Xᶜ = reshape(XYZᶜ[1, :], size(rgrid.X)...)
    Yᶜ = reshape(XYZᶜ[2, :], size(rgrid.X)...)
    Zᶜ = reshape(XYZᶜ[3, :], size(rgrid.X)...)
    return Xᶜ, Yᶜ, Zᶜ
end

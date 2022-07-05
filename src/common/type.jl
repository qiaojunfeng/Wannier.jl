using StaticArrays


# Frequently-used array types
const Mat3{T} = SMatrix{3, 3, T, 9} where T
const Vec3{T} = SVector{3, T} where T

# High-symmetry kpoint, label => coordinates
const StrVec3{T} = Pair{String,Vec3{T}} where T

# Kpath for band, Vector: list of segments, Segment: start and end point
const Kpath{T} = Vector{Tuple{StrVec3{T}, StrVec3{T}}} where T

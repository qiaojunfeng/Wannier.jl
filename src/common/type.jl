using StaticArrays

# Frequently-used array types
const Mat3{T} = SMatrix{3,3,T,9} where {T}
const Vec3{T} = SVector{3,T} where {T}

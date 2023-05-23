using StaticArrays

# Frequently-used array types
const Mat3{T} = SMatrix{3,3,T,9} where {T}
const Vec3{T} = SVector{3,T} where {T}

AbstractArray3{T} = AbstractArray{T,3} where {T}

# some small type piracy?
Base.zeros(m::AbstractArray{T}) where {T} = fill!(similar(m), zero(T))
Base.parse(::Type{Vec3{T}}, s) where {T} = Vec3(parse.(T, s)...)

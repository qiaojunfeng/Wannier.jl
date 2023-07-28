using StaticArrays
using WannierIO: Vec3, Mat3

const AbstractArray3{T} = AbstractArray{T,3} where {T}

# some small type piracy?
# TODO: check if these are actually needed
# Base.zeros(m::AbstractArray{T}) where {T} = fill!(similar(m), zero(T))
# Base.parse(::Type{Vec3{T}}, s) where {T} = Vec3(parse.(T, s)...)

"""
mutable 3 x 3 matrix type.

For lattice and recip_lattice.
"""
const MMat3{T} = MMatrix{3,3,T,9} where {T}

"""
Mutable length-3 vector type.

For atom posistions, kpoints, etc.
"""
const MVec3{T} = MVector{3,T} where {T}

"""
`Vector{Vector}` -> `MMat3`
"""
MMat3(A::AbstractVector) = MMat3(reduce(hcat, A))

"""
`MMat3` -> `Vec3{Vec3}`
"""
MVec3(A::MMat3) = MVec3(eachcol(A))

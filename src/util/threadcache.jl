using Base.Threads

struct ThreadCache{T}
    caches::Vector{T}
    ThreadCache(orig::T) where {T} =
        new{T}([deepcopy(orig) for i = 1:nthreads()])
end

@inline cache(t::ThreadCache) =
    t.caches[threadid()]

for f in (:getindex, :setindex!, :copyto!, :size, :length, :iterate, :sum, :view, :fill!)
    @eval Base.$f(t::ThreadCache{<:AbstractArray}, i...) = Base.$f(cache(t), i...)
end

for op in (:+, :-, :*, :/)
    @eval Base.$op(t::ThreadCache{T}, v::T) where {T} = $op(cache(t), v)
    @eval Base.$op(v::T, t::ThreadCache{T}) where {T} = $op(v, cache(t))
end

fillall!(t::ThreadCache{<:AbstractArray{T}}, v::T)   where {T} =
    fill!.(t.caches, (v,))

Base.reduce(op, t::ThreadCache; kwargs...) = reduce(op, t.caches; kwargs...)

LinearAlgebra.mul!(t1::ThreadCache{T}, v::T, t2::ThreadCache{T}) where {T<:AbstractArray} =
    mul!(cache(t1), v, cache(t2))

LinearAlgebra.mul!(t1::T, v::T, t2::ThreadCache{T}) where {T<:AbstractArray} =
    mul!(t1, v, cache(t2))

LinearAlgebra.mul!(t1::ThreadCache{T}, t2::ThreadCache{T}, t3::ThreadCache{T}) where {T<:AbstractArray} =
    mul!(cache(t1), cache(t2), cache(t3))

LinearAlgebra.adjoint!(t1::ThreadCache{T}, v::T) where {T} =
    adjoint!(cache(t1), v)

Base.ndims(::Type{ThreadCache{T}}) where {T<:AbstractArray} =
    ndims(T)
Base.Broadcast.broadcastable(tc::ThreadCache{<:AbstractArray}) =
    cache(tc)


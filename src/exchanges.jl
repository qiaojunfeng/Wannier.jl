# Magnetic matrices
using Base: @propagate_inbounds
using AtomsBase
using LinearAlgebra.LAPACK: BlasInt
import LinearAlgebra: eigen, eigen!

div1(x, y) = div(x - 1, y) + 1

phases(kpoints::Vector{<:Vec3}, R::Vec3) = exp.(-2im * π .* dot.(kpoints, (R,)))

struct ThreadCache{T}
    caches::Vector{T}
    ThreadCache(orig::T) where {T} =
        new{T}([deepcopy(orig) for i = 1:Threads.nthreads()])
end

@inline cache(t::ThreadCache) =
    t.caches[Threads.threadid()]

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

abstract type Spin end
struct Up <: Spin end
struct Down <: Spin end
   
"Represents a magnetic Hamiltonian matrix with the block structure [up updown;downup down]"
abstract type AbstractMagneticMatrix{T} <: AbstractMatrix{T} end

data(m::AbstractMatrix) = m
data(m::AbstractMagneticMatrix) = m.data

Base.similar(::Type{M}, i::NTuple{2,Int}) where {M <: AbstractMagneticMatrix} =
    M(Matrix{M.parameters[1]}(undef, i))

for f in (:length, :size, :setindex!, :elsize)
    @eval @inline @propagate_inbounds Base.$f(c::AbstractMagneticMatrix, args...) =
    Base.$f(c.data, args...)
end

Base.pointer(c::AbstractMagneticMatrix, i::Integer) = pointer(c.data, i)

"Magnetic block dimensions"
blockdim(c::AbstractMatrix) = div(size(data(c), 2), 2)

up(c::AbstractMatrix) =   (d = blockdim(c); view(data(c), 1:d, 1:d))
down(c::AbstractMatrix) = (d = blockdim(c); r = d + 1:2 * d; view(data(c), r, r))

# Standard getindex behavior
for f in (:view, :getindex)
    @eval @inline @propagate_inbounds Base.$f(c::AbstractMagneticMatrix, args...) =
        $f(c.data, args...)
    @eval @inline @propagate_inbounds Base.$f(c::AbstractMagneticMatrix, r::Union{Colon, AbstractUnitRange}, i::Int) =
        MagneticVector(Base.$f(c.data, r, i))
    @eval @inline @propagate_inbounds Base.$f(c::AbstractMagneticMatrix, i::Int, r::Union{Colon, AbstractUnitRange}) =
        MagneticVector(Base.$f(c.data, i, r))
end

Base.similar(c::M, args::AbstractUnitRange...) where {M <: AbstractMagneticMatrix} =
    M(similar(c.data), args...)

Base.iterate(c::AbstractMagneticMatrix, args...) = iterate(c.data, args...)

"""
    ColinMatrix{T, M <: AbstractMatrix{T}} <: AbstractMagneticMatrix{T}

Defines a Hamiltonian Matrix with [up zeros; zeros down] structure.
It is internally only storing the up and down block.
"""
struct ColinMatrix{T,M <: AbstractMatrix{T}} <: AbstractMagneticMatrix{T}
    data::M
end

function ColinMatrix(up::AbstractMatrix, down::AbstractMatrix)
    @assert size(up) == size(down)
    return ColinMatrix([up down])
end

Base.Array(c::ColinMatrix{T}) where T =
    (d = blockdim(c); [c[Up()] zeros(T, d, d); zeros(T, d, d) c[Down()]])

down(c::ColinMatrix) = (d = blockdim(c); view(c.data, 1:d, d + 1:2 * d))
blockdim(c::ColinMatrix) = size(c.data, 1)

function LinearAlgebra.diag(c::ColinMatrix)
    d = blockdim(c)
    r = LinearAlgebra.diagind(d, d)
    [c[r];c[r.+last(r)]]
end

"""
    NonColinMatrix{T, M <: AbstractMatrix{T}} <: AbstractMagneticMatrix{T}


Defines a Hamiltonian Matrix with [up updown;downup down] structure.
Since atomic projections w.r.t spins are defined rather awkwardly in Wannier90 for exchange calculations,
i.e. storing the up-down parts of an atom sequentially,
a NonColinMatrix reshuffles the entries of a matrix such that it follows the above structure. 
"""
struct NonColinMatrix{T,M <: AbstractMatrix{T}} <: AbstractMagneticMatrix{T}
    data::M
end


"Reshuffles standard Wannier90 up-down indices to the ones for the structure of a NonColinMatrix."
function Base.convert(::Type{NonColinMatrix}, m::M) where {M <: AbstractMatrix}
    @assert iseven(size(m, 1)) "Error, dimension of the supplied matrix is odd, i.e. it does not contain both spin components."
    data = similar(m)
    d    = blockdim(m)
    for i in 1:2:size(m, 1), j in 1:2:size(m, 2) 
        up_id1 = div1(i, 2) 
        up_id2 = div1(j, 2) 
        data[up_id1, up_id2] = m[i, j] 
        data[up_id1 + d, up_id2] = m[i + 1, j] 
        data[up_id1, up_id2 + d] = m[i, j + 1] 
        data[up_id1 + d, up_id2 + d] = m[i + 1, j + 1]
    end
    return NonColinMatrix(data)
end

function NonColinMatrix(up::AbstractMatrix{T}, down::AbstractMatrix{T}) where {T}
    @assert size(up) == size(down)
    return NonColinMatrix([up zeros(T, size(up));zeros(T, size(up)) down])
end

Base.Array(c::NonColinMatrix) = copy(c.data)

#TODO Index with atoms
uprange(a::AtomsBase.Atom) = range(a)
    
## Indexing ##
Base.IndexStyle(::AbstractMagneticMatrix) = IndexLinear()
for f in (:view, :getindex)
    @eval function Base.$f(c::ColinMatrix, a1::T, a2::T) where {T <: Atom}
        projrange1 = range(a1)
        projrange2 = range(a2)

        return ColinMatrix($f(c, projrange1, projrange2), $f(c, projrange1, projrange2 .+ blockdim(c)))
    end
    @eval function Base.$f(c::NonColinMatrix, a1::T, a2::T) where {T <: Atom}
        up_range1 = uprange(a1)
        up_range2 = uprange(a2)
        d = blockdim(c)
        dn_range1 = up_range1 .+ d
        dn_range2 = up_range2 .+ d
        return NonColinMatrix([$f(c, up_range1, up_range2) $f(c, up_range1, dn_range2)
                               $f(c, dn_range1, up_range2) $f(c, dn_range1, dn_range2)])
    end

    @eval Base.$f(c::AbstractMagneticMatrix, a1::T) where {T <: Atom} =
        $f(c, a1, a1)

    @eval Base.$f(c::ColinMatrix, a1::T, a2::T, ::Up) where {T <: Atom} =
        $f(c, range(a1), range(a2))
    
    @eval Base.$f(c::NonColinMatrix, a1::T, a2::T, ::Up) where {T <: Atom} =
        $f(c, uprange(a1), uprange(a2))

    @eval Base.$f(c::ColinMatrix, a1::T, a2::T, ::Down) where {T <: Atom} =
        $f(c, range(a1), range(a2) .+ blockdim(c))
        
    @eval Base.$f(c::NonColinMatrix, a1::T, a2::T, ::Down) where {T <: Atom} =
        $f(c, uprange(a1) .+ blockdim(c), uprange(a2) .+ blockdim(c))
        
    @eval Base.$f(c::NonColinMatrix, a1::T, a2::T, ::Up, ::Down) where {T <: Atom} =
        $f(c, uprange(a1), uprange(a2) .+ blockdim(c))

    @eval Base.$f(c::NonColinMatrix, a1::T, a2::T, ::Down, ::Up) where {T <: Atom} =
        $f(c, uprange(a1) .+ blockdim(c), uprange(a2))
        
    @eval Base.$f(c::NonColinMatrix, ::Up, ::Down) =
        (s = size(c,1); $f(c, 1:div(s, 2), div(s, 2)+1:s))

    @eval Base.$f(c::NonColinMatrix, ::Down, ::Up) =
        (s = size(c,1); $f(c, div(s, 2)+1:s, 1:div(s, 2)))
        
    @eval Base.$f(c::AbstractMatrix, a1::T, a2::T, ::Up) where {T<:Atom} =
        $f(c, range(a1), range(a2))

    @eval Base.$f(c::AbstractMatrix, a1::T, a2::T, ::Down) where {T<:Atom} =
        $f(c, range(a1), range(a2) .+ blockdim(c))

    @eval Base.$f(c::AbstractMagneticMatrix, ::Up) =
        (r = 1:blockdim(c); $f(c, r, r))
    @eval Base.$f(c::AbstractMagneticMatrix, ::Up, ::Up) =
        (r = 1:blockdim(c); $f(c, r, r))

    @eval Base.$f(c::ColinMatrix, ::Down) =
        (d = blockdim(c); r = 1:d; $f(c, r, r .+ d))
    @eval Base.$f(c::ColinMatrix, ::Down, ::Down) =
        (d = blockdim(c); r = 1:d; $f(c, r, r .+ d))
        
    @eval Base.$f(c::NonColinMatrix, ::Down) =
        (d = blockdim(c); r = d+1 : 2*d; $f(c, r, r))
    @eval Base.$f(c::NonColinMatrix, ::Down, ::Down) =
        (d = blockdim(c); r = d+1 : 2*d; $f(c, r, r))

    @eval Base.$f(c::AbstractMatrix, ::Up) =
        (r = 1:blockdim(c); $f(c, r, r))
    @eval Base.$f(c::AbstractMatrix, ::Up, ::Up) =
        (r = 1:blockdim(c); $f(c, r, r))

    @eval Base.$f(c::AbstractMatrix, ::Down) =
        (d = blockdim(c); r = d + 1:2 * d; $f(c, r, r))
    @eval Base.$f(c::AbstractMatrix, ::Down, ::Down) =
        (d = blockdim(c); r = d + 1:2 * d; $f(c, r, r))
    
end

for op in (:*, :-, :+, :/)
    @eval @inline Base.$op(c1::ColinMatrix, c2::ColinMatrix) =
        ColinMatrix($op(c1[Up()], c2[Up()]), $op(c1[Down()], c2[Down()]))
    @eval @inline Base.$op(c1::NonColinMatrix, c2::NonColinMatrix) =
        NonColinMatrix($op(c1.data, c2.data))
end

    # BROADCASTING
Base.BroadcastStyle(::Type{T}) where {T<:AbstractMagneticMatrix} =
    Broadcast.ArrayStyle{T}()

Base.ndims(::Type{<:AbstractMagneticMatrix}) =
    2

Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{T}}, ::Type{ElType}) where {T<:AbstractMagneticMatrix,ElType} =
    Base.similar(T, axes(bc))

Base.axes(c::AbstractMagneticMatrix) =
    Base.axes(c.data)

@inline @propagate_inbounds Base.broadcastable(c::AbstractMagneticMatrix) =
    c

@inline @propagate_inbounds Base.unsafe_convert(::Type{Ptr{T}}, c::AbstractMagneticMatrix{T}) where {T} =
    Base.unsafe_convert(Ptr{T}, c.data)

for (elty, cfunc) in zip((:ComplexF32, :ComplexF64), (:cgemm_, :zgemm_))
    @eval @inline function LinearAlgebra.mul!(C::ColinMatrix{$elty}, A::ColinMatrix{$elty}, B::ColinMatrix{$elty})
        dim = blockdim(C)
        dim2 = dim * dim
        ccall((LinearAlgebra.LAPACK.@blasfunc($(cfunc)), liblapack), Cvoid,
                        (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                         Ref{BlasInt}, Ref{$elty}, Ptr{$elty}, Ref{BlasInt},
                         Ptr{$elty}, Ref{BlasInt}, Ref{$elty}, Ptr{$elty},
                         Ref{BlasInt}),
                         'N', 'N', dim, dim,
                         dim, one($elty), A, dim,
                         B, dim, zero($elty), C, dim)
        ccall((LinearAlgebra.LAPACK.@blasfunc($(cfunc)), liblapack), Cvoid,
                        (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                         Ref{BlasInt}, Ref{$elty}, Ptr{$elty}, Ref{BlasInt},
                         Ptr{$elty}, Ref{BlasInt}, Ref{$elty}, Ptr{$elty},
                         Ref{BlasInt}),
                         'N', 'N', dim, dim,
                         dim, one($elty), pointer(A, dim2 + 1), dim,
                         pointer(B, dim2 + 1), dim, zero($elty), pointer(C, dim2 + 1), dim)

        return C
    end
end

@inline function LinearAlgebra.adjoint(c::AbstractMagneticMatrix)
    out = similar(c)
    adjoint!(out, c)
end

@inline @inbounds function LinearAlgebra.adjoint!(out::ColinMatrix, in1::ColinMatrix)
    dim = blockdim(out)
    for i in 1:dim, j in 1:dim
        out[j, i] = in1[i, j]'
        out[j, i + dim] = in1[i, j + dim]'
    end
    return out
end

@inline LinearAlgebra.adjoint!(out::NonColinMatrix, in1::NonColinMatrix) =
    adjoint!(out.data, in1.data)
    
@inline LinearAlgebra.tr(c::AbstractMagneticMatrix) =
    tr(c[Up()]) + tr(c[Down()])

"Vector following the same convention as the in AbstractMagneticMatrix, i.e. first half of the indices contain the up part, second the down part"
struct MagneticVector{T, VT<:AbstractVector{T}} <: AbstractVector{T}
    data::VT
end

for f in (:length, :size, :setindex!, :elsize)
    @eval @inline @propagate_inbounds Base.$f(c::MagneticVector, args...) =
    Base.$f(c.data, args...)
end

up(c::MagneticVector) = view(c.data, 1:div(length(c), 2))
down(c::MagneticVector) = (lc = length(c); view(c.data, div(lc, 2):lc))

# Standard getindex behavior
@inline @propagate_inbounds Base.getindex(c::MagneticVector, args...) =
    getindex(c.data, args...)

for f in (:view, :getindex)
    @eval @inline @propagate_inbounds Base.$f(c::MagneticVector, args::AbstractUnitRange...) =
    Base.$f(c.data, args...)
end

Base.similar(v::MagneticVector) = MagneticVector(similar(v.data))

"Reshuffles standard Wannier90 up-down indices to the ones for the structure of a MagneticVector."
function Base.convert(::Type{MagneticVector}, v::V) where {V <: AbstractVector}
    @assert iseven(length(v)) "Error, dimension of the supplied matrix is odd, i.e. it does not contain both spin components."
    data = similar(v)
    vl = length(v)
    d    =div(vl, 2)
    
    for i in 1:2:vl
        up_id1 = div1(i, 2) 
        data[up_id1] = v[i] 
        data[up_id1 + d] = v[i + 1] 
    end
    return NonColinMatrix(data)
end

## Indexing with atoms and spins ##
# TODO index
for f in (:view, :getindex)
    @eval function Base.$f(c::MagneticVector, a1::T) where {T <: Atom}
        projrange1 = uprange(a1)
        return MagneticVector([$f(c, projrange1); $f(c, projrange1 .+ div(length(c), 2))])
    end
    @eval Base.$f(c::MagneticVector, a1::T, ::Up) where {T <: Atom} =
        $f(c, uprange(a1))
    
    @eval Base.$f(c::MagneticVector, a1::T, ::Down) where {T <: Atom} =
        $f(c, range(a1), uprange(a1) + div(length(c), 2))
        
    @eval Base.$f(c::MagneticVector, ::Up) =
        $f(c, 1:div(length(c), 2))

    @eval Base.$f(c::MagneticVector, ::Down) =
        (lc = length(c); $f(c, div(lc, 2) + 1:lc))
end

for op in (:*, :-, :+, :/)
    @eval @inline Base.$op(c1::MagneticVector, c2::MagneticVector) =
        MagneticVector($op(c1.data, c2.data))
end

import FastLapackInterface: HermitianEigenWs
HermitianEigenWs(c::ColinMatrix) = HermitianEigenWs(c[Up()])
HermitianEigenWs(c::NonColinMatrix) = HermitianEigenWs(c.data)

@inline function eigen!(vals::AbstractVector, vecs::ColinMatrix, c::HermitianEigenWs)
    n = div(length(vals),2)
    n2 = div(length(vecs),2)
    te = eigen!(up(vecs), c)
    copyto!(vals, te.values)
    copyto!(vecs, te.vectors)
    te = eigen!(down(vecs), c)
    copyto!(vals, n+1, te.values, 1, n)
    copyto!(vecs, n2+1, te.vectors, 1, n2)
    return Eigen(vals, vecs)
end

@inline function eigen!(vals::AbstractVector, vecs::NonColinMatrix, c::HermitianEigenWs)
    c.w = vals.data
    te = eigen!(vecs.data, c)
    return Eigen(vals, NonColinMatrix(te.vectors))
end

  
@inline function eigen(vecs::AbstractMatrix{T}, c::HermitianEigenWs{T}) where {T}
    return eigen!(copy(vecs), c)
end

@inline function eigen(vecs::AbstractMagneticMatrix{T}, c::HermitianEigenWs{T}) where {T}
    out  = copy(vecs)
    vals = MagneticVector(similar(out, T <: AbstractFloat ? T : T.parameters[1], size(out, 2)))
    return eigen!(vals, out, c)
end

@inline function eigen(h::AbstractMagneticMatrix)
    return eigen(h, HermitianEigenWs(h))
end

function Base.Matrix(e::Eigen{CT,T,<:ColinMatrix{CT}}) where {CT, T}
    d = size(e.vectors, 1)
    return ColinMatrix([e.vectors[1:d, 1:d] * diagm(0 => e.values[1:d]) * e.vectors[1:d, 1:d]' e.vectors[1:d, d + 1:2d] * diagm(0 => e.values[d + 1:2d]) * e.vectors[1:d, d + 1:2d]'])
end

Base.Array(e::Eigen{CT,T,<:ColinMatrix{CT}}) where {CT, T} = Matrix(e)

# TODO see if this is not duplicated
# KGrids
function uniform_shifted_kgrid(::Type{T}, nkx::Integer, nky::Integer,
                               nkz::Integer, gamma_center = false) where {T}

    t = [Vec3{T}(kx, ky, kz) for kx in 0:nkx-1, ky in 0:nky-1, kz in 0:nkz-1]
    s = Vec3(nkx, nky, nkz)
    t = map(t) do x
        (x .+ 0.5) ./ s .- 0.5
    end
    if gamma_center
        shift = 0.5 .* ((s.+ 1) .% 2)./s
        t = map(t) do x
            x .+ shift
        end
    end
                               
    return reshape(t, nkx * nky * nkz)
end

function uniform_shifted_kgrid(nkx::Integer, nky::Integer, nkz::Integer, gamma_center=false)
    return uniform_shifted_kgrid(Float64, nkx, nky, nkz, gamma_center)
end

using SpecialFunctions, SpecialPolynomials

function ω_grid(ωh, n_ωh, offset = 0.001)
    p = 13
    x, ws= SpecialPolynomials.gauss_nodes_weights(Legendre, n_ωh)
    R= (offset-ωh)/2.0
    R0= (offset+ ωh)/2.0
    y1 = -log(1+p*pi)
    y2 = 0
    y   = (y1 - y2)/2 .* x .- (y2+y1)/2
    phi = (exp.(y) .- 1) ./ p
    path  = R0 .+ R .* exp.(1.0im.*phi)
    return path
end

function integrate_simpson(f, x)
    dx = diff(x)
    N = length(x) 
    result = zero(f(1))
    
    for i = 2:2:length(dx)
        xpx = dx[i] + dx[i-1]
        result += f(i) * (dx[i]^3 + dx[i - 1]^3
                          + 3. * dx[i] * dx[i - 1] * xpx) / (6 * dx[i] * dx[i - 1])
        result += f(i - 1) * (2. * dx[i - 1]^3 - dx[i]^3
                              + 3. * dx[i] * dx[i - 1]^2) / (6 * dx[i - 1] * xpx)
        result += f(i + 1) * (2. * dx[i]^3 - dx[i - 1]^3
                              + 3. * dx[i - 1] * dx[i]^2) / (6 * dx[i] * xpx)
    end

    if length(x) % 2 == 0
        result += f(N) * (2 * dx[end - 1]^2
                  + 3. * dx[end - 2] * dx[end - 1])/ (6 * (dx[end - 2] + dx[end - 1]))
        result += f(N - 1) * (dx[end - 1]^2
                      + 3*dx[end - 1] * dx[end - 2]) / (6 * dx[end - 2])
        result -= f(N - 2) * dx[end - 1]^3/ (6 * dx[end - 2] * (dx[end - 2] + dx[end - 1]))
    end
    return result
end

# struct ExchangeKGrid{T,MT, VT} <: AbstractKGrid{T}
#     hamiltonian_kgrid::HamiltonianKGrid{T,MT, VT}
#     phases::Vector{Complex{T}}
#     D::Matrix{Complex{T}}
# end

# core_kgrid(x::ExchangeKGrid) = core_kgrid(x.hamiltonian_kgrid)
# eigvecs(x::ExchangeKGrid) = eigvecs(x.hamiltonian_kgrid)
# eigvals(x::ExchangeKGrid) = eigvals(x.hamiltonian_kgrid)

# function ExchangeKGrid(hami::TBHamiltonian, kpoints::Vector{Vec3{T}}, R = zero(Vec3{T})) where {T}
#     D = ThreadCache(zeros_block(hami))
#     hami_kpoints = HamiltonianKGrid(hami, kpoints, x -> D .+= x)
#     nk = length(hami_kpoints)
#     Ds = reduce(+, D)
#     return ExchangeKGrid(hami_kpoints, phases(kpoints, R),
#                          Array(Ds[Up()] - Ds[Down()]) / nk)
# end

function calc_greens_functions(ω_grid, kpoints, μ::T) where {T}
    g_caches = [ThreadCache(fill!(similar(eigvecs(kpoints)[1]), zero(Complex{T})))
                for i in 1:3]
    Gs = [fill!(similar(eigvecs(kpoints)[1]), zero(Complex{T})) for i in 1:length(ω_grid)]
    function iGk!(G, ω)
        fill!(G, zero(Complex{T}))
        return integrate_Gk!(G, ω, μ, kpoints, cache.(g_caches))
    end
    p = Progress(length(ω_grid), 1, "Calculating contour G(ω)...")
    Threads.@threads for j in 1:length(ω_grid)
        iGk!(Gs[j], ω_grid[j])
        next!(p)
    end
    return Gs
end

function integrate_Gk!(G::AbstractMatrix, ω, μ, kpoints, caches)
    dim = blockdim(G)
    cache1, cache2, cache3 = caches

    @inbounds for ik in 1:length(kpoints)
        # Fill here needs to be done because cache1 gets reused for the final result too
        fill!(cache1, zero(eltype(cache1)))
        for x in 1:2dim
            cache1[x, x] = 1.0 / (μ + ω - eigvals(kpoints)[ik][x])
        end
        # Basically Hvecs[ik] * 1/(ω - eigvals[ik]) * Hvecs[ik]'

        mul!(cache2,     eigvecs(kpoints)[ik], cache1)
        adjoint!(cache3, eigvecs(kpoints)[ik])
        mul!(cache1,     cache2, cache3)
        t = kpoints.phases[ik]
        tp = t'
        for i in 1:dim, j in 1:dim
            G[i, j]         += cache1[i, j] * t
            G[i+dim, j+dim] += cache1[i+dim, j+dim] * tp
            G[i+dim, j]      = cache1[i+dim, j]
            G[i, j+dim]      = cache1[i, j+dim]
        end
    end
    return G ./= length(kpoints)
end

function integrate_Gk!(G::ColinMatrix, ω, μ, kpoints, caches)
    dim = size(G, 1)
    cache1, cache2, cache3 = caches

    @inbounds for ik in 1:length(kpoints)
        # Fill here needs to be done because cache1 gets reused for the final result too
        fill!(cache1, zero(eltype(cache1)))
        for x in 1:dim
            cache1[x, x]     = 1.0 / (μ + ω - kpoints.hamiltonian_kgrid.eigvals[ik][x])
            cache1[x, x+dim] = 1.0 / (μ + ω - kpoints.hamiltonian_kgrid.eigvals[ik][x+dim])
        end
        # Basically Hvecs[ik] * 1/(ω - eigvals[ik]) * Hvecs[ik]'

        mul!(cache2, kpoints.hamiltonian_kgrid.eigvecs[ik], cache1)
        adjoint!(cache3, kpoints.hamiltonian_kgrid.eigvecs[ik])
        mul!(cache1, cache2, cache3)
        t = kpoints.phases[ik]
        tp = t'
        for i in 1:dim, j in 1:dim
            G[i, j]     += cache1[i, j] * t
            G[i, j+dim] += cache1[i, j+dim] * tp
        end
    end
    return G ./= length(kpoints)
end

function integrate_Gk!(G_forward::ThreadCache, G_backward::ThreadCache, ω, μ, Hvecs, Hvals,
                       R, kgrid, caches)
    dim = size(G_forward, 1)
    cache1, cache2, cache3 = caches

    @inbounds for ik in 1:length(kgrid)
        # Fill here needs to be done because cache1 gets reused for the final result too
        fill!(cache1, zero(eltype(cache1)))
        for x in 1:dim
            cache1[x, x] = 1.0 / (μ + ω - Hvals[ik][x])
        end
        # Basically Hvecs[ik] * 1/(ω - eigvals[ik]) * Hvecs[ik]'
        mul!(cache2, Hvecs[ik], cache1)
        adjoint!(cache3, Hvecs[ik])
        mul!(cache1, cache2, cache3)
        t = exp(2im * π * dot(R, kgrid[ik]))
        
        G_forward  .+= cache1 .* t
        G_backward .+= cache1 .* t'
    end
    G_forward.caches ./= length(kgrid)
    return G_backward.caches ./= length(kgrid)
end

abstract type Exchange{T<:AbstractFloat} end
Base.eltype(::Exchange{T}) where {T} = T
Base.eltype(::Type{Exchange{T}}) where {T} = T

function (::Type{E})(at1::Atom, at2::Atom; site_diagonal::Bool = false) where {E<:Exchange}
    l1 = length(uprange(at1))
    l2 = length(uprange(at2))
    return site_diagonal ? E(zeros(Float64, l1, l2), at1, at2) :
           E(zeros(Float64, l1, l1), at1, at2)
end

"""
    Exchange2ndOrder{T <: AbstractFloat}

This holds the exhanges between different orbitals and calculated sites.
Projections and atom datablocks are to be found in the corresponding wannier input file.
It turns out the ordering is first projections, then atom order in the atoms datablock.
"""
mutable struct Exchange2ndOrder{T<:AbstractFloat} <: Exchange{T}
    J::Matrix{T}
    atom1::Atom
    atom2::Atom
end

function Base.show(io::IO, e::Exchange)
    print(io, "atom1:")
    println(io, "name: $(e.atom1.atomic_symbol), pos: $(e.atom1.position)")
    print(io, " atom2:")
    println(io, "name: $(e.atom2.atomic_symbol), pos: $(e.atom2.position)")

    return print(io, " J: $(tr(e.J))")
end

"""
    Exchange4thOrder{T <: AbstractFloat}

This holds the exhanges between different orbitals and calculated sites.
Projections and atom datablocks are to be found in the corresponding wannier input file.
It turns out the ordering is first projections, then atom order in the atoms datablock.
"""
mutable struct Exchange4thOrder{T<:AbstractFloat} <: Exchange{T}
    J::Matrix{T}
    atom1::Atom
    atom2::Atom
end

"""
    calc_exchanges(hamiltonian::TBHamiltonian, atoms::Vector{<:Atom}, fermi, exchange_type; kwargs...)

Calculates the magnetic exchange parameters between the `atoms`. `exchange_type` can be [`Exchange2ndOrder`](@ref) or [`Exchange4thOrder`](@ref). The `kwargs` control various numerical parameters for the calculation:
- `nk = (10,10,10)`: the amount of _k_-points to be used for the uniform interpolation grid.
- `R = (0,0,0)`: the unit cell index to which the exchange parameters are calculated.
- `ωh = -30.0`: the lower bound of the energy integration
- `ωv = 0.15`: the height of the contour in complex space to integrate the Green's functions
- `n_ωh = 3000`: number of integration points along the horizontal contour direction
- `n_ωv = 500`: number of integration points along the vertical contour direction
- `site_diagonal = false`: if `true` the hamiltonians and `Δ` will diagonalized on-site and the
returned exchange matrices hold the exchanges between well-defined orbitals. If this is not done,
the exchange matrices entries don't mean anything on themselves and a trace should be performed to
find the exchange between the spins on sites `i` and `j`.
""" 
"""
    calc_exchanges(hamiltonian::TBHamiltonian, atoms::Vector{<:Atom}, fermi, exchange_type; kwargs...)

Calculates the magnetic exchange parameters between the `atoms`. `exchange_type` can be [`Exchange2ndOrder`](@ref) or [`Exchange4thOrder`](@ref). The `kwargs` control various numerical parameters for the calculation:
- `nk = (10,10,10)`: the amount of _k_-points to be used for the uniform interpolation grid.
- `R = (0,0,0)`: the unit cell index to which the exchange parameters are calculated.
- `ωh = -30.0`: the lower bound of the energy integration
- `ωv = 0.15`: the height of the contour in complex space to integrate the Green's functions
- `n_ωh = 3000`: number of integration points along the horizontal contour direction
- `n_ωv = 500`: number of integration points along the vertical contour direction
- `site_diagonal = false`: if `true` the hamiltonians and `Δ` will diagonalized on-site and the
returned exchange matrices hold the exchanges between well-defined orbitals. If this is not done,
the exchange matrices entries don't mean anything on themselves and a trace should be performed to
find the exchange between the spins on sites `i` and `j`.
""" 
function calc_exchanges(hami, structure::DFControl.Structures.Structure, atsyms::Vector{Symbol}, fermi::T, ::Type{E} = Exchange2ndOrder;
                        nk::NTuple{3,Int} = (10, 10, 10),
                        R = Vec3(0, 0, 0),
                        ωh::T = T(-30.0), # starting energy
                        n_ωh::Int = 100,
                        emax::T = T(0.001)) where {T<:AbstractFloat,E<:Exchange}
    R_     = Vec3(R...)
    μ      = fermi
    ω_grid = ω_grid(ωh, n_ωh, emax)

    exchanges = E{T}[]
    for at1 in Iterators.flatten(map(x->structure[element(x)], atsyms))
        for at2 in Iterators.flatten(map(x->structure[element(x)], atsyms))
            at2_ = deepcopy(at2)
            at2_.position_cart =  at2_.position_cart .+ structure.cell * Vec3(R...)
            push!(exchanges, E(at1, at2_))
        end
    end
    kpoints = ExchangeKGrid(hami, uniform_shifted_kgrid(nk...), R_)

    calc_exchanges!(exchanges, μ, ω_grid, kpoints, kpoints.D)

    return exchanges
end

function calc_exchanges!(exchanges::Vector{<:Exchange{T}},
                         μ::T,
                         ω_grid::AbstractArray{Complex{T}},
                         kpoints,
                         D::Matrix{Complex{T}}) where {T<:AbstractFloat}
    dim = size(kpoints.hamiltonian_kgrid.eigvecs[1])
    d2 = div(dim[1], 2)
    J_caches = [ThreadCache(zeros(T, size(e.J))) for e in exchanges]
    Gs = calc_greens_functions(ω_grid, kpoints, μ)
    for j in 1:length(exchanges)
        J_caches[j] .+= imag.(integrate_simpson(i -> Jω(exchanges[j], D, Gs[i]), ω_grid))
    end
    for (eid, exch) in enumerate(exchanges)
        exch.J = -1e3 / 4π * reduce(+, J_caches[eid])
    end
end

spin_sign(D) = -sign(real(tr(D))) # up = +1, down = -1. If D_upup > D_dndn, onsite spin will be down and the tr(D) will be positive. Thus explaining the - in front of this.
spin_sign(D::Vector) = sign(real(sum(D))) # up = +1, down = -1

@inline function Jω(exch, D, G)
    if size(D, 1) < size(G, 1)
        ra1 = uprange(exch.atom1)
        ra2 = uprange(exch.atom2)
    else
        ra1 = range(exch.atom1)
        ra2 = range(exch.atom2)
    end
    D_site1    = view(D, ra1, ra1)
    D_site2    = view(D, ra2, ra2)
    s1         = spin_sign(D_site1)
    s2         = spin_sign(D_site2)
    t          = zeros(ComplexF64, size(exch.J))
    G_forward  = view(G, exch.atom1, exch.atom2, Up())
    G_backward = view(G, exch.atom2, exch.atom1, Down())
    for j in 1:size(t, 2), i in 1:size(t, 1)
        t[i, j] = s1 * s2 *
                  D_site1[i,i] *
                       G_forward[i, j] *
                       D_site2[j, j] *
                       G_backward[j, i]
    end
    return t
end

mutable struct AnisotropicExchange2ndOrder{T<:AbstractFloat} <: Exchange{T}
    J::Matrix{Matrix{T}}
    atom1::Atom
    atom2::Atom
end

function AnisotropicExchange2ndOrder(at1::Atom, at2::Atom)
    return AnisotropicExchange2ndOrder{Float64}([zeros(length(range(at1)),
                                                       length(range(at1)))
                                                 for i in 1:3, j in 1:3], at1, at2)
end

function calc_anisotropic_exchanges(hami, atoms, fermi::T;
                                    nk::NTuple{3,Int} = (10, 10, 10),
                                    R = Vec3(0, 0, 0),
                                    ωh::T = T(-30.0), # starting energy
                                    ωv::T = T(0.1), # height of vertical contour
                                    n_ωh::Int = 3000,
                                    n_ωv::Int = 500,
                                    temp::T = T(0.01)) where {T<:AbstractFloat}
    μ         = fermi
    k_grid    = uniform_shifted_kgrid(nk...)
    ω_grid    = setup_ω_grid(ωh, ωv, n_ωh, n_ωv)
    exchanges = setup_anisotropic_exchanges(atoms)

    Hvecs, Hvals, D = DHvecvals(hami, k_grid, atoms)

    calc_anisotropic_exchanges!(exchanges, μ, R, k_grid, ω_grid, Hvecs, Hvals, D)
    return exchanges
end

function calc_anisotropic_exchanges!(exchanges::Vector{AnisotropicExchange2ndOrder{T}},
                                     μ::T,
                                     R::Vec3,
                                     k_grid::AbstractArray{Vec3{T}},
                                     ω_grid::AbstractArray{Complex{T}},
                                     Hvecs::Vector{Matrix{Complex{T}}},
                                     Hvals::Vector{Vector{T}},
                                     D::Vector{Vector{Matrix{Complex{T}}}}) where {T<:AbstractFloat}
    dim = size(Hvecs[1])
    J_caches = [ThreadCache([zeros(T, size(e.J[i, j])) for i in 1:3, j in 1:3])
                for e in exchanges]
    g_caches = [ThreadCache(zeros(Complex{T}, dim)) for i in 1:3]
    G_forward, G_backward = [ThreadCache(zeros(Complex{T}, dim)) for i in 1:2]

    function iGk!(ω)
        fill!(G_forward, zero(Complex{T}))
        fill!(G_backward, zero(Complex{T}))
        return integrate_Gk!(G_forward, G_backward, ω, μ, Hvecs, Hvals, R, k_grid, g_caches)
    end

    for j in 1:length(ω_grid[1:end-1])
        ω  = ω_grid[j]
        dω = ω_grid[j+1] - ω
        iGk!(ω)
        # The two kind of ranges are needed because we calculate D only for the projections we care about 
        # whereas G is calculated from the full Hamiltonian, the is needed. 
        for (eid, exch) in enumerate(exchanges)
            rm = range(exch.atom1)
            rn = range(exch.atom2)
            for i in 1:3, j in 1:3 # x,y,z 
                J_caches[eid][i, j] .+= imag.((view(D[eid][i], 1:length(rm), 1:length(rm)) *
                                               view(G_forward, rm, rn) *
                                               view(D[eid][j], 1:length(rn), 1:length(rn)) *
                                               view(G_backward, rn, rm)) .* dω)
            end
        end
    end

    for (eid, exch) in enumerate(exchanges)
        exch.J = 1e3 / 2π * reduce(+, J_caches[eid])
    end
end

function setup_anisotropic_exchanges(atoms::Vector{Atom})
    exchanges = AnisotropicExchange2ndOrder{Float64}[]
    for (i, at1) in enumerate(atoms), at2 in atoms[i:end]
        push!(exchanges, AnisotropicExchange2ndOrder(at1, at2))
    end
    return exchanges
end

@doc raw"""
    DHvecvals(hami::TBHamiltonian{T, Matrix{T}}, k_grid::Vector{Vec3{T}}, atoms::Atom{T}) where T <: AbstractFloat

Calculates $D(k) = [H(k), J]$, $P(k)$ and $L(k)$ where $H(k) = P(k) L(k) P^{-1}(k)$.
`hami` should be the full Hamiltonian containing both spin-diagonal and off-diagonal blocks.
"""
function DHvecvals(hami, k_grid::AbstractArray{Vec3{T}},
                   atoms::Vector{Atom}) where {T<:AbstractFloat}
    nk = length(k_grid)
    Hvecs = [zeros_block(hami) for i in 1:nk]
    Hvals = [Vector{T}(undef, blocksize(hami, 1)) for i in 1:nk]
    δH_onsite = ThreadCache([[zeros(Complex{T}, 2length(range(at)), 2length(range(at)))
                              for i in 1:3] for at in atoms])
    calc_caches = [HermitianEigenWs(block(hami[1])) for i in 1:nthreads()]
    for i in 1:nk
        # for i=1:nk
        tid = threadid()
        # Hvecs[i] is used as a temporary cache to store H(k) in. Since we
        # don't need H(k) but only Hvecs etc, this is ok.
        Hk!(Hvecs[i], hami, k_grid[i])

        for (δh, at) in zip(δH_onsite, atoms)
            rat = range(at)
            lr  = length(rat)
            δh .+= commutator.((view(Hvecs[i], rat, rat),), at[:operator_block].J) # in reality this should be just range(at)
            # δh .+= commutator.(([Hvecs[i][rat, rat] zeros(Complex{T},lr, lr); zeros(Complex{T}, lr, lr) Hvecs[i][div(blocksize(hami, 1), 2) .+ rat, div(blocksize(hami, 1), 2) .+ rat]],), at[:operator_block].J) #in reality this should be just range(at)
        end
        eigen!(Hvals[i], Hvecs[i], calc_caches[tid])
    end
    return Hvecs, Hvals, reduce(+, δH_onsite) ./ nk
end

commutator(A1, A2) = A1 * A2 - A2 * A1

function totocc(Hvals, fermi::T, temp::T) where {T}
    totocc = zero(Complex{T})
    for i in 1:length(Hvals)
        totocc += reduce(+, 1 ./ (exp.((Hvals[i] .- fermi) ./ temp) .+ 1))
    end
    return totocc / length(Hvals)
end

"""
    TBBlock

Building block for [`TBOperator`](@ref). It holds the matrix elements of the operator between
central and a shifted unit cell. Upon construction, the wigner-seitz shifts are taken into
account to create the correct matrix elements between the Wannierfunctions, stored in
`tb_block`. The `block` field is basically `tb_block` but with each element divided by
the amount of Wigner-Seitz degeneracies and shifts which speeds up later _k_-point interpolation.
"""
struct TBBlock{T<:AbstractFloat,LT,M<:AbstractMatrix{Complex{T}}}
    R_cryst  :: Vec3{Int}
    R_cart   :: Vec3{LT}
    block    :: M
    tb_block :: M
end

"""
    TBOperator

Alias for a `Vector` of [`TBBlocks`](@ref TBBlock). Indexing with `NTuple{3,Int}` or `Vec3`
is supported which allows for easily retrieving the [`TBBlock`](@ref) that corresponds
to the shifted unit cell.
Aliases: `TBHamiltonian`, `TBSpin`
"""
const TBOperator{T,LT,M} = Vector{TBBlock{T,LT,M}}
const TBHamiltonian = TBOperator
const TBSpin = TBOperator
block(x::TBBlock) = x.block
for f in (:getindex, :size, :similar)
    @eval Base.$f(h::TBBlock, args...) = $f(block(h), args...)
end

LinearAlgebra.eigen(h::TBBlock) = eigen(block(h))

Base.getindex(h::TBOperator, R::Vec3{Int}) = getfirst(x -> x.R_cryst == R, h)

zeros_block(h::TBOperator) = zeros(block(h[1]))

similar_block(h::TBOperator) = similar(block(h[1]))

blocksize(h::TBOperator, args...) = size(block(h[1]), args...)

n_wann(h::TBOperator) = blocksize(h, 1)
n_rvecs(h::TBOperator) = filter(x -> !iszero(sum(x.tb_block)), h)

for op in (:+, :-, :*, :/)
    @eval function Base.$op(t::TBBlock{T}, v::T) where {T}
        return TBBlock(t.R_cart, t.R_cryst, $op(block(t), v), $op(t.tb_block, v))
    end
    @eval function Base.$op(v::T, t::TBBlock{T}) where {T}
        return TBBlock(t.R_cart, t.R_cryst, $op(v, block(t)), $op(v, t.tb_block))
    end
    @eval function Base.$op(t::TBBlock{T,M}, v::M) where {T,M}
        return TBBlock(t.R_cart, t.R_cryst, $op(block(t), v), $op(t.tb_block, v))
    end
    @eval function Base.$op(v::M, t::TBBlock{T,M}) where {T,M}
        return TBBlock(t.R_cart, t.R_cryst, $op(v, block(t)), $op(v, t.tb_block))
    end
    @eval function Base.$op(t::TBBlock{T,M}, v::TBBlock{T,M}) where {T,M}
        return TBBlock(t.R_cart, t.R_cryst, $op(block(t), block(v)),
                       $op(t.tb_block, v.tb_block))
    end
end

Base.isapprox(b1::TBBlock, b2::TBBlock, args...;kwargs...) = all(x -> isapprox(x[1], x[2], args...;kwargs...), zip(b1.block, b2.block))


function HR_ws(Hq, kpts, R_cryst_ws, n_wann)
    HR = [zeros(eltype(Hq[1]), n_wann, n_wann) for i = 1:length(R_cryst_ws)]

    fourier(kpts, R_cryst_ws) do iR, ik, phase
        @inbounds HR[iR] .+= phase .* Hq[ik]
    end

    HR ./= length(kpts)
end

function TBHamiltonian(model::Model, rvectors::RVectors; kwargs...)
    kpts = model.kpoints
    Hq = get_Hk(model.E, model.U)
    R_cryst = rvectors.R

    HR = HR_ws(Hq, kpts, R_cryst, model.n_wann)
    
    return map(enumerate(HR)) do (iR, H)
        rcryst = R_cryst[iR]
        rcart  = model.lattice * rcryst
        return TBBlock(rcryst, rcart, H ./ rvectors.N[iR], H)
    end
end

function mdrs_v1tov2(Oᴿ::Vector{MT}, rvectors::RVectorsMDRS) where {MT <: AbstractMatrix}
    
    R_cryst = rvectors.R̃vectors.R
    R_cryst_ws = rvectors.Rvectors.R

    n_r̃vecs = rvectors.n_r̃vecs
    n_wann = size(Oᴿ[1], 1)

    # This is for when there's no "real" tb block for a given mdrs R
    zeros_mat = zeros(Oᴿ[1])
    return map(enumerate(R_cryst)) do (iR, rcryst)
        rcart  = rvectors.Rvectors.lattice * rcryst

        H = zeros(Oᴿ[1])
        @inbounds for (ir, m, n, it) in rvectors.R̃_RT[iR]
            fac = rvectors.Nᵀ[ir][m, n]
            fac *= rvectors.N[ir]

            H[m, n] += Oᴿ[ir][m, n] / fac
        end
        rid = findfirst(isequal(rcryst), R_cryst_ws)
        return TBBlock(rcryst, rcart, H, rid === nothing ? zeros_mat : Oᴿ[rid])
    end
end

function TBHamiltonian(model::Model, rvectors::RVectorsMDRS; kwargs...)

    kpts = model.kpoints

    Hq = get_Hk(model.E, model.U)
    
    HR = HR_ws(Hq, kpts, rvectors.Rvectors.R, model.n_wann)

    return mdrs_v1tov2(HR, rvectors)
end




# TODO: see DFWannier for wigner seitz stuff
# function TBHamiltonian(model::Model; kwargs...)
#     R = wigner_seitz_R(model; kwargs...)

#     kpts = model.kpoints

#     Hq = get_Hk(model.E, model.U)
    
#     HR = [zeros(ComplexF64, model.n_wann, model.n_wann) for i = 1:length(R.cryst)]
    
#     fourier(kpts, R.cryst) do iR, ik, phase
#         HR[iR] .+= phase .* Hq[ik]
#     end

#     HR ./= size(kpts, 2)

#     return generate_TBBlocks(HR, model.lattice, R; kwargs...)
# end

function compute_MU_UtMU!(
        MU,
        UtMU,
        bvectors,
        M,
        U::AbstractVector{<:AbstractMatrix{<:Number}},
    )
    kpb_k = bvectors.kpb_k
    n_bvecs = length(kpb_k[1])

    @inbounds for ik in 1:length(U)
        Ut = U[ik]'
        for ib in 1:n_bvecs
            MUkb = MU[ik][ib]
            ikpb = kpb_k[ik][ib]
            mul!(MUkb, M[ik][ib], U[ikpb])
            mul!(UtMU[ik][ib], Ut, MUkb)
        end
    end
    return MU, UtMU
end

function compute_MU_UtMU!(MU, UtMU, bvectors, M, U::AbstractArray{<:Number, 3})
    kpb_k = bvectors.kpb_k
    n_bvecs = length(kpb_k[1])

    @inbounds for ik in axes(U, 3)
        Ut = view(U, :, :, ik)'
        for ib in 1:n_bvecs
            MUkb = MU[ik][ib]
            ikpb = kpb_k[ik][ib]
            Ukpb = view(U, :, :, ikpb)
            mul!(MUkb, M[ik][ib], Ukpb)
            mul!(UtMU[ik][ib], Ut, MUkb)
        end
    end
    return MU, UtMU
end

function compute_MU_UtMU!(
        MU::AbstractArray{<:Number, 4},
        UtMU::AbstractArray{<:Number, 4},
        bvectors,
        M,
        U::AbstractVector{<:AbstractMatrix{<:Number}},
    )
    kpb_k = bvectors.kpb_k
    n_bvecs = length(kpb_k[1])

    @inbounds for ik in 1:length(U)
        Ut = U[ik]'
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            MUkb = view(MU, :, :, ib, ik)
            Nkb = view(UtMU, :, :, ib, ik)
            mul!(MUkb, M[ik][ib], U[ikpb])
            mul!(Nkb, Ut, MUkb)
        end
    end
    return MU, UtMU
end

function compute_MU_UtMU!(
        MU::AbstractArray{<:Number, 4},
        UtMU::AbstractArray{<:Number, 4},
        bvectors,
        M,
        U::AbstractArray{<:Number, 3},
    )
    kpb_k = bvectors.kpb_k
    n_bvecs = length(kpb_k[1])

    @inbounds for ik in axes(U, 3)
        Ut = view(U, :, :, ik)'
        for ib in 1:n_bvecs
            ikpb = kpb_k[ik][ib]
            Ukpb = view(U, :, :, ikpb)
            MUkb = view(MU, :, :, ib, ik)
            Nkb = view(UtMU, :, :, ib, ik)
            mul!(MUkb, M[ik][ib], Ukpb)
            mul!(Nkb, Ut, MUkb)
        end
    end
    return MU, UtMU
end

function compute_MU_UtMU!(cache, bvectors, M, U)
    return compute_MU_UtMU!(getfield(cache, :MU), getfield(cache, :UtMU), bvectors, M, U)
end

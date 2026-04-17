function compute_MU_UtMU!(
        MU::AbstractArray{<:Number, 4},
        UtMU::AbstractArray{<:Number, 4},
        bvectors,
        M::AbstractArray{<:Number, 4},
        U::AbstractArray{<:Number, 3},
    )
    kpb_k = bvectors.kpb_k
    n_bvecs = size(kpb_k, 1)

    @inbounds for ik in axes(U, 3)
        Ut = view(U, :, :, ik)'
        for ib in 1:n_bvecs
            ikpb = kpb_k[ib, ik]
            Ukpb = view(U, :, :, ikpb)
            MUkb = view(MU, :, :, ib, ik)
            Nkb = view(UtMU, :, :, ib, ik)
            mul!(MUkb, view(M, :, :, ib, ik), Ukpb)
            mul!(Nkb, Ut, MUkb)
        end
    end
    return MU, UtMU
end

function compute_MU_UtMU!(workspace, bvectors, M, U)
    return compute_MU_UtMU!(workspace.MU, workspace.UtMU, bvectors, M, U)
end

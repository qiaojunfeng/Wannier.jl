@testitem "position_transl_inv_full" begin
    using LinearAlgebra
    using Wannier.Datasets
    model = read_w90_with_chk(
        dataset"Si2_coarse/Si2", dataset"Si2_coarse/outputs/Si2.chk"
    )
    Rspace = Wannier.generate_Rspace(model)

    # incompatible with the Im-log band-diagonal, as in wannier90
    @test_throws ErrorException TBPosition(
        Rspace, model; transl_inv_full = true, imlog_diag = true
    )
    @test_throws ErrorException TBPosition(
        Rspace, model; transl_inv_full = true, force_hermiticity = true
    )

    # `transl_inv_full = true` alone is enough: the other two are implied
    A = TBPosition(Rspace, model; transl_inv_full = true)
    # the R-vectors are untouched, only the operator changes
    A0 = TBPosition(Rspace, model; force_hermiticity = false)
    @test A.Rspace.Rvectors == A0.Rspace.Rvectors

    # the R = 0 band-diagonal is pinned to the Wannier centres
    iR0 = findfirst(iszero, A.Rspace.Rvectors)
    centers = omega(model).r
    nwann = n_wannier(model)
    @test all(A.operator[iR0][i, i] ≈ centers[i] for i in 1:nwann)
end

@testitem "position_transl_inv_full translational invariance" begin
    using LinearAlgebra
    using Wannier.Datasets
    model = read_w90_with_chk(
        dataset"Si2_coarse/Si2", dataset"Si2_coarse/outputs/Si2.chk"
    )
    Rspace = Wannier.generate_Rspace(model)

    # A rigid shift of every Wannier function by `t` sends
    # M^{k,b} -> exp(-i b⋅t) M^{k,b}, under which <0i|r|Rj> must gain exactly
    # t δ_ij δ_{R,0} and nothing else. This is the property `transl_inv_full`
    # exists for.
    t = [0.37, -0.21, 0.13]
    kstencil = model.kstencil
    recip_lattice = Wannier.reciprocal_lattice(kstencil)
    shifted = map(1:n_kpoints(kstencil)) do ik
        map(1:n_bvectors(kstencil)) do ib
            b = recip_lattice * (
                kstencil.kpoints[kstencil.kpb_k[ik][ib]] +
                    kstencil.kpb_G[ik][ib] - kstencil.kpoints[ik]
            )
            cis(-(b ⋅ t)) .* model.overlaps[ik][ib]
        end
    end
    # the transformation really is a rigid shift of the centres
    @test omega(kstencil, shifted, model.gauges).r ≈ omega(model).r .+ Ref(t)

    # `Model(model, kstencil, overlaps)` reuses the lattice and atomic data;
    # the positional form would put `entangled_bands` into the `frozen_bands`
    # slot
    model_shifted = Wannier.Model(model, kstencil, shifted)

    function shift_error(; kwargs...)
        A = TBPosition(Rspace, model; force_hermiticity = false, kwargs...)
        B = TBPosition(Rspace, model_shifted; force_hermiticity = false, kwargs...)
        iR0 = findfirst(iszero, A.Rspace.Rvectors)
        worst = 0.0
        for iR in eachindex(A.operator)
            for i in axes(A.operator[iR], 1), j in axes(A.operator[iR], 2)
                expected = (iR == iR0 && i == j) ? t : zero(t)
                worst = max(
                    worst,
                    maximum(abs.(B.operator[iR][i, j] - A.operator[iR][i, j] - expected)),
                )
            end
        end
        return worst
    end

    @test shift_error(transl_inv_full = true) < 1.0e-12
    # the naive discretisation is not translationally invariant
    @test shift_error(transl_inv_full = false) > 1.0e-3
end

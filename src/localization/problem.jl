abstract type AbstractLocalizationTerm end

export VarianceTerm, CenterConstraintTerm

struct VarianceTerm <: AbstractLocalizationTerm end

struct CenterConstraintTerm{T <: Real} <: AbstractLocalizationTerm
    r0::Vector{Vec3{T}}
    λ::T
end

struct LocalizationProblem{TT, M, O}
    terms::TT
    model::M
    parameterization::Symbol
    solver_options::O
end

function LocalizationProblem(
        terms::Tuple,
        model,
        parameterization::Symbol;
        kwargs...,
    )
    return LocalizationProblem(terms, model, parameterization, (; kwargs...))
end

function build_fg!(problem::LocalizationProblem)
    if problem.parameterization == :maxloc
        return _build_fg_maxloc(problem)
    elseif problem.parameterization == :disentangle
        return _build_fg_disentangle(problem)
    elseif problem.parameterization == :mag_disentangle
        return _build_fg_mag_disentangle(problem)
    elseif problem.parameterization == :mag_disentangle_center
        return _build_fg_mag_disentangle_center(problem)
    else
        error("Unsupported localization parameterization: $(problem.parameterization)")
    end
end

@inline _is_variance_term(::VarianceTerm) = true
@inline _is_variance_term(::AbstractLocalizationTerm) = false
@inline _is_center_term(::CenterConstraintTerm) = true
@inline _is_center_term(::AbstractLocalizationTerm) = false

function _has_variance_term(terms)
    return any(_is_variance_term, terms)
end

function _has_center_term(terms)
    return any(_is_center_term, terms)
end

# Keep the common constrained-maxloc pair on a dedicated path.
#
# Although summing terms compositionally is mathematically equivalent,
# finite-precision differences in value/gradient accumulation can shift
# short LBFGS trajectories (e.g. fixed small `max_iter` tests). This helper
# lets us reuse the legacy center-penalty evaluation ordering for
# `(VarianceTerm(), CenterConstraintTerm(...))` to preserve historical
# numerical behavior.
function _is_variance_plus_center_terms(terms)
    if length(terms) != 2
        return false
    end
    return _has_variance_term(terms) && _has_center_term(terms)
end

function _find_center_term(terms)
    for term in terms
        if term isa CenterConstraintTerm
            return term
        end
    end
    return nothing
end

function _accumulate_spread_terms_value!(terms, Ωbase::Spread)
    # Preserve legacy objective evaluation for the dominant two-term case.
    if _is_variance_plus_center_terms(terms)
        center_term = _find_center_term(terms)
        return omega_center(Ωbase; r₀ = center_term.r0, λ = center_term.λ).Ωt
    end

    Ω = zero(Ωbase.Ω)
    @inbounds for term in terms
        if term isa VarianceTerm
            Ω += Ωbase.Ω
        elseif term isa CenterConstraintTerm
            Ω += omega_center(Ωbase; r₀ = term.r0, λ = term.λ).Ωc
        end
    end
    return Ω
end

function _accumulate_spread_terms_grad!(
        Gacc,
        terms,
        cache,
        Gbase,
        Gtmp,
        kstencil,
        overlaps,
    )
    fill!(Gacc, 0)

    # Preserve legacy gradient path for numerical compatibility with
    # short-iteration localization trajectories.
    if _is_variance_plus_center_terms(terms)
        center_term = _find_center_term(terms)
        cache.G = Gacc
        omega_grad!(center_penalty(center_term.r0, center_term.λ), cache, kstencil, overlaps)
        return Gacc
    end

    need_base = _has_variance_term(terms) || _has_center_term(terms)
    if need_base
        cache.G = Gbase
        omega_grad!(cache, kstencil, overlaps)
    end

    @inbounds for term in terms
        if term isa VarianceTerm
            Gacc .+= Gbase
        elseif term isa CenterConstraintTerm
            cache.G = Gtmp
            omega_grad!(center_penalty(term.r0, term.λ), cache, kstencil, overlaps)
            Gacc .+= (Gtmp .- Gbase)
        end
    end

    cache.G = Gacc
    return Gacc
end

function _build_fg_maxloc(problem::LocalizationProblem)
    terms = problem.terms
    model = problem.model
    cache = Cache(model)
    Gbase = similar(cache.G)
    Gtmp = similar(cache.G)

    function fg!(F, G, U)
        compute_MU_UtMU!(cache, model.kstencil, model.overlaps, U)

        Ωbase = nothing
        if F !== nothing || _has_center_term(terms)
            Ωbase = omega!(cache, model.kstencil, model.overlaps)
        end

        if G !== nothing
            cache.G = G
            _accumulate_spread_terms_grad!(
                G,
                terms,
                cache,
                Gbase,
                Gtmp,
                model.kstencil,
                model.overlaps,
            )
        end
        if F !== nothing
            return _accumulate_spread_terms_value!(terms, Ωbase)
        end
    end

    return fg!
end

function _build_fg_disentangle(problem::LocalizationProblem)
    terms = problem.terms
    model = problem.model
    cache = Cache(model)
    Gbase = similar(cache.G)
    Gtmp = similar(cache.G)
    Gacc = similar(cache.G)

    function fg!(Ω, G, XY)
        X, Y = XY_to_X_Y!(cache.X, cache.Y, XY)
        U = X_Y_to_U!(cache.U, X, Y)
        compute_MU_UtMU!(cache, model.kstencil, model.overlaps, U)

        Ωbase = nothing
        if Ω !== nothing || _has_center_term(terms)
            Ωbase = omega!(cache, model.kstencil, model.overlaps)
        end

        if G !== nothing
            G_ = _accumulate_spread_terms_grad!(
                Gacc,
                terms,
                cache,
                Gbase,
                Gtmp,
                model.kstencil,
                model.overlaps,
            )
            GX, GY = GU_to_GX_GY(G_, X, Y, model.frozen_bands)

            n = n_wannier(model)^2

            @inbounds for ik in 1:n_kpoints(model)
                gxk = view(GX, :, :, ik)
                gyk = view(GY, :, :, ik)
                for i in eachindex(gxk)
                    G[i, ik] = gxk[i]
                end
                for i in eachindex(gyk)
                    G[n + i, ik] = gyk[i]
                end
            end
        end
        if Ω !== nothing
            return _accumulate_spread_terms_value!(terms, Ωbase)
        end
    end

    return fg!
end

function _build_fg_mag_disentangle(problem::LocalizationProblem)
    model = problem.model
    λ = get(problem.solver_options, :lambda, 1.0)

    nb = n_bands(model.up)
    nw = n_wannier(model.up)
    nk = n_kpoints(model.up)
    n_inner = nb * nw + nw^2

    function f(XY)
        XY = reshape(XY, (2 * n_inner, nk))
        XYup = @view XY[1:n_inner, :]
        XYdn = @view XY[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, nb, nw)
        Xdn, Ydn = XY_to_X_Y(XYdn, nb, nw)
        Ωup = omega(model.up.kstencil, model.up.overlaps, Xup, Yup).Ω
        Ωdn = omega(model.dn.kstencil, model.dn.overlaps, Xdn, Ydn).Ω
        if λ == 0
            Ωupdn = 0
        else
            Ωupdn = omega_updn(model, X_Y_to_U(Xup, Yup), X_Y_to_U(Xdn, Ydn))
        end
        return Ωup + Ωdn + λ * Ωupdn
    end

    function g!(G, XY)
        XY = reshape(XY, (2 * n_inner, nk))
        XYup = @view XY[1:n_inner, :]
        XYdn = @view XY[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, nb, nw)
        Xdn, Ydn = XY_to_X_Y(XYdn, nb, nw)
        GXup, GYup = omega_grad(
            model.up.kstencil, model.up.overlaps, Xup, Yup, model.up.frozen_bands
        )
        GXdn, GYdn = omega_grad(
            model.dn.kstencil, model.dn.overlaps, Xdn, Ydn, model.dn.frozen_bands
        )

        if λ != 0
            GOXup, GOYup, GOXdn, GOYdn = omega_updn_grad(model, Xup, Yup, Xdn, Ydn)
            GXup += λ * GOXup
            GYup += λ * GOYup
            GXdn += λ * GOXdn
            GYdn += λ * GOYdn
        end

        n = nw^2
        for ik in 1:nk
            for (i, v) in enumerate(view(GXup, :, :, ik))
                G[i, ik] = v
            end
            for (i, v) in enumerate(view(GYup, :, :, ik))
                G[n + i, ik] = v
            end
            for (i, v) in enumerate(view(GXdn, :, :, ik))
                G[n_inner + i, ik] = v
            end
            for (i, v) in enumerate(view(GYdn, :, :, ik))
                G[n_inner + n + i, ik] = v
            end
        end

        return nothing
    end

    return f, g!
end

function _build_fg_mag_disentangle_center(problem::LocalizationProblem)
    center_term = _find_center_term(problem.terms)
    isnothing(center_term) && error("CenterConstraintTerm is required for :mag_disentangle_center")
    model = problem.model
    λs = get(problem.solver_options, :lambda, 1.0)

    nb = n_bands(model.up)
    nw = n_wannier(model.up)
    nk = n_kpoints(model.up)
    n_inner = nb * nw + nw^2

    function f(XY)
        XY = reshape(XY, (2 * n_inner, nk))
        XYup = @view XY[1:n_inner, :]
        XYdn = @view XY[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, nb, nw)
        Xdn, Ydn = XY_to_X_Y(XYdn, nb, nw)
        Ωup = omega_center(
            omega(model.up.kstencil, model.up.overlaps, Xup, Yup);
            r₀ = center_term.r0,
            λ = center_term.λ,
        ).Ωt
        Ωdn = omega_center(
            omega(model.dn.kstencil, model.dn.overlaps, Xdn, Ydn);
            r₀ = center_term.r0,
            λ = center_term.λ,
        ).Ωt
        if λs == 0
            Ωupdn = 0
        else
            Ωupdn = omega_updn(model, X_Y_to_U(Xup, Yup), X_Y_to_U(Xdn, Ydn))
        end
        return Ωup + Ωdn + λs * Ωupdn
    end

    function g!(G, XY)
        XY = reshape(XY, (2 * n_inner, nk))
        XYup = @view XY[1:n_inner, :]
        XYdn = @view XY[(n_inner + 1):end, :]
        Xup, Yup = XY_to_X_Y(XYup, nb, nw)
        Xdn, Ydn = XY_to_X_Y(XYdn, nb, nw)
        GXup, GYup = omega_grad(
            center_penalty(center_term.r0, center_term.λ),
            model.up.kstencil,
            model.up.overlaps,
            Xup,
            Yup,
            model.up.frozen_bands,
        )
        GXdn, GYdn = omega_grad(
            center_penalty(center_term.r0, center_term.λ),
            model.dn.kstencil,
            model.dn.overlaps,
            Xdn,
            Ydn,
            model.dn.frozen_bands,
        )

        if λs != 0
            GOXup, GOYup, GOXdn, GOYdn = omega_updn_grad(model, Xup, Yup, Xdn, Ydn)
            GXup += λs * GOXup
            GYup += λs * GOYup
            GXdn += λs * GOXdn
            GYdn += λs * GOYdn
        end

        n = nw^2
        for ik in 1:nk
            G[1:n, ik] = vec(view(GXup, :, :, ik))
            G[(n + 1):n_inner, ik] = vec(view(GYup, :, :, ik))
            G[(n_inner + 1):(n_inner + n), ik] = vec(view(GXdn, :, :, ik))
            G[(n_inner + n + 1):end, ik] = vec(view(GYdn, :, :, ik))
        end

        return nothing
    end

    return f, g!
end

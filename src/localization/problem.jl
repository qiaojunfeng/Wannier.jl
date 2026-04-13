struct LocalizationProblem{P <: AbstractPenalty, M, O}
    penalty::P
    model::M
    mode::Symbol
    solver_options::O
end

function LocalizationProblem(p::AbstractPenalty, model, mode::Symbol; kwargs...)
    return LocalizationProblem(p, model, mode, (; kwargs...))
end

function build_fg!(problem::LocalizationProblem)
    if problem.mode == :maxloc
        return _build_fg_maxloc(problem)
    elseif problem.mode == :disentangle
        return _build_fg_disentangle(problem)
    elseif problem.mode == :mag_disentangle
        return _build_fg_mag_disentangle(problem)
    elseif problem.mode == :mag_disentangle_center
        return _build_fg_mag_disentangle_center(problem)
    else
        error("Unsupported localization mode: $(problem.mode)")
    end
end

function _build_fg_maxloc(problem::LocalizationProblem)
    p = problem.penalty
    model = problem.model
    cache = Cache(model)

    function fg!(F, G, U)
        compute_MU_UtMU!(cache, model.kstencil, model.overlaps, U)

        if G !== nothing
            cache.G = G
            omega_grad!(p, cache, model.kstencil, model.overlaps)
        end
        if F !== nothing
            return omega!(p, cache, model.kstencil, model.overlaps).Ω
        end
    end

    return fg!
end

function _build_fg_disentangle(problem::LocalizationProblem)
    p = problem.penalty
    model = problem.model
    cache = Cache(model)

    function fg!(Ω, G, XY)
        X, Y = XY_to_X_Y!(cache.X, cache.Y, XY)
        U = X_Y_to_U!(cache.U, X, Y)
        compute_MU_UtMU!(cache, model.kstencil, model.overlaps, U)

        if G !== nothing
            G_ = omega_grad!(p, cache, model.kstencil, model.overlaps)
            GX, GY = GU_to_GX_GY(G_, X, Y, model.frozen_bands)

            n = n_wannier(model)^2

            @inbounds for ik in 1:n_kpoints(model)
                for i in eachindex(GX[ik])
                    G[i, ik] = GX[ik][i]
                end
                for i in eachindex(GY[ik])
                    G[n + i, ik] = GY[ik][i]
                end
            end
        end
        if Ω !== nothing
            return omega!(p, cache, model.kstencil, model.overlaps).Ω
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
            for (i, v) in enumerate(GXup[ik])
                G[i, ik] = v
            end
            for (i, v) in enumerate(GYup[ik])
                G[n + i, ik] = v
            end
            for (i, v) in enumerate(GXdn[ik])
                G[n_inner + i, ik] = v
            end
            for (i, v) in enumerate(GYdn[ik])
                G[n_inner + n + i, ik] = v
            end
        end

        return nothing
    end

    return f, g!
end

function _build_fg_mag_disentangle_center(problem::LocalizationProblem)
    p = problem.penalty
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
        Ωup = omega(p, model.up.kstencil, model.up.overlaps, Xup, Yup).Ωt
        Ωdn = omega(p, model.dn.kstencil, model.dn.overlaps, Xdn, Ydn).Ωt
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
            p, model.up.kstencil, model.up.overlaps, Xup, Yup, model.up.frozen_bands
        )
        GXdn, GYdn = omega_grad(
            p, model.dn.kstencil, model.dn.overlaps, Xdn, Ydn, model.dn.frozen_bands
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
            G[1:n, ik] = vec(GXup[ik])
            G[(n + 1):n_inner, ik] = vec(GYup[ik])
            G[(n_inner + 1):(n_inner + n), ik] = vec(GXdn[ik])
            G[(n_inner + n + 1):end, ik] = vec(GYdn[ik])
        end

        return nothing
    end

    return f, g!
end

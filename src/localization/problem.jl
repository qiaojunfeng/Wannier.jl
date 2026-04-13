struct LocalizationProblem{P <: AbstractPenalty, M <: Model, O}
    penalty::P
    model::M
    mode::Symbol
    solver_options::O
end

function LocalizationProblem(p::AbstractPenalty, model::Model, mode::Symbol; kwargs...)
    return LocalizationProblem(p, model, mode, (; kwargs...))
end

function build_fg!(problem::LocalizationProblem)
    if problem.mode == :maxloc
        return _build_fg_maxloc(problem)
    elseif problem.mode == :disentangle
        return _build_fg_disentangle(problem)
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

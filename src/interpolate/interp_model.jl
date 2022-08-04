struct InterpolationModel{T<:Real}
    # model storing results of Wannierization
    model::Model{T}

    # R vectors for Fourier transform
    Rvectors::Union{RVectors{T},RVectorsMDRS{T}}

    # kpoint path for band structure
    kpath::KPath
end

function Base.getproperty(x::InterpolationModel, sym::Symbol)
    if sym ∈ fieldnames(Model)
        return getfield(x.model, sym)
    else
        # fallback to getfield
        getfield(x, sym)
    end
end

function Base.show(io::IO, model::InterpolationModel)
    show(io, model.model)
    println(io, "\n")

    show(io, model.Rvectors)
    println(io, "\n")

    show(io, "text/plain", model.kpath)
    return nothing
end

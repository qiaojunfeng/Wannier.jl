
struct Cache{M<:AbstractMatrix}
    isolated::Bool
    M::Vector{Vector{M}}
    X::Vector{M}
    Y::Vector{M}
    U::Vector{M}
    # n_bands x n_wann x n_kpts
    G::Array{Complex{T},3}
    r::Vector{Vec3}
    UtMU::Vector{Vector{M}}
    MU::Vector{Vector{M}}
end

function Cache(model::Model) end

# TODO check this after refactor of interpolation code
function getfirst(f::Function, itr)
    id = findfirst(f, itr)
    return id === nothing ? id : itr[id]
end

include("const.jl")
include("type.jl")
include("rgrid.jl")

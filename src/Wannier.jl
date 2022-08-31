module Wannier

using Requires

include("common/include.jl")
include("util/include.jl")
include("bvector.jl")
include("model.jl")
include("spread.jl")
include("io/include.jl")
include("wannierize/include.jl")
include("interpolate/include.jl")
include("realspace/include.jl")
include("cli/main.jl")

function __init__()
    @require PlotlyJS = "f0f68f2c-4968-5e81-91da-67840de0976a" begin
        include("plot/plotlyjs.jl")
    end

    @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
        include("plot/makie.jl")
    end
end

end

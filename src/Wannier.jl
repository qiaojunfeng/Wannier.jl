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

        # 3D plot of realspace WFs
        # @require StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91" begin
        #     @require GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326" begin
        #         @require Meshing = "e6723b4c-ebff-59f1-b4b7-d97aa5274f73" begin
        #             include("plot/realspace/wf_plotlyjs.jl")
        #         end
        #     end
        # end
    end

    # @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
    #     include("plot/makie.jl")
    # end
end

end

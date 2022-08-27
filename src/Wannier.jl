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

# function __init__()
#     @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
#         include("plot/include.jl")
#     end
# end

end

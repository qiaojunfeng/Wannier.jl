using Reexport
using WannierIO

#=
Note for some functions, e.g., `read_nnkp`, they are already defined in WannierIO.jl,
however, in Wannier.jl, we have some wrappers for these functions, which return a
`BVectors` instead of `NamedTuple`, so that they are better integrated with the rest of
the Wannier.jl.
Luckily, it seems if I export the new definitions, the old ones from WannierIO.jl are
hidden although I reexport them here.
=#
@reexport using WannierIO

include("w90/nnkp.jl")
include("w90/amn.jl")
include("w90/chk.jl")
include("w90/model.jl")
include("w90/band.jl")
include("w90/tb.jl")

include("truncate.jl")

include("volume/xsf.jl")
include("volume/cube.jl")

module Wannier

using Requires

include("common/const.jl")
include("common/type.jl")
include("util/io.jl")

include("util/linalg.jl")
export get_recip_lattice, get_lattice
export orthonorm_lowdin, eyes_A, rotate_M
export isunitary

include("util/kpoint.jl")
export get_kpoint_mappings

include("util/element.jl")
include("util/center.jl")

include("bvector.jl")
export get_bvectors

include("model.jl")

include("spread.jl")
export Spread, omega, omega_grad, center

include("io/w90.jl")

include("io/model.jl")
export write_w90

include("io/truncate.jl")
export truncate_w90

include("io/cube.jl")
export read_cube
include("io/wavefunction.jl")

include("wannierize/max_localize.jl")
export max_localize

include("wannierize/disentangle.jl")
export disentangle

include("wannierize/opt_rotate.jl")
export opt_rotate, rotate_A

include("wannierize/parallel_transport/parallel_transport.jl")
export parallel_transport

include("wannierize/split.jl")
export split_wannierize, split_unk, split_model

include("wannierize/constraint_center/max_localize.jl")
export max_localize_center

include("wannierize/constraint_center/disentangle.jl")
export disentangle_center

include("cli/main.jl")

include("interpolate/rvectors.jl")
include("interpolate/kpath.jl")
include("interpolate/band.jl")

# function __init__()
#     @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
include("plot.jl")
export plot_band, plot_band!

# include("plot/parallel_transport.jl")
#     end
# end

end

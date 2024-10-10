using Comonicon

# I have to use a different name for the CLI subcommands,
# it seems Comonicon returns wrong doc when there are multiple
# functions having the same name.

include("max_localize.jl")
include("disentangle.jl")
include("parallel_transport.jl")
include("opt_rotate.jl")
include("split_wannierize.jl")
include("truncate.jl")
include("band.jl")
include("fermisurf.jl")

"""
Julia package for Wannier functions.
"""
@Comonicon.main

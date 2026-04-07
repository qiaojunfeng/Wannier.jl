"""
    bandplot(kpath, eigenvals; kwargs...)

Plot the band structure given a path of k-points and the corresponding eigenvalues.
The kpath can be specified in three ways:
1. As a `KPath` object
2. As a vector of cumulative distances along the kpath (in units of 1/L, where L is unit of length),
    together with the indices and labels of high-symmetry kpoints
"""
function bandplot end

function bandplot! end

"""
Current Makie recipe does not allow setting up Axis labels and ticks.
This is a custom wrapper function to create the FigAxisPlot object and set the
labels and ticks, then call the bandplot! recipe to plot the bands.
"""
function get_bandplot end


"""
    projbandplot(kpath, eigenvals, projs, labels; kwargs...)

The `projs` is the projectability onto each orbital, i.e., a real-valued vector
of matrices, can be obtained from the gauge matrices by
```julia
nprojs = size(gauges[1], 2)
projs = Wannier.projectability(gauges)
```
"""
function projbandplot end

function projbandplot! end

function get_projbandplot end

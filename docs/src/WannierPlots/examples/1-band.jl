# # 1. Band structure of silicon

#=
```@meta
CurrentModule = WannierPlots
```
=#

#=
In this tutorial, we will use plot the band structure.

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`1-band.ipynb`](./1-band.ipynb)
    - Julia script: [`1-band.jl`](./1-band.jl)
=#

# ## Preparation
# Load the package
using Wannier: read_w90_band
using WannierPlots

# Path of current tutorial
PWD = "3-band"

# Read bands data files
recip_lattice = [
    -1.15701 1.15701 1.15701
    1.15701 -1.15701 1.15701
    1.15701 1.15701 -1.15701
]
kpi, E = read_w90_band("$PWD/si2", recip_lattice)
#=
!!! tip

    Here I pass a `recip_lattice` to the `read_w90_band` function,
    so that it will return a tuple of `(KPathInterpolant, Matrix)`.
    You can also call the `read_w90_band` function without `recip_lattice`,
    however, this "raw" version will return rather verbose outputs,
    not very handy for usage. See the API [`read_w90_band`](@ref) for details.
=#

# Now let's have a look at the interpolated band
P = plot_band(kpi, E)
Main.HTMLPlot(P, 500) # hide

# Finally, do the comparison
# here I generate a fake band structure by shifting it upwards by 0.1 eV
E_shift = E .+ 0.1
P = plot_band_diff(kpi, E, E_shift)
Main.HTMLPlot(P, 500) # hide

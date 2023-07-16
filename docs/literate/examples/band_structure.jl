# # 3. Band structure

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In this tutorial, we will use Wananier interpolation to calculate the band structure of silicon valence and conduction manifold. A bit different from previous tutorials, we will

1. construct a [`InterpModel`](@ref) by reading the `win`, `mmn`, `eig`, and `chk.fmt` files
2. run [`interpolate`](@ref) on the `InterpModel` to compute band structure
3. read `Wannier90` interpolated `band.dat` and compare with `Wannier.jl` interpolated bands

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`tutorial.ipynb`](./tutorial.ipynb)
    - Julia script: [`tutorial.jl`](./tutorial.jl)
=#

# ## Preparation
# Load the package
using Wannier
using WannierPlots
#=
The `WannierPlots` is a companion package of `Wannier.jl`, for plotting band structures,
real space WFs, etc.
This separation reduces the compiling latency of `Wannier.jl`, and avoids
introduing too many package dependencies that are irrelevant to Wannierization.
=#

#=
!!! tip

    Use the `run.sh` script which automate the scf, nscf, pw2wannier90 steps.
=#

#=
## Model generation

### from `chk.fmt` file

We will use the [`read_w90_interp`](@ref) function to read the
`win`, `mmn`, `eig`, and `chk.fmt` files, and construct an
[`InterpModel`](@ref) that is used for interpolation purpose.
=#
model = read_w90_interp("si2")

#=
!!! tip

    To avoid a bloated `Model` that contains everything, and to
    ["Do One Thing And Do It Well"](https://en.wikipedia.org/wiki/Unix_philosophy#Do_One_Thing_and_Do_It_Well),
    we separate on purpose the `Model` that is solely for
    Wannierization, and the `InterpModel`, that is only used
    for Wannier interpolation of operators.
    This is convenient for developers to focus on the
    the Wannierization or interpolation algorithm without
    being distracted by the other part.
=#
#=
### from `amn` file

You can also use the [`read_w90`](@ref) function to read the
`amn` file and generate a [`Model`](@ref),
=#
m = read_w90("si2")
#=
Then use the [`InterpModel`](@ref) constructor to
construct an `InterpModel` from an existing `Model`,
=#
m = Wannier.InterpModel(m)
#=
However, there are some differences:
1. the `amn` gauge is used, instead of that from `chk`
2. the kpoint path for band structure is auto generated from
    the lattice, instead of using that in `win` file (if found)

So, it is recommended to use the `read_w90_interp` function,
or you run a `max_localize` or `disentangle` on the `Model`
to smooth the gauge, then construct an `InterpModel`.
=#

#=
### Tips on kpath

The `read_w90_interp` function will parse the `kpoint_path` block in the `win` file:
1. if the `kpoint_path` block is found, it will use that path
2. if the `kpoint_path` block is not found, the `InterpModel`
    constructor will auto generate a kpath from the lattice,
    by using the [`get_kpath`](@ref) function

!!! tip

    The [`get_kpath`](@ref) works for arbitrary lattice, e.g.,
    either conventional or primitive, it will generate the correct kpath.

During the band interpolation, the [`interpolate`](@ref)
function will call the [`interpolate_w90`](@ref) function to
generate an exactly identical kpath to that of `Wannier90`.
=#

#=
### Tips on interpolation algorithm

There are two interpolation algorithms
1. Wigner-Seitz (WS) interpolation
2. Minimal-distance replica selection (MDRS) method

The `read_w90_interp` function will parse the `use_ws_distance` parameter in
the `win` file:
1. if `use_ws_distance` is `true`, it will use the MDRS
2. if `use_ws_distance` is `false`, it will use the WS
3. if not found, it will use the MDRS

You can also control the interpolation algorithm by the constructor
[`InterpModel(model::Model; mdrs::Bool=true)`](@ref InterpModel(model::Model; mdrs::Bool=true)).

Moreover, there are two versions of MDRS, i.e., `mdrsv1` and `mdrsv2` in the
jargon of `Wannier.jl`. The `v2` is faster than `v1` so it is used by default.
For details, see the API for interpolation,
e.g., [`fourier`](@ref), [`invfourier`](@ref),
[`_fourier_mdrs_v1`](@ref), [`_fourier_mdrs_v2`](@ref), [`_invfourier_mdrs_v1`](@ref), and
[`_invfourier_mdrs_v2`](@ref), etc.

If you are tech-savvy, you can even control the generation of `R` vectors for interpolation,
this is useful for developing new
interpolation algorithms, see
[`get_Rvectors_mdrs`](@ref), [`get_Rvectors_ws`](@ref), etc.
=#

#=
## Band interpolation

Computing band structure is very easy, by calling
the [`interpolate`](@ref) function,
it will return a `KPathInterpolant`, and the
band eigen energies,
=#
kpi, E = interpolate(model)

#=
The `kpi` stores the specific kpoint coordinates along the kpath,
while the `InterpModel.kpath` only stores the high-symmetry kpoints
their labels, that's why we need to return the `kpi` object
=#
kpi

#=
## Plotting band structure

### Saving interpolated band

You can save the result to the same format
as `Wannier90` `band.dat`, by
=#
write_w90_band("wjl", kpi, E)
#=
where `wjl` is the seedname of the output,
i.e., written files are
- `wjl_band.kpt`
- `wjl_band.dat`
- `wjl_band.labelinfo.dat`
=#

#=
### Visualization in the Julia world

Instead of saving, you can also plot the band structure by
calling the `WannierPlots.plot_band` function,

!!! note

    This requires you having executed `using PlotlyJS`
=#
P = plot_band(kpi, E)
Main.HTMLPlot(P, 500) # hide
#=
Or, you can use the plotting functions provided by
[`Brillouin.jl`](https://thchr.github.io/Brillouin.jl/stable/kpaths/#Band-structure),
but requires a little bit transposation to size `(n_kpts, n_bands)`
=#
P = plot(kpi, E')
Main.HTMLPlot(P, 500) # hide

#=
## Comparing band structures

Now we load the `Wannier90` interpolated band,
to compare between the two codes,
=#
kpi_w90, E_w90 = read_w90_band("si2", model.recip_lattice)
#=
!!! tip

    Here I pass a `recip_lattice` to the `read_w90_band` function,
    so that it will return a tuple of `(KPathInterpolant, Matrix)`.
    You can also call the `read_w90_band` function without `recip_lattice`,
    however, this "raw" version will return rather verbose outputs,
    not very handy for usage. See the API [`read_w90_band`](@ref) for details.
=#

# Now let's have a look at `Wannier90` interpolated band
P = plot_band(kpi_w90, E_w90)
Main.HTMLPlot(P, 500) # hide

# Finally, do the comparison
P = plot_band_diff(kpi, E_w90, E)
Main.HTMLPlot(P, 500) # hide

#=
As expected, the two band structures exactly overlaps. 🥳
=#

#=
Well done! You have finished the first Wannier interpolation.
=#

# # Wannier interpolation of band structure

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In this tutorial, we will use Wananier interpolation to compute the band structure
of silicon valence + conduction bands. A bit different from previous tutorials, we will

1. Construct a [`TBHamiltonian`](@ref) by reading the `win`, `mmn`, `eig`, and `chk` files
2. Run [`interpolate`](@ref) on the `TBHamiltonian` to compute band structure
3. Read wannier90 interpolated `band.dat` and compare with Wannier.jl interpolated bands
=#

# ## Preparation
# Load the package
using Wannier
using Wannier.Datasets
using WannierPlots
#=
The `WannierPlots` is a companion package of `Wannier.jl`, for plotting band structures,
real space WFs, etc.
=#

#=
## Model construction

We will use the [`read_w90_with_chk`](@ref) function to read the `win`, `mmn`,
`eig`, and `chk` files, so that the U matrix is already the optimized one.
=#
model = read_w90_with_chk(dataset"Si2/Si2", dataset"Si2/reference/Si2.chk")

# and check the spread to make sure our `Model` is sensible
omega(model)

# Now construct a [`TBHamiltonian`](@ref) which contains the ``\mathbf{R}``-space
# Wannier Hamiltonian, ``H(\mathbf{R})``.
H = Wannier.TBHamiltonian(model)

#=
!!! tip

    To avoid a bloated `Model` that contains everything, and to
    ["Do One Thing And Do It Well"](https://en.wikipedia.org/wiki/Unix_philosophy#Do_One_Thing_and_Do_It_Well),
    we separate on purpose the `Model` that is solely for Wannierization, and
    the `TBHamiltonian`, that is only used for Wannier interpolation of
    tight-binding Hamiltonian.
    This is also convenient for developers to focus on the the Wannierization or
    interpolation algorithm without being distracted by one or the other.

## Band-structure kpoint path

There are two possible ways to generate a kpath for band-structure interpolation.

### From `win` file `kpoint_path` block

First read the `win` file,
=#
win = read_win(dataset"Si2/Si2.win")
# the returned `win` is a `NamedTuple` that contains all the input tags in the `win` file.
# Then generate a `KPath` based on crystal structure and `kpoint_path` block,
kpath = get_kpath(win.unit_cell_cart, win.kpoint_path)

#=
### Auto generate kpath from lattice

Another approach is to auto generate a kpath from the lattice, which can be
either conventional or primitive, the function [`get_kpath`](@ref) will generate
a correct kpath.
=#
kpath_auto = get_kpath(model)

#=
### Set kpoint spacing along the kpath

To help comparing with the `band.dat` file, we provide a function
[`interpolate_w90`](@ref) which will return the exact same kpoints as that
in wannier90 `prefix_band.kpt` file.
The function returns a `KPathInterpolant` object, containing a list of kpoint
coordinates to be interpolated on,
=#
kpi = Wannier.interpolate_w90(kpath, 100)

#=
!!! tip

    In comparison, the `KPath` only stores the high-symmetry kpoints and their labels.

## Band interpolation

### Tips on interpolation algorithm

There are two interpolation algorithms
1. Wigner-Seitz (WS) interpolation
2. Minimal-distance replica selection (MDRS) method

In wannier90 input file, the `use_ws_distance` parameter determines which algorithm
to use, and the default is `true`:
1. if `use_ws_distance` is `true`, use MDRS
2. if `use_ws_distance` is `false`, use WS

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
Computing band structure is very easy, by calling
the [`interpolate`](@ref) function,
it will return a `KPathInterpolant`, and the
band eigen energies,
=#
eigenvalues = interpolate(H, kpi)

#=
## Plotting band structure

### Saving interpolated band

You can save the result to the same format
as `Wannier90` `band.dat`, by
=#
write_w90_band("wjl", kpi, eigenvalues)
#=
where `wjl` is the prefix of the output,
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
P = plot_band(kpi, eigenvalues)
Main.HTMLPlot(P, 500) # hide
#=
Or, you can use the plotting functions provided by
[`Brillouin.jl`](https://thchr.github.io/Brillouin.jl/stable/kpaths/#Band-structure),
but requires a little bit transposation to size `(n_kpts, n_bands)`
=#
P = plot(kpi, eigenvalues)
Main.HTMLPlot(P, 500) # hide

#=
## Comparing band structures

Now we load the `Wannier90` interpolated band,
to compare between the two codes,
=#
kpi_w90, eigenvalues_w90 = read_w90_band(dataset"Si2/reference/Si2", model.recip_lattice)
#=
!!! tip

    Here I pass a `recip_lattice` to the `read_w90_band` function,
    so that it will return a tuple of `(KPathInterpolant, Matrix)`.
    You can also call the `read_w90_band` function without `recip_lattice`,
    however, this "raw" version will return rather verbose outputs,
    not very handy for usage. See the API [`read_w90_band`](@ref) for details.
=#

# and compare the two band structures,
P = plot_band_diff(kpi, eigenvalues_w90, eigenvalues)
Main.HTMLPlot(P, 500) # hide

# Finally, we can also compare with DFT bands
using WannierIO
qe = WannierIO.read_qe_xml(dataset"Si2/reference/qe_bands.xml")
P = plot_band_diff(kpi, qe.eigenvalues, eigenvalues)
Main.HTMLPlot(P, 500) # hide
#=
As expected, the two band structures exactly overlaps. 🥳
=#

#=
Well done! You have finished the first Wannier interpolation.
=#

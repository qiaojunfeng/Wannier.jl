# # Wannier interpolation of band structure

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In this tutorial, we will use Wananier interpolation to compute the band structure
of silicon valence + conduction bands. A bit different from previous tutorials, we will

1. Construct a [`Model`](@ref) by reading the `win`, `mmn`, `eig`, and `chk` files
2. Construct a tight-binding Hamiltonian [`HamiltonianRspace`](@ref) from the `Model`
2. Run [`interpolate`](@ref) on the `TBHamiltonian` to compute band structure
3. Read wannier90 interpolated `band.dat` and compare with our interpolated bands
=#

# ## Preparation
# Load the package
using Wannier
using Wannier.Datasets

#=
## Model construction

We will use the [`read_w90_with_chk`](@ref) function to read the `win`, `mmn`,
`eig`, and `chk` files, so that the U matrix is already the optimized one.
=#
model = read_w90_with_chk(dataset"Si2/Si2", dataset"Si2/reference/Si2.chk")

# and check the spread to make sure our `Model` is sensible
omega(model)

# Now construct a tight-binding Hamiltonian, the [`TBHamiltonian`](@ref) function
# returns a [`TBOperator`](@ref) struct, which contains the ``\mathbf{R}``-space
# Wannier Hamiltonian, ``H(\mathbf{R})``.
H = TBHamiltonian(model)

#=
!!! tip

    To avoid a bloated `Model` that contains everything, and to
    ["Do One Thing And Do It Well"](https://en.wikipedia.org/wiki/Unix_philosophy#Do_One_Thing_and_Do_It_Well),
    we separate on purpose the `Model` that is solely for Wannierization, and
    the `TBOperator`, that is only used for Wannier interpolation of
    tight-binding operators.
    This is also convenient for developers to focus on the the Wannierization or
    interpolation algorithm without being distracted by one or the other.

## Band-structure kpoint path

There are two possible ways to generate a kpath for band-structure interpolation.

### From `win` file `kpoint_path` block

First read the `win` file,
=#
win = read_win(dataset"Si2/Si2.win")
# the returned `win` is a `NamedTuple` that contains all the input tags in the `win` file.
#
# Then generate a `KPath` based on crystal structure and `kpoint_path` block,
kpath = generate_kpath(win.unit_cell_cart, win.kpoint_path)

#=
### Auto generate kpath from lattice

Another approach is to auto generate a kpath from the lattice, which can be
either conventional or primitive, the function [`generate_kpath`](@ref) will generate
a correct kpath.
=#
kpath_auto = generate_kpath(model)

#=
### Set kpoint spacing along the kpath

To help comparing with the `band.dat` file, we provide a function
[`generate_w90_kpoint_path`](@ref) which will return the exact same kpoints as that
in wannier90 `prefix_band.kpt` file.
The function returns a `KPathInterpolant` object, containing a list of kpoint
coordinates to be interpolated on,
=#
kpi = Wannier.generate_w90_kpoint_path(kpath)

# you can also directly pass the inputs in `win` file to directly genereate
# the kpoints,
kpi = generate_w90_kpoint_path(win.unit_cell_cart, win.kpoint_path)

#=
!!! tip

    In comparison, the `KPath` only stores the high-symmetry kpoints and their
    labels, while the `KPathInterpolant` stores the kpoint coordinates.

## Band interpolation

Computing band structure is very easy, we first construct a
[`HamiltonianInterpolator`] from the `hamiltonian`,
=#
interp = HamiltonianInterpolator(hamiltonian)

#=
the returned `interp` is a functor, i.e., under the hood is a Julia `struct`
but can be called as a function: in our case, we can pass kpoint coordinates
to the interpolator and it will return the interpolated eigenvalues and eigenvectors.
We can either pass a vector of 3-vectors for fractional coordinates, or directly
a `KPathInterpolant` object,
=#
E, V = interp(kpi)

#=
## Plotting band structure

### Saving interpolated band

You can save the result to the same format
as `Wannier90` `band.dat`, by
=#
write_w90_band("wjl", kpi, E)
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
calling the [`Wannier.plot_band`](@ref) function.

To activate the `plot_band` function, we need to first load `PlotlyJS` package,
=#
using PlotlyJS

# then we can plot the band structure by
P = plot_band(kpi, E; win.fermi_energy)
Main.HTMLPlot(P, 500) # hide
#=
Or, you can use the plotting functions provided by
[`Brillouin.jl`](https://thchr.github.io/Brillouin.jl/stable/kpaths/#Band-structure),
but requires a different memory layout: a vector of length-`n_kpoints`, each
elmenet is a length-`n_bands` vector
=#
E_t = eachrow(reduce(hcat, E))
P = plot(kpi, E_t)
Main.HTMLPlot(P, 500) # hide

#=
## Comparing band structures

Now we load the `Wannier90` interpolated band,
to compare between the two codes,
=#
kpi_w90, E_w90 = read_w90_band(dataset"Si2/reference/MDRS/Si2", reciprocal_lattice(model))
#=
!!! tip

    Here I pass a `recip_lattice` to the `read_w90_band` function,
    so that it will return a tuple of `(KPathInterpolant, Matrix)`.
    You can also call the `read_w90_band` function without `recip_lattice`,
    however, this "raw" version will return rather verbose outputs,
    not very handy for usage. See the API [`read_w90_band`](@ref) for details.
=#

# and compare the two band structures,
P = plot_band_diff(kpi, E_w90, E)
Main.HTMLPlot(P, 500) # hide

# Finally, we can also compare with DFT bands
using WannierIO
qe = WannierIO.read_qe_xml(dataset"Si2/reference/qe_bands.xml")
P = plot_band_diff(kpi, qe.eigenvalues, E)
Main.HTMLPlot(P, 500) # hide
#=
As expected, the Wannier-interpolated band structures nicely reproduce
DFT bands.
=#

#=
Well done! You have finished the first Wannier interpolation 🥳.
=#

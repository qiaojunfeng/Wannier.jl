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
2. Construct a [`TBHamiltonian`](@ref) from the `Model`
3. Run [`HamiltonianInterpolator`](@ref) to compute band structure
4. Read wannier90 interpolated `band.dat` and compare with our interpolated bands
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
model = read_w90_with_chk(dataset"Si2/Si2", dataset"Si2/outputs/Si2.chk")

# and check the spread to make sure our `Model` is sensible
spread(model)

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
# the returned `win` is an `OrderedDict` that contains all the input tags in the `win` file.
#
# Then generate a `KPath` based on crystal structure and `kpoint_path` block:
# `KSegment` holds the segment definitions, `KPath` samples the kpoints on them.
# Passing `default_w90_kpath_num_points()` (=100) matches Wannier90's spacing.
kseg = KSegment(reciprocal_lattice(win["unit_cell_cart"]), win["kpoint_path"])
kpath = KPath(kseg, default_w90_kpath_num_points())

#=
### Auto generate kpath from lattice

Another approach is to auto-generate a kpath from the crystal structure
stored in the `Model` (uses spacegroup symmetry via `Brillouin.jl` and `Spglib.jl`):
=#
using Spglib
kpath_auto = KPath(KSegment(model), default_w90_kpath_num_points())

#=
## Band interpolation

Computing band structure is very easy, we first construct a
[`HamiltonianInterpolator`](@ref) from the Hamiltonian,
=#
interp = HamiltonianInterpolator(H)

#=
the returned `interp` is a functor that accepts a vector of fractional kpoint
coordinates and returns `(eigenvalues, eigenvectors)`.
`collect(kpath)` materialises the `KPath` into that vector.
The eigenvalues matrix has shape `n_bands × n_kpoints`; `eachcol` converts it
to the `Vector{Vector{Float64}}` format expected by plotting and IO functions.
=#
E_mat, V = interp(collect(kpath))
E = collect(eachcol(E_mat))

#=
## Plotting band structure

### Saving interpolated band

You can save the result to the same format
as `Wannier90` `band.dat`, by
=#
write_w90_band("wjl", kpath, E)
#=
where `wjl` is the prefix of the output,
i.e., written files are
- `wjl_band.kpt`
- `wjl_band.dat`
- `wjl_band.labelinfo.dat`
=#

#=
### Visualization in the Julia world

Instead of saving, you can also plot the band structure using the
[`get_bandplot`](@ref) function, which is activated by loading any
Makie backend:
=#
using WGLMakie

# then plot the band structure by
fig, ax, plt = get_bandplot(kpath, E; fermi_energy = win["fermi_energy"])
fig

#=
## Comparing band structures

Now we load the `Wannier90` interpolated band,
to compare between the two codes,
=#
kpath_w90, E_w90 = read_w90_band(dataset"Si2/outputs/MDRS/Si2", reciprocal_lattice(model))
#=
!!! tip

    Passing `recip_lattice` to [`read_w90_band`](@ref) returns a
    `(KPath, Vector{Vector{Float64}})` pair ready for plotting.
    Without it, the function returns verbose raw data.

The two-argument form of [`get_bandplot`](@ref) overlays both band structures:
=#
fig, ax, plt = get_bandplot(kpath_w90, E_w90, E;
    kwargs1 = (label = "Wannier90",),
    kwargs2 = (label = "Wannier.jl", linestyle = :dash),
)
fig

# Finally, compare with DFT bands from the QE XML output
using QuantumEspressoIO
qe = QuantumEspressoIO.read_pw_xml(dataset"Si2/outputs/qe_bands.xml")
fig, ax, plt = get_bandplot(kpath_w90, qe.eigenvalues, E;
    kwargs1 = (label = "QE DFT",),
    kwargs2 = (label = "Wannier.jl", linestyle = :dash),
)
fig
#=
As expected, the Wannier-interpolated band structures nicely reproduce
DFT bands.
=#

#=
Well done! You have finished the first Wannier interpolation 🥳.
=#

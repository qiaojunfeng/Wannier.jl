# # 9. Parallel transport using custom b-vectors

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In this tutorial, we will Wannierize a single isolated band of CuBr2,
by constructing the parallel transport gauge (PTG).

However, the CuBr2 system has a special set of b-vectors: some nearest neighbors
are not included. This breaks the requirement of the PTG, since it needs the
overlap matrices of nearest neighbors along 6 directions. Therefore, we need to
generate a custom set of b-vectors, write them to a `nnkp` file, rerun
`pw2wannier90.x`, and use the new `mmn` file to construct the PTG.

!!! warning

    The generated b-vectors do not satisfy the completeness condition, therefore
    it should not be used to compute the spread and other physical operators.

## Outline

1. construct a [`Model`](@ref), by reading the `win`, `amn`, `mmn`, and `eig` files
2. generate a new set of b-vectors containing only the 6 directions
3. write the new b-vectors to a `nnkp` file
4. launch `pw2wannier90.x` to compute a new `mmn` file
5. read the new `mmn` file, and construct the PTG
6. use the PTG gauge matrices in the `Model`
7. maxiaml localize to smoothen the gauge
8. interpolate band structure

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`9-parallel_transport_bvectors.ipynb`](./9-parallel_transport_bvectors.ipynb)
    - Julia script: [`9-parallel_transport_bvectors.jl`](./9-parallel_transport_bvectors.jl)
=#

# ## Preparation
# Load the package
using Wannier
using WannierPlots

# Path of current tutorial
CUR_DIR = "9-parallel_transport_bvectors"

#=
## Model generation

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct a [`Model`](@ref)

!!! tip

    Before you do the actual calculations, if you don't know the kpath of your
    structure, you can use the following functions to auto generate a kpath
    ```julia
    win = read_win("CuBr2.win")
    kp = get_kpath(win.unit_cell, win.atoms_frac, win.atom_labels)
    kpi = Wannier.interpolate_w90(kp, 100)
    Wannier.write_w90_kpt_label("CuBr2", kpi)
    ```
    The outputs are two files: `CuBr2_band.kpt` and `CuBr2_band.labelinfo.dat`,
    which are exactly in the same format as the `Wannier90` output bands, you
    can inspect these two files and use them to populate the `kpoint_path` block
    in the `win` file, or the `K_POINTS` block in the `pw.x` bands calculation.
=#
model = read_w90("$CUR_DIR/CuBr2")

#=
!!! note

    This is a simple system, you can even get good WFs with a simple
    `max_localize(model)`. However, we will use this tutorial to demonstrate
    the method to overcome the b-vector problem in PTG.
=#

# The initial spread is
omega(model)

#=
## Custom b-vectors

As shown here, the b-vectors are a bit strange that some nearest neighbors are not included,
=#
using PlotlyJS
WannierPlots.plot(model.bvectors)

# thus we need to use a different set of b-vectors, which simply includes all the
# 6 nearest neighbors.
bvectors_nn = Wannier.get_bvectors_nearest(model.kpoints, model.recip_lattice)
WannierPlots.plot(bvectors_nn)

# now write to a new `nnkp` file
n_wann = model.n_wann
exclude_bands = collect(1:16)
Wannier.write_nnkp("$CUR_DIR/CuBr2_nn.nnkp", bvectors_nn, n_wann, exclude_bands)

#=
now we can run `pw2wannier90.x` to generate a new `mmn` file.
A template input file for `pw2wannier90.x` is provided here

```fortran
&INPUTPP
    outdir = './out/'
    prefix = 'cubr2'
    seedname = 'CuBr2_nn'
    write_amn = .false.
    write_eig = .false.
    scdm_proj = .true.  ! this is useless, but we need it to avoid p2w error
/
```
=#

#=
## Parallel transport

Now, we load the new `mmn` file, and construct a new [`Model`](@ref)
=#
M_nn = WannierIO.read_mmn("$CUR_DIR/CuBr2_nn.mmn")[1];
model_nn = Wannier.Model(
    model.lattice,
    model.atom_positions,
    model.atom_labels,
    model.kgrid,
    model.kpoints,
    bvectors_nn,
    model.frozen_bands,
    M_nn,
    model.U,
    model.E,
)

# finally, we run parallel transport
U, _ = parallel_transport(model_nn)

#=
Since the b-vectors are not complete, computing spread on top of `model_nn`
is meaningless, we assign back the gauge matrices to `model`
=#
model.U .= U;

# and the new spread is
omega(model)

# then further maximal localize it
U = max_localize(model)
model.U .= U;
omega(model)

#=
## Band interpolation

Finally, let's have a look at the band interpolation.

First load the QE band structure
=#
kpoints_qe, E_qe = WannierIO.read_qe_band("$CUR_DIR/qe_bands.dat");
# the Fermi energy from scf calculation
ef = 4.6459

# Force using `kpoint_path` in `win` file
win = read_win("$CUR_DIR/CuBr2.win")
kpath = Wannier.KPath(win.unit_cell, win.kpoint_path)

interp_model = Wannier.InterpModel(model; kpath=kpath)

# interpolate band structure
# the QE bands use 50 points per segment, so we use 50 here as well
kpi = Wannier.interpolate_w90(kpath, 50)
E = Wannier.interpolate(interp_model, kpi)

# plot band difference
fig = plot_band_diff(kpi, E_qe, E; fermi_energy=ef)
fig.plot.layout.width = 1000
fig.plot.layout.height = 1000
fig.plot.layout.autosize = false
fig

#=
That's all!
=#

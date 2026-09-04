# # Manifold-remixed Wannier functions of TiO$_2$

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In this tutorial, we will use the [manifold-remixing method]()
to construct MLWFs for the valence and the conduction manifolds of TiO$_2$,
separately.

!!! warning "Tutorial pending migration"
    This tutorial still references several pre-rewrite helpers that no longer
    exist in the current API (`Wannier.get_U`, `Wannier.Tools.mrwf`,
    `Wannier.Tools.merge_gauge`, `has_cubic_neighbors`, `write_nnkp_cubic`, …).
    Read it as conceptual overview until the example is fully migrated.

## Outline

1. construct a [`Model`](@ref) for a Wannier90-Wannierized
    valence + conduction calculation of TiO$_2$
2. split the `Model` into two groups for valence and conduction, respectively,
    and get initial gauge matrices for the two groups through [`mrwf`](@ref)
3. maximally localize the two groups of WFs
=#

# ## Preparation
# Load the package
using Wannier
using Wannier.Datasets
using WannierPlots
# For plotting band structures
using PlotlyJS

#=
## Model construction

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct a [`Model`](@ref) that abstracts the calculation
=#
path = "/home/jqiao/git/WannierDatasets/datasets"
model = read_w90_with_chk("$path/TiO2/TiO2", "$path/TiO2/outputs/TiO2.chk")

# Always good to check the band structure
qe = WannierIO.read_qe_xml("$path/TiO2/outputs/qe_bands.xml");
kpi, w90 = read_w90_band("$path/TiO2/outputs/TiO2", model.recip_lattice);
P = plot_band_diff(kpi, qe.eigenvalues, w90; qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

#=
TiO$_2$ is a special case, where the default set of ``b``-vectors does not
contain 6 cubic nearest neighbors. This can be checked by
=#
Wannier.has_cubic_neighbors(model.kstencil)

# Alternatively, one can also pass a filename of a `nnkp` file to the function
Wannier.has_cubic_neighbors("$path/TiO2/outputs/TiO2.nnkp")

#=
Therefore, we need to generate a special `nnkp` file containing those 6 neighbors,
run `pw2wannier90.x` to generate another `mmn` file, and feed that to the
parallel transport algorithm.

The following function reads a `win` file and write a `nnkp` file containing 6
cubic neighbors.
=#
Wannier.write_nnkp_cubic("TiO2_cubic.nnkp", "$path/TiO2/TiO2.win")

#=
Since the 6 neighbors are not complete bvectors, we won't be able to compute a
proper weights for them. Therefore, we explicitly pass a `nothing` to the
2nd argument of `read_nnkp`, to construct a `KSpaceStencil` object with zero weights.
=#
kstencil_cubic = read_nnkp("$path/TiO2/inputs/TiO2_cubic.nnkp", nothing)
# And we read the overlap matrices for the cubic neighbors
overlaps_cubic = read_mmn("$path/TiO2/TiO2_cubic.mmn")[1];
#=
Finally, we can construct a special `Model` that can be used for parallel transport.
This model will reuse the lattice, atomic positions, and gauge matrices from the
original model, only the `KSpaceStencil` and the overlap matrices are changed.
=#
model_cubic = Wannier.Model(model, kstencil_cubic, overlaps_cubic)

#=
!!! warning

    The `model_cubic` is only used for parallel transport, one cannot use it
    to compute spreads since it does not contain complete set of `b`-vectors.

## Manifold-remixed Wannier functions

Then we split the model into two subgroups.
The first subgroup is the valence bands, with band indices starting from 1 to 16,
the second subgroup is the conduction bands, with band indices starting from 17
to 26 (which is the number of WFs).
=#
indices = [1:16, 17:26]

#=
Apart from calling the [`split_model`](@ref) function, we can also use a
convenience function [`mrwf`](@ref),
=#
(model_v, U_v), (model_c, U_c) = Wannier.Tools.mrwf(model, indices, model_cubic);
#=
Returned are two separated `Model`s, and the two corresponding gauge transformation
matrices from the total manifold to the two separated submanifolds.
This is useful if you want to further rotate the gauge of some operators,
or if you want to rotate the `UNK` files for plotting WFs,
in that case you need to pass these two matrices to the function
[`split_unk`](@ref).

The valence bands
=#
model_v
# and the gauge transformation matrix
size(U_v[1])
# The conduction bands
model_c
# and the gauge transformation matrix
size(U_c[1])

# and take a look at the spread
spread(model_v)

# and the conduction bands
spread(model_c)

#=
Since after parallel transport, the WFs are still not maximally localized WFs,
thus here their spreads are large. Later on we will run maximal localizations
to smoothen the gauge and get MLWFs.

The `U_v` or `U_c` matrices are the gauge transformation matrices from the
`mmn`/`eig` file of the total manifold to the `mmn`/`eig` files of the valence
or conduction manifolds in the respective folders, `val` or `cond`.
This can be verified by
=#
model_v_test = transform_gauge(model, U_v)
# and we need to copy the parallel transport gauge matrices
model_v_test.gauges .= model_v.gauges;
# the fields of the two models are the same
model_v_test ≈ model_v
# and compare the spreads with that of `model_v`
spread(model_v_test)

# and the conduction bands
model_c_test = transform_gauge(model, U_c)
# and we need to copy the parallel transport gauge matrices
model_c_test.gauges .= model_c.gauges;
# the fields of the two models are the same
model_c_test ≈ model_c
# and compare the spreads with that of `model_c`
spread(model_c_test)

#=
!!! tip

    Note that if you prefer working with files instead of `Model`s, you can
    also pass prefix for valence+conduction calculation, and the `mmn` file name
    of the cubic `mmn` file. The function will read the `chk` file of the
    valence+conduction calculation, and write the `amn`/`mmn`/`eig` files of
    valence or conduction in the `val` or `cond` folder, respectively.
=#
outdirs = ["val", "cond"]
Wannier.Tools.mrwf("$path/TiO2/TiO2", indices, outdirs, "$path/TiO2/TiO2_cubic.mmn")

#=
Specifically, the returned gauge matrices `U_v` and `U_c` are written to the
`val/TiO2_split.amn` and `cond/TiO2_split.amn` files, respectively.
=#
U_v = read_amn("$path/TiO2/val/TiO2_split.amn");
U_c = read_amn("$path/TiO2/cond/TiO2_split.amn");

#=
Now, we can call [`localize`](@ref) on the `model_v` and `model_c`, respectively.
But here as an demonstration, we run Wannier90 inside the `val` directory and
the `cond` directory. The `chk` file contains the final Wannierized
valence and conduction.
=#
U_v2 = Wannier.get_U(read_chk("$path/TiO2/outputs/val/TiO2.chk"));
U_c2 = Wannier.get_U(read_chk("$path/TiO2/outputs/cond/TiO2.chk"));

#=
Note that since each run of [`mrwf`](@ref) might give different gauge matrices,
we use `U_v` from previous result so that we can compare spreads wrt. that from Wannier90.
=#
model_v = transform_gauge(model, U_v);
spread(model_v, U_v2)

# here is the Wannier90 `wout`
open("$path/TiO2/outputs/val/TiO2.wout") do io
    print(join(readlines(io)[(end - 52):(end - 30)], "\n"))
end

# and the conduction bands
model_c = transform_gauge(model, U_c);
spread(model_c, U_c2)

# here is the Wannier90 `wout`
open("$path/TiO2/outputs/cond/TiO2.wout") do io
    print(join(readlines(io)[(end - 48):(end - 30)], "\n"))
end

# you can save the gauge matrices to an `amn` file
write_amn("val/TiO2.chk.amn", U_v2)
write_amn("cond/TiO2.chk.amn", U_c2)

#=
If you want to get the gauge matrices that convert the valence+conduction
`mmn`/`eig` to the final (maximally localized) valence MLWFs,
you can
=#
U_v_tot = merge_gauge(U_v, U_v2);
size(U_v_tot[1])
#=
i.e., from 26 Bloch states to 16 valence MLWFs.
For conduction bands, it's the same
=#
U_c_tot = merge_gauge(U_c, U_c2);
size(U_c_tot[1])
#=
i.e., from 26 Bloch states to 10 conduction MLWFs.

If you prefer to work with the W90's `u.mat` files, you can read them through
[`read_u_mat`](@ref) function, and they should be identical to that from `chk` files
=#
U_v2_test = WannierIO.read_u_mat("$path/TiO2/outputs/val/TiO2_u.mat")[1];
U_v2_test ≈ U_v2
# similarly for conduction
U_c2_test = WannierIO.read_u_mat("$path/TiO2/outputs/cond/TiO2_u.mat")[1];
U_c2_test ≈ U_c2

#=
To get a full gauge matrix for the valence and conduction MLWFs together,
i.e., gauge matrices in the block-diagonal form, we can stack them together, by
=#
U_tot = map(zip(U_v_tot, U_c_tot)) do (v, c)
    [v c]
end;
size(U_tot[1])
# and let's check the spreads--the 1st 16 are valence MLWFs, and the last
# 10 are conduction MLWFs
spread(model, U_tot)
#=
!!! tip

    For convenience, we also provide a function [`Wannier.Tools.merge_gauge`](@ref)
    that merges the valence and conduction `TiO2_split.amn` and `TiO2.chk` files
    into one block-diagonal matrix `TiO2_mrwf.amn`, that maps from the
    valence+conduction `mmn/eig` files directly to the final MRWFs.
=#
outdirs = ["$path/TiO2/val", "$path/TiO2/cond"];
Utot_test = Wannier.Tools.merge_gauge("TiO2", outdirs);
Utot_test ≈ U_tot
# and the file header contains a hint for this
open("TiO2_mrwf.amn") do io
    println(readline(io))
    println(readline(io))
end

#=
In summary:
- `TiO2_split.amn` contains the gauge matrix from valence+conduction
    `mmn`/`eig` to the valence or conduction `mmn`/`eig`
- `TiO2.chk.amn` contains the gauge matrix from valence or conduction
    `mmn`/`eig` to the final MLWFs

If one wants to get the gauge transformation start from the valence+conduction
`mmn`/`eig` files to the final MLWFs of valence or conduction, the two matrices
should be combined together.

Finally, let's compare band structures as sanity check.
=#
w90_v = read_w90_band("$path/TiO2/outputs/val/TiO2", model_v.recip_lattice)[2];
P = plot_band_diff(kpi, w90, w90_v; qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

# and conduction,
w90_c = read_w90_band("$path/TiO2/outputs/cond/TiO2", model_c.recip_lattice)[2];
P = plot_band_diff(kpi, w90, w90_c; qe.fermi_energy)
Main.HTMLPlot(P, 500)  # hide

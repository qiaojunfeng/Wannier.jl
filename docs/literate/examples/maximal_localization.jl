# # Maximal localization of isolated manifold

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In the first tutorial, we run the maximal localization algorithm on
the silicon valence bands to construct WFs for isolated manifold. We need to

1. compute the initial guess, the overlap matrices, and the engery eigenvalues
    from a density function theory (DFT) calculation
    - These files, following wannier90's convention, are called
        `amn`, `mmn`, and `eig`, respectively. In addition, wannier90 defines
        a `win` file containing some input parameters, e.g., the crystal structure.
    - We have pre-computed them using Quantum ESPRESSO (QE) code, and the
        results are stored hosted in
        [WannierDatasets](https://github.com/qiaojunfeng/WannierDatasets) repo.
        For more details, please refer to [WannierDatasets](@ref) section.
2. construct a `Wannier.jl` [`Model`](@ref), by reading the `win`, `amn`, `mmn`,
    and `eig` files
3. run `Wannier.jl` [`max_localize`](@ref) on the `Model` to minimize the spread
4. store the maximal-localized gauge matrices
=#

# ## Preparation
# Load the package
using Wannier
using Wannier.Datasets

#=
## Model generation

We have pre-computed the `amn`/`mmn`/`eig` files using QE, so we can simply
load the datasets by
=#
model = load_dataset("Si2_valence")

#=
Internally, it uses the [`read_w90`](@ref) function to read the
`win`/`amn`/`mmn`/`eig` files, and return a [`Model`](@ref) that stores the
matrices and parameters.

The dataset contains the following files:
=#
show_dataset("Si2_valence")

# Using Julia's [artifact system](https://pkgdocs.julialang.org/v1/artifacts/),
# they are located in the folder
dataset"Si2_valence"

# Therefore, you can equivalently load it by calling the [`read_w90`](@ref) function
read_w90(dataset"Si2_valence/Si2_valence")

#=
!!! note

    Here we are using four bond-centered ``s`` orbitals as initial projections,
    see [Compute bond centers of silicon](@ref) section for a code snippet
    on how to compute the bond center coordinates.
=#

#=
## Maximal localization

Maximal localization can be easily done by calling the
[`max_localize`](@ref) function, which returns the gauge matrices `U`,
=#
U = max_localize(model);

# The initial spread is
omega(model)

# The final spread is
omega(model, U)

#=
Since we have very good initial guess, the spreads only decrease a little bit.

!!! note

    The convergence thresholds is determined by the
    keyword arguments of [`max_localize`](@ref), e.g., `f_tol` for
    the tolerance on spread, and `g_tol` for the tolerance on the
    norm of spread gradient, etc. You can use stricter thresholds
    to further minimize a bit the spread.
=#

#=
## Save the new gauge

You can assign the new gauge to the `Model` by
=#
model.U .= U;

#=
or save the new gauge to an `amn` file, which can be loaded again in
`Wannier.jl`, or used as the a initial guess for wannier90.
=#
write_amn("Si2.maxloc.amn", U)

#=
Voilà! We have finished the first tutorial! 🎉🎉🎉
=#

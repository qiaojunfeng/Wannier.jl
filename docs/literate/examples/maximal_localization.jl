# # Maximal localization of isolated manifold

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In the first tutorial, we maximally localize the spread of silicon valence bands
to construct maximally localized Wannier functions (MLWFs) for an isolated manifold.

In general, a Wannierization of a material consists of the following steps:

1. Compute the initial guess, the overlap matrices, and the engery eigenvalues
    from a density function theory (DFT) calculation
    - These files, following wannier90's convention, are called
        `amn`, `mmn`, and `eig`, respectively. In addition, wannier90 defines
        a `win` file containing some input parameters, e.g., the crystal structure.
    - We have pre-computed them using Quantum ESPRESSO (QE) code, and the
        results are stored hosted in
        [WannierDatasets](https://github.com/qiaojunfeng/WannierDatasets) repo.
        For more details, please refer to [WannierDatasets](@ref) section.
2. Construct a `Wannier.jl` [`Model`](@ref), by reading the `win`, `amn`, `mmn`,
    and `eig` files
3. Run `Wannier.jl` [`max_localize`](@ref) on the `Model` to minimize the spread
4. Store the maximal-localized gauge matrices
=#

# ## Preparation
# Load the package
using Wannier
using Wannier.Datasets
#=
The `Wannier` package contains functions for Wannierization and Wannier
interpolation.
In addition, we also provide a submodule `Wannier.Datasets`, which help us to
load some pre-computed datasets, so that we don't need to run time-consuming
DFT calculations to obtain the `amn`/`mmn`/`eig` files.

## Model construction

We can simply load the datasets by
=#
model = load_dataset("Si2_valence")

#=
Internally, it uses the [`read_w90`](@ref) function to read the
`win`/`amn`/`mmn`/`eig` files, and return a [`Model`](@ref) that stores the
matrices and parameters.

Usually a dataset contains the following files:
=#
show_dataset("Si2")

#=
!!! note

    Here we are peeking at the `Si2` dataset instead of the loaded `Si_valence`
    dataset, since the latter contains wavefunction `UNK` files, which will print
    a large amount of distracting text.
    Of course you can still inspect the loaded dataset files by running
    ```julia
    show_dataset("Si2_valence")
    ```
=#

# Using Julia's [artifact system](https://pkgdocs.julialang.org/v1/artifacts/),
# the `dataset`-string is just a string showing the folder where the files are located,
dataset"Si2_valence"

#=
!!! tip

    Julia's artifact system will automatically download the dataset from the
    [WannierDatasets](https://github.com/qiaojunfeng/WannierDatasets) repo and
    place them under the folder `~/.julia/artifacts/`.

You can equivalently load it by calling the [`read_w90`](@ref) function
=#
read_w90(dataset"Si2_valence/Si2_valence")

#=
Here we need to use `dataset"Si2_valence/Si2_valence"`, where the second
`Si2_valence` is the prefix (or called `seedname` in wannier90) of the files,
i.e., `Si2_valence.amn`/`Si2_valence.mmn` inside the `dataset"Si2_valence"` folder.
You can also use the `dataset`-string to get the actual path of the files, e.g.,
printing the first 5 lines of the `win` file,
=#
open(dataset"Si2_valence/Si2_valence.win") do io
    for i in 1:5
        println(readline(io))
    end
end
#=
!!! tip

    Here we are using four bond-centered ``s`` orbitals as initial projections,
    see [Compute bond centers of silicon](@ref) section for a code snippet
    on how to compute the bond center coordinates.
=#

#=
## Maximal localization

Maximal localization can be easily achieved by calling the
[`max_localize`](@ref) function, which returns the maximally-localized gauge
matrices `U`,
=#
U = max_localize(model);
#=
Here we append a semicolon `;` to suppress the printing of the content of `U`,
which are many complex numbers.

We can compute the initial spread by
=#
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

You can assign the new gauge back to the `Model` by
=#
model.U .= U;
#=
!!! tip

    You need to use `.=` to assign the `U` to the `model`,  because the `model`
    is an immutable Julia `struct`, so it is not allowed to use `model.U = U`.

Or save the new gauge to an `amn` file, which can be loaded again in
`Wannier.jl`, or used as an initial guess in wannier90.
=#
write_amn("Si2.maxloc.amn", U)

#=
Voilà! We have finished the first tutorial! 🎉🎉🎉
=#

# # 1. Maximal localization of isolated manifold

#=
```@meta
CurrentModule = Wannier
```
=#

#=
In the first tutorial, we run the maximal localization algorithm on
the silicon valence bands. We need to

1. generate the `amn`, `mmn`, and `eig` files by using `Quantum ESPRESSO` (QE)
2. construct a [`Model`](@ref) for `Wannier.jl`, by reading the `win`, `amn`, `mmn`, and `eig` files
3. run `Wannier.jl` [`max_localize`](@ref) on the `Model` to minimize the spread
4. write the maximal localized gauge to a new `amn` file

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`1-maxloc.ipynb`](./1-maxloc.ipynb)
    - Julia script: [`1-maxloc.jl`](./1-maxloc.jl)
=#

# ## Preparation
# Load the package
using Wannier
using Printf  # for pretty print

# Path of current tutorial
CUR_DIR = "1-maxloc"

# ### Tips on initial guess
#=
Before we running the calculation, one thing worth mentioning is the
initial projection block in the `win` file for valence bands of silicon.
In such cases, we cannot use the `s` and `p` projections, but using
bond-centered `s` orbitals as guides for the bonding WFs.

The bond centers can be calculated by using the [`find_nearests`](@ref) function
in `Wannier.jl`
=#
win = read_win("$CUR_DIR/si2.win")

#=
!!! tip

    Refer to the [`API`](@ref Input/Output) documentation for more details
    about the function usage.
=#

#=
We will find the 4 nearest atoms to
the 1st atom, i.e. 4 of the periodic images of
the 2nd atom, then calculate the middle point
between each of these 4 atoms and the 1st atom.

!!! note

    In Julia, array index starts from 1, and array is column-major

=#

# fractional coordiantes of the 1st Si atom
atom1 = win.atoms_frac[:, 1]
# cartesian coordiantes of `atom1`
atom1_cart = win.unit_cell * atom1
# find 4 nearest neighbors of `atom1`,
# I use 5 because the 1st returned neighbor is the `atom1` itself
distances, indexes, translations = Wannier.find_nearests(
    atom1, 5, win.unit_cell, win.atoms_frac
)
# print the nearest atom and bond center, in cartesian coordinates
for i in 2:5
    idx = indexes[i]
    nn_frac = win.atoms_frac[:, idx] + translations[:, i]
    nn_cart = win.unit_cell * nn_frac
    @printf(
        "nearest atom: \n  index = %d, distance = %.5f, (x, y, z) = (%.5f, %.5f, %.5f)\n",
        idx,
        distances[i],
        nn_cart...
    )
    c = (atom1_cart + nn_cart) / 2
    @printf("  bond center = (%.5f, %.5f, %.5f)\n", c...)
end

#=
Now we can use these outputs for writing `projection` block in `win` file,
and run the QE calculations for `amn`, `mmn`, and `eig` files.

!!! tip

    Use the `run.sh` script which automate the scf, nscf, pw2wannier90 steps.
=#

#=
## Model generation

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct a [`Model`](@ref Model)
that abstracts the calculation
=#
model = read_w90("$CUR_DIR/si2")

#=
## Maximal localization

Maximal localization can be easily done by calling the
[`max_localize`](@ref) function, which returns the gauge matrices `U`,
=#
U = max_localize(model)

# The initial spread is
omega(model)

# The final spread is
omega(model, U)

#=
Since we have very good initial guess, the
spread only decreases a little bit.

!!! note

    The convergence thresholds is determined by the
    keyword arguments of [`max_localize`](@ref), e.g., `f_tol` for
    the tolerance on spread, and `g_tol` for the tolerance on the
    norm of spread gradient, etc. You can use stricter thresholds
    to further minimize a bit the spread.
=#

#=
## Save the new gauge

We can further do band interpolation in `Wannier.jl` with
the new `U` matrices, however, for this very first tutorial,
we will just save the new gauge to an `amn` file,
which can be used as the new initial guess for `Wannier90`,
or reuse it in `Wannier.jl`.
=#
write_amn("$CUR_DIR/si2.maxloc.amn", U)

#=
Voilà! We have finished the first tutorial! 🎉🎉🎉
=#

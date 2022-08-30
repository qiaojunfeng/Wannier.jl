# # 1. Maximal localization of isolated manifold

#=
In the first tutorial, we run the maximal localization algorithm on
the silicon valence bands. We need to

1. generate the `amn`, `mmn`, and `eig` files by using `Quantum ESPRESSO` (QE)
2. construct a [`Model`](@ref Model) for `Wannier.jl`, by reading the `win`, `amn`, `mmn`, and `eig` files
3. run `Wannier.jl` [`max_localize`](@ref max_localize) on the `Model` to minimize the spread
4. write the maximal localized gauge to a new `amn` file

!!! note

    These tutorials assume you have already been familiar with the
    Wannierization workflow using `QE` and `Wannier90`, a good starting
    point could be the tutorials of
    [`Wannier90`](https://github.com/wannier-developers/wannier90).
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

The bond centers can be calculated by using the [`find_nearests`](@ref find_nearests) function
in `Wannier.jl`
=#
win = read_win("$CUR_DIR/si2.win")

#=
!!! tip

    Refer to the [`API`](@ref) documentation for more details about the function usage.
=#

#=
We will find the 4 nearest atoms to
the 2nd atom, i.e. 4 of the periodic images of
the 1st atom, then calculate the middle point
between these 4 atoms and the 2nd atom.

!!! note

    In Julia, array index starts from 1, and array is column-major

=#

# fractional coordiantes of the 2nd Si atom
atom2 = win.atoms_frac[:, 2]
# cartesian coordiantes of `atom2`
atom2_cart = win.unit_cell * atom2
# find 4 nearest neighbors of `atom2`,
# I use 5 because the 1st returned neighbor is the `atom2` itself
distances, indexes, translations = Wannier.find_nearests(
    atom2, 5, win.unit_cell, win.atoms_frac
);
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
    c = (atom2_cart + nn_cart) / 2
    @printf("  bond center = (%.5f, %.5f, %.5f)\n", c...)
end

#=
Now we can use these outputs for writing `projection` block in `win` file, and run the QE calculations for `amn`, `mmn`, and `eig` files.

!!! tip

    Use the `run.sh` script which automate the scf, nscf, pw2wannier90 steps.
=#

#=
## Model generation

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct a [`Model`](@ref Model) that abstracts the calculation
=#
model = read_w90("$CUR_DIR/si2")

#=
## Maximal localization

Maximal localization can be easily done by calling the [`max_localize`](@ref max_localize) function, which
returns the gauge matrices `A`
=#
A = max_localize(model)

#=
## Save the new gauge

We can further do band interpolation in `Wannier.jl` with
the new `A` matrices, however, for this very first tutorial,
we will just save the new gauge to an `amn` file,
which can be used as the new initial guess for `Wannier90`,
or reuse it in `Wannier.jl`.
=#
write_amn("$CUR_DIR/si2.maxloc.amn", A)

#=
Voilà! We have finished the first tutorial! 🎉🎉🎉
=#

# # Tips and code snippets

#=
Here are some small code snippets that might be useful.
=#

# ## Compute bond centers of silicon

#=
When Wannierizing the valence manifold of silicon, a good initial guess is to
put four ``s`` orbitals at the four bond centers.

The bond centers can be calculated by using the [`find_nearests`](@ref) function:
=#
using Printf
using Wannier
using Wannier.Datasets

win = read_win(dataset"Si2/Si2.win")

#=
We will find the 4 nearest atoms to the 1st atom, i.e. 4 of the periodic images of
the 2nd atom, then compute the middle point between each of these 4 atoms and the 1st atom.

!!! note

    In Julia, array index starts from 1, and array is column-major
=#

# the `atoms_frac` tag in the `win` file is a nested data structure
win.atoms_frac

# convert it to a Vector of atomic labels
atom_labels = [first(x) for x in win.atoms_frac]
# and a Vector of fractional coordinates
atom_positions = [last(x) for x in win.atoms_frac]

# the fractional coordinates of the 1st Si atom is
atom1 = atom_positions[1]
# convert to Cartesian coordinates in Å
lattice = win.unit_cell_cart
atom1_cart = lattice * atom1

# find 4 nearest neighbors of `atom1`,
# I use 5 because the 1st returned neighbor is the `atom1` itself
distances, indices, translations = Wannier.find_neighbors(atom1, 5, lattice, atom_positions)
# print the nearest atom and bond center, in Cartesian coordinates
for i in 2:5
    idx = indices[i]
    nn_frac = atom_positions[idx] + translations[i]
    nn_cart = lattice * nn_frac
    @printf(
        "nearest atom:\n  index = %d, distance = %.5f Å, (x, y, z) = (%.5f, %.5f, %.5f) Å\n",
        idx,
        distances[i],
        nn_cart...
    )
    c = (atom1_cart + nn_cart) / 2
    @printf("  bond center = (%.5f, %.5f, %.5f) Å\n", c...)
end

#=
Now we can use these centers for writing `projection` block in `win` file,
and run the QE pw2wannier90 calculation for computing `amn`, `mmn`, and `eig` files.
=#

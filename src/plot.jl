export bandplot, bandplot!

"""
    bandplot(kpath, eigenvals; kwargs...)
    bandplot(x, symm_indices, symm_labels, eigenvals; kwargs...)

Plot the band structure given a path of k-points and the corresponding eigenvalues.
The kpath can be specified in three ways:
1. As a `RecipPath` object
2. As a vector of cumulative distances along the kpath (in units of 1/L, where L is unit of length),
    together with the indices and labels of high-symmetry kpoints
"""
function bandplot end

function bandplot! end

export n_atoms, n_kpoints, n_bvectors, n_bands, n_wannier, n_Rvectors

# This file list the function interfaces and docstrings for getting the
# size(s) of some data structures.

"""
    n_atoms(::Model)

Number of atoms.
"""
function n_atoms end

"""
    n_kpoints(::Model)
    n_kpoints(::KspaceStencil)
    n_kpoints(::KspaceStencilShells)

Number of kpoints.
"""
function n_kpoints end

"""
    n_bvectors(::Model)
    n_bvectors(::KspaceStencil)
    n_bvectors(::KspaceStencilShells)

Number of ``\\mathbf{b}``-vectors.
"""
function n_bvectors end

"""
    n_bands(::Model)

Number of bands.
"""
function n_bands end

"""
    n_wannier(::Model)

Number of Wannier functions.
"""
function n_wannier end

"""
    n_Rvectors(::AbstractRspace)
    n_Rvectors(::TBOperator)

Number of ``\\mathbf{R}``-vectors.
"""
function n_Rvectors end

export read_amn_ortho

"""
    $(SIGNATURES)

Read `amn` file and Lowdin orthonormalize the unitary matrices.

The `U` matrix for Wannier functions must be unitary or semi-unitary.
Thus, in most cases, this function should be used instead of `WannierIO.read_amn`,
where the latter one just parse the `amn` file and return whatever is in it.
"""
function read_amn_ortho(filename::AbstractString)
    U = read_amn(filename)
    U .= orthonorm_lowdin(U)
    @info "Lowdin orthonormalization applied to U matrices"
    return U
end

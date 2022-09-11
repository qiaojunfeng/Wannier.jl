export read_orthonorm_amn

"""
    read_orthonorm_amn(filename::AbstractString)

Read and orthonormalize the `amn` file.

Wrapper function to read `amn` and Lowdin orthonormalize it.
The `A` matrix needs to be unitary or semi-unitary,
so in most cases this function should be used instead of [`WannierIO.read_amn`](@ref).

See also [`WannierIO.read_amn`](@ref).
"""
function read_orthonorm_amn(filename::AbstractString)
    A = read_amn(filename)
    A .= orthonorm_lowdin(A)
    return A
end

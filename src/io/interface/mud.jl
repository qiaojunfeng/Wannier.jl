
@doc raw"""
    compute_mud(dir_up, dir_dn)

Compute the overlap matrix between spin up and down from UNK files.

```math
M^{\uparrow\downarrow}_{m n \bm{k}} =
\langle u^{\uparrow}_{m \bm{k}} | u^{\downarrow}_{n \bm{k}} \rangle
```

This function compute `mud` matrix in real space (thus much slower),
checked with QE `pw2wannier90.x` function which computes `mud` in reciprocal
space, the difference is on the order of 1e-13.

!!! warning
    This only works for norm-conserving pseudopotentials since in that case
    the overlap operator is just identity; for ultrasoft or PAW, a simple
    dot product is not enough.
    Also I assume the UNK files are written with a normalization factor of
    ``N`` (total points of the FFT grid) over the unit cell, i.e.,
    the UNK generated by QE `pw2wannier90.x`. If the UNK files are normalized
    to 1, the result should be multiplied by ``N``.

# Arguments
- `dir_up`: directory of spin up UNK files
- `dir_dn`: directory of spin down UNK files
"""
function compute_mud(dir_up::AbstractString, dir_dn::AbstractString)
    unk_up = filter(x -> startswith(x, "UNK"), readdir(dir_up))
    unk_up = sort(unk_up; by=x -> parse(Int, x[4:(end - 2)]))

    unk_dn = filter(x -> startswith(x, "UNK"), readdir(dir_dn))
    unk_dn = sort(unk_dn; by=x -> parse(Int, x[4:(end - 2)]))

    @assert length(unk_up) == length(unk_dn)
    n_kpts = length(unk_up)

    ik, ψk = read_unk("up/$(unk_up[1])")
    n_bands = size(ψk, 4)
    N = prod(size(ψk)[1:3])

    Mud = zeros(ComplexF64, n_bands, n_bands, n_kpts)

    for ik in 1:n_kpts
        _, ψu = read_unk("$(dir_up)/$(unk_up[ik])")
        _, ψd = read_unk("$(dir_dn)/$(unk_dn[ik])")
        for i in 1:n_bands
            for j in 1:n_bands
                Mud[i, j, ik] = sum(conj.(ψu[:, :, :, i]) .* ψd[:, :, :, j])
            end
        end
    end
    Mud ./= N
    return Mud
end

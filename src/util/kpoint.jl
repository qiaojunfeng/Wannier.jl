
"""
kpoints: in fractional coordinates, 3 * n_kpts
"""
function get_kpoint_mappings(kpoints::Matrix{T}, kgrid::AbstractVector{Int}) where {T<:Real}
    n_kpts = prod(kgrid)
    n_kx, n_ky, n_kz = kgrid
    dkx, dky, dkz = 1 / n_kx, 1 / n_ky, 1 / n_kz

    kpts_int = kpoints ./ [dkx; dky; dkz]
    kpts_int = round.(Int, kpts_int)

    for ik in 1:n_kpts
        kpts_int[1, ik] = mod(kpts_int[1, ik], 0:(n_kx - 1)) + 1
        kpts_int[2, ik] = mod(kpts_int[2, ik], 0:(n_ky - 1)) + 1
        kpts_int[3, ik] = mod(kpts_int[3, ik], 0:(n_kz - 1)) + 1
    end

    k_xyz = Vector{Vec3{Int}}(undef, n_kpts)
    xyz_k = Array{Int,3}(undef, n_kx, n_ky, n_kz)

    for ik in 1:n_kpts
        k_xyz[ik] = kpts_int[:, ik]
        xyz_k[k_xyz[ik]...] = ik
    end

    return k_xyz, xyz_k
end

@doc raw"""
Make a supercell of kpoints by translating it along 3 directions.
Input and returned kpoints are in fractional coordinates.

repeat: number of repetitions along ±x, ±y, ±z directions, on output
there are (2*repeat + 1)^3 cells.
"""
function make_supercell(kpoints::Matrix{T}, replica::R=5) where {T<:Real,R<:Integer}
    size(kpoints, 1) ≉ 3 && error("kpoints must be 3 * n_kpts")
    n_kpts = size(kpoints, 2)

    n_cell = (2 * replica + 1)^3

    supercell = Matrix{T}(undef, 3, n_cell * n_kpts)
    translations = Matrix{Int}(undef, 3, n_cell * n_kpts)

    counter = 1
    for ix in (-replica):replica
        for iy in (-replica):replica
            for iz in (-replica):replica
                for ik in 1:n_kpts
                    supercell[:, counter] = kpoints[:, ik] + [ix, iy, iz]
                    translations[:, counter] = [ix, iy, iz]
                    counter += 1
                end
            end
        end
    end

    return supercell, translations
end

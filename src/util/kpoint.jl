
"""
kpoints: in fractional coordinates, 3 * n_kpts
"""
function get_kpoint_mappings(kpoints::Matrix{T}, kgrid::AbstractVector{Int}) where {T<:Real}
    n_kpts = prod(kgrid)
    n_kx, n_ky, n_kz = kgrid
    dkx, dky, dkz = 1 / n_kx, 1 / n_ky, 1 / n_kz

    kpts_int = kpoints ./ [dkx; dky; dkz]
    kpts_int = round.(Int, kpts_int)

    for ik = 1:n_kpts
        kpts_int[1, ik] = mod(kpts_int[1, ik], 0:n_kx-1) + 1
        kpts_int[2, ik] = mod(kpts_int[2, ik], 0:n_ky-1) + 1
        kpts_int[3, ik] = mod(kpts_int[3, ik], 0:n_kz-1) + 1
    end

    k_xyz = Vector{Vec3{Int}}(undef, n_kpts)
    xyz_k = Array{Int,3}(undef, n_kx, n_ky, n_kz)

    for ik = 1:n_kpts
        k_xyz[ik] = kpts_int[:, ik]
        xyz_k[k_xyz[ik]...] = ik
    end

    k_xyz, xyz_k
end

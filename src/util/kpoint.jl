
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
function make_supercell(
    kpoints::Matrix{T}, replica::AbstractVector{R}
) where {T<:Number,R<:AbstractRange}
    size(kpoints, 1) ≉ 3 && error("kpoints must be 3 * n_kpts")
    n_kpts = size(kpoints, 2)

    rep_x, rep_y, rep_z = replica
    n_cell = length(rep_x) * length(rep_y) * length(rep_z)

    supercell = Matrix{T}(undef, 3, n_cell * n_kpts)
    translations = Matrix{Int}(undef, 3, n_cell * n_kpts)

    counter = 1
    for ix in rep_x
        for iy in rep_y
            for iz in rep_z
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

function make_supercell(kpoints::Matrix{T}, replica::R=5) where {T<:Number,R<:Integer}
    return make_supercell(
        kpoints, [(-replica):replica, (-replica):replica, (-replica):replica]
    )
end

@doc raw"""
Generate list of kpoint coordinates from kpoint grid.
Work just like `kmesh.pl` of W90.

kgrid: contains a N1 * N2 * N3 mesh
return an explicit list of kpoints in fractional coordinates
"""
function get_kpoints(kgrid::AbstractVector{Int})
    n_pts = prod(kgrid)
    nx, ny, nz = kgrid

    kpoints = zeros(Float64, 3, n_pts)
    idx = 1
    for x in 0:(nx - 1)
        for y in 0:(ny - 1)
            for z in 0:(nz - 1)
                kpoints[:, idx] = [x / nx, y / ny, z / nz]
                idx += 1
            end
        end
    end

    return kpoints
end

"""
Sort kpoints such that z increases the fastest, then y, then x.
"""
function sort_kpoints(kpoints::AbstractMatrix{T}) where {T<:Real}
    return sortslices(kpoints; dims=2)
end

@doc raw"""
Guess kgrid from list of kpoint coordinates.

kpoints: contains a N1 * N2 * N3 kpoints, in fractional coordinates
"""
function get_kgrid(kpoints::AbstractMatrix{T}) where {T<:Real}
    kgrid = zeros(Int, 3)
    # 3 directions
    for i in 1:3
        uniq_kpt = unique(kpoints[i, :])
        kgrid[i] = length(uniq_kpt)
    end

    if prod(kgrid) != size(kpoints, 2)
        error("kgrid and kpoints do not match")
    end

    kpoints_recovered = get_kpoints(kgrid)
    kpoints_sorted = sort_kpoints(kpoints)
    if !all(isapprox.(kpoints_recovered, kpoints_sorted; atol=1e-5))
        error("cannot convert kpoints to a kgrid")
    end

    return kgrid
end


"""
    get_kpoint_mappings(kpoints, kgrid)

Get the mappings between kpoint indexes and kpoint coordiantes.

Return a tuple of `(k_xyz, xyz_k)`:
- `k_xyz[ik]` maps kpoint `kpoints[:, ik]` to kpoint coordinates `[ikx, iky, ikz]`
- `xyz_k[ikx, iky, ikz]` maps kpoint coordinates `[ikx, iky, ikz]` to kpoint index `ik`
- the kpoint fractional coordinates is `[(ikx - 1)/nkx, (iky - 1)/nky, (ikz - 1)/nkz]`

# Arguments
- `kpoints`: `3 * n_kpts`, in fractional coordinates
- `kgrid`: `3`, number of kpoints along each reciprocal lattice vector
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

"""
    make_supercell(kpoints, replica)

Make a supercell of kpoints by translating it along 3 directions.

On output there are `(2*replica + 1)^3` cells, in fractional coordinates.

# Arguments
- `kpoints`: `3 * n_kpts`, in fractional coordinates
- `replica`: `3`, number of repetitions along ±x, ±y, ±z directions
"""
function make_supercell(
    kpoints::AbstractMatrix{<:Number}, replica::AbstractVector{<:AbstractRange}
)
    size(kpoints, 1) ≉ 3 && error("kpoints must be 3 * n_kpts")
    n_kpts = size(kpoints, 2)

    rep_x, rep_y, rep_z = replica
    n_cell = length(rep_x) * length(rep_y) * length(rep_z)

    supercell = similar(kpoints, 3, n_cell * n_kpts)
    translations = zeros(Int, 3, n_cell * n_kpts)

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

"""
    make_supercell(kpoints, replica=5)

Make a supercell of kpoints by translating it along 3 directions.

# Arguments
- `replica`: integer, number of repetitions along ±x, ±y, ±z directions
"""
function make_supercell(kpoints::AbstractMatrix{<:Number}, replica::Integer=5)
    return make_supercell(
        kpoints, [(-replica):replica, (-replica):replica, (-replica):replica]
    )
end

"""
    get_kpoints(kgrid; fractional=true, endpoint=false)

Generate list of kpoint coordinates from kpoint grid.

# Arguments
- `kgrid`: vector of 3 integers specifying a `nkx * nky * nkz` mesh

# Keyword Arguments
- `fractional`: return an explicit list of kpoints in fractional coordinates, else integers
- `endpoint`: include the endpoint of the grid, only for fractional case. E.g., if true, 1.0 is included

!!! note

    If the default keyword arguments are used, this function works just like `kmesh.pl` of `Wannier90`.
"""
function get_kpoints(
    kgrid::AbstractVector{<:Integer}; fractional::Bool=true, endpoint::Bool=false
)
    endpoint && !fractional && error("endpoint can only be used for fractional coordinates")

    n_pts = prod(kgrid)
    nx, ny, nz = kgrid

    if fractional
        kpoints = zeros(Float64, 3, n_pts)
    else
        kpoints = zeros(Int, 3, n_pts)
    end

    idx = 1
    for x in 0:(nx - 1)
        for y in 0:(ny - 1)
            for z in 0:(nz - 1)
                kpoints[:, idx] = [x, y, z]
                idx += 1
            end
        end
    end

    if fractional
        if endpoint
            kpoints[1, :] ./= (nx - 1)
            kpoints[2, :] ./= (ny - 1)
            kpoints[3, :] ./= (nz - 1)
        else
            kpoints[1, :] ./= nx
            kpoints[2, :] ./= ny
            kpoints[3, :] ./= nz
        end
    end

    return kpoints
end

"""
    sort_kpoints(kpoints)

Sort kpoints such that z increases the fastest, then y, then x.

# Arguments
- `kpoints`: `3 * n_kpts`
"""
function sort_kpoints(kpoints::AbstractMatrix{T}) where {T<:Real}
    return sortslices(kpoints; dims=2)
end

"""
    get_kgrid(kpoints)

Guess kgrid from list of kpoint coordinates.

Input `kpoints` has size `3 * n_kpts`, where `n_kpts = nkx * nky *  nkz`,
output `[nkx, nky, nkz]`.

# Arguments
- `kpoints`: `3 * n_kpts`, fractional coordiantes
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

    if T <: Integer
        kpoints_recovered = get_kpoints(kgrid; fractional=false)
    else
        kpoints_recovered = get_kpoints(kgrid)
    end
    kpoints_sorted = sort_kpoints(kpoints)
    if !all(isapprox.(kpoints_recovered, kpoints_sorted; atol=1e-5))
        error("cannot convert kpoints to a kgrid")
    end

    return kgrid
end

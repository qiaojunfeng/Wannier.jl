"""
    $(SIGNATURES)

Get the mappings between kpoint indices and kpoint coordiantes.

# Arguments
- `kpoints`: length-`n_kpts` vector of fractional coordinates, e.g. `Vec3`.
    The kpoint fractional coordinates should be
    `[(ikx - 1)/n_kx, (iky - 1)/n_ky, (ikz - 1)/n_kz]` where `n_kx, n_ky, n_kz = kgrid_size`
    are the number of kpoints along each reciprocal lattice vector, and
    `ikx`, `iky`, `ikz` are integers in `1:n_kx`, `1:n_ky`, `1:n_kz`, respectively.
- `kgrid_size`: length-`3` vector, number of kpoints along each reciprocal lattice vector

# Return
- `k_xyz`: mapping from linear index to Cartesian index, such that
    `ix, iy, iz = k_xyz[ik]` maps kpoint `kpoints[ik]` to kpoint fractional coordinate
    `[ix, iy, iz] ./ kgrid`
- `xyz_k`: `xyz_k[ikx, iky, ikz]` maps kpoint coordinates `[ikx, iky, ikz]` to kpoint index `ik`
"""
function get_kpoint_mappings(kpoints::AbstractVector, kgrid_size::AbstractVector)
    n_kpts = prod(kgrid_size)
    n_kx, n_ky, n_kz = kgrid_size
    dk = Vec3(1 / n_kx, 1 / n_ky, 1 / n_kz)
    kgrid = Vec3(kgrid_size...)

    kpoints_int = map(kpoints) do k
        t = round.(Int, k ./ dk)
        # I assume the fractional kpoints are in [0, 1), so plus 1 to make sure
        # the integer indices are in 1:n_k
        t = mod.(t, kgrid) .+ 1
        return t
    end

    k_xyz = zeros(Vec3{Int}, n_kpts)
    xyz_k = zeros(Int, n_kx, n_ky, n_kz)

    for ik in 1:n_kpts
        k_xyz[ik] = Vec3(kpoints_int[ik]...)
        xyz_k[k_xyz[ik]...] = ik
    end

    return k_xyz, xyz_k
end

"""
    $(SIGNATURES)

Make a supercell of points by translating it along 3 directions.

# Arguments
- `points`: length-`n_pts` vector, each element is a length-`3` vector of
    coordinates, e.g. `Vec3`
- `replica`: length-`3` vector of range, or an integer; number of repetitions
    along ±x, ±y, ±z directions

# Return
- `supercell`: a length-`n_cells * n_pts` vector of `Vec3` for each point in the
    supercell, where `n_cells` is the number of replica and `n_pts` is the number
    of points in `points`. e.g., if `replica` is a number, then there are
    `(2*replica + 1)^3` cells.
- `translations`: the translation vectors (w.r.t. the original point in `points`)
    for each point in the `supercell`
"""
function make_supercell(
    points::AbstractVector, replica::AbstractVector{R}
) where {R<:AbstractRange}
    n_pts = length(points)
    @assert n_pts > 0 "points is empty"
    @assert length(replica) == 3 "replica must be length-3 vector"

    rep_x, rep_y, rep_z = replica
    n_cells = length(rep_x) * length(rep_y) * length(rep_z)

    T = eltype(points[1])
    supercell = zeros(Vec3{T}, n_cells * n_pts)
    translations = zeros(Vec3{Int}, n_cells * n_pts)

    counter = 1
    for ix in rep_x
        for iy in rep_y
            for iz in rep_z
                for i in 1:n_pts
                    supercell[counter] = Vec3(points[i]) + Vec3(ix, iy, iz)
                    translations[counter] = Vec3(ix, iy, iz)
                    counter += 1
                end
            end
        end
    end

    return supercell, translations
end

function make_supercell(points::AbstractVector, replica::Integer=5)
    return make_supercell(
        points, [(-replica):replica, (-replica):replica, (-replica):replica]
    )
end

"""
    $(SIGNATURES)

Get the kpoint integer indices from a kpoint grid.

# Arguments
- `kgrid_size`: length-`3` vector, number of kpoints along each reciprocal lattice vector
- `contain_zero`: whether to include 0 or 1 in the returned indices

# Return
- `kpoint_indices`: length-`n_kpts` vector of `Vec3{Int}` for each kpoint index
    in the kpoint grid. ``k_z`` increases the fastest, then ``k_y``, then ``k_x``.

See also [`get_kpoints`](@ref).
"""
function get_kpoint_indices(kgrid_size::AbstractVector; contain_zero::Bool=false)
    n_kx, n_ky, n_kz = kgrid_size
    n_kpts = prod(kgrid_size)

    kpoint_indices = zeros(Vec3{Int}, n_kpts)
    counter = 1
    if contain_zero
        x_range = 0:(n_kx - 1)
        y_range = 0:(n_ky - 1)
        z_range = 0:(n_kz - 1)
    else
        x_range = 1:n_kx
        y_range = 1:n_ky
        z_range = 1:n_kz
    end
    for ikx in x_range
        for iky in y_range
            for ikz in z_range
                kpoint_indices[counter] = Vec3(ikx, iky, ikz)
                counter += 1
            end
        end
    end

    return kpoint_indices
end

"""
    $(SIGNATURES)

Generate list of kpoint coordinates from kpoint grid.

# Arguments
- `kgrid_size`: length-`3` vector, number of kpoints along each reciprocal lattice vector
- `endpoint`: whether to the last point is 1.0 or not, e.g. `[0.0, 1.0]` or
    `[0.0, 0.5]`

See also [`get_kpoint_indices`](@ref).

!!! note

    If the default keyword arguments are used, this function works just like `kmesh.pl` of wannier90.
"""
function get_kpoints(kgrid_size::AbstractVector; endpoint::Bool=false)
    any(isone, kgrid_size) &&
        endpoint &&
        error("cannot have endpoint when kgrid contains 1")

    kpoint_indices = get_kpoint_indices(kgrid_size; contain_zero=true)
    denom = Vec3(kgrid_size)
    if endpoint
        denom = Vec3(kgrid_size) .- 1
    end

    return map(kpoint_indices) do k
        k ./ denom
    end
end

"""
    $(SIGNATURES)

Sort points such that z increases the fastest, then y, then x.

# Arguments
- `points`: (often kpoints) length-`n_kpts` vector of fractional coordinates,
    e.g. `Vec3`
"""
function sort_points(points::AbstractVector)
    return sort(points)
end

"""
    $(SIGNATURES)

Guess kgrid_size from list of kpoint coordinates.

Input `kpoints` has size `3 * n_kpts`, where `n_kpts = nkx * nky *  nkz`,
output `[nkx, nky, nkz]`.

# Arguments
- `kpoints`: length-`n_kpts`, fractional coordiantes

# Keyword arguments
- `atol`: absolute tolerance for comparing kpoint coordinates
"""
function guess_kgrid_size(kpoints::AbstractVector; atol=1e-5)
    @assert length(kpoints) > 0 "kpoints is empty"

    kgrid_size = zeros(Int, 3)
    # 3 directions
    for i in 1:3
        uniq_kpt = unique(k -> k[i], kpoints)
        kgrid_size[i] = length(uniq_kpt)
    end

    if prod(kgrid_size) != length(kpoints)
        error("kgrid_size and kpoints do not match")
    end

    # do some sanity check to be very safe
    if eltype(kpoints[1]) <: Integer
        kpoints_recovered = get_kpoint_indices(kgrid_size; contain_zero=true)
    else
        kpoints_recovered = get_kpoints(kgrid_size; endpoint=false)
    end
    kpoints_sorted = sort_points(kpoints)
    kpoints_recovered_sorted = sort_points(kpoints_recovered)

    # if they are shifted by a constant, then our guess is fine
    diff = kpoints_sorted .- kpoints_recovered_sorted
    diff .-= Ref(diff[1])
    if !all(isapprox.(diff, Ref([0, 0, 0]); atol))
        error("cannot guess kgrid_size from kpoint coordinates")
    end

    return kgrid_size
end

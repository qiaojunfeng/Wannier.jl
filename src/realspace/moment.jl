export moment, center, omega, position_op

"""
    moment(rgrid::RGrid, W::AbstractArray, n)

Compute WF moment to arbitrary order in real space.

# Arguments
- `rgrid`: real space grid on which `W` is defined
- `W`: WFs, `nx * ny * nz * n_wann`, or `nx * ny * nz` for single WF
- `n`: order of moment, e.g., 1 for WF center, 2 for variance, etc.

Returned value in Cartesian coordinates.

!!! note

    The WFs are defined in a supercell that is `n_kpts` times unit cell,
    however, usually we only calculate realspace WFs in a smaller supercell
    that is 2^3 or 3^3 times unit cell (as defined by the `n_supercells` argument
    of [`read_realspace_wf`](@ref read_realspace_wf)).
    Some times this is not sufficient if the WFs are
    truncated by the boundries of the smaller supercell, thus the center
    calculated by this function is inexact. In principle, we should calculate
    centers in the `n_kpts` supercell, however, this is memory-consuming.
"""
function moment(rgrid::RGrid, W::AbstractArray{T,3}, n::U) where {T<:Complex,U<:Integer}
    Xᶜ, Yᶜ, Zᶜ = cartesianize_xyz(rgrid)
    x = sum(conj(W) .* Xᶜ .^ n .* W)
    y = sum(conj(W) .* Yᶜ .^ n .* W)
    z = sum(conj(W) .* Zᶜ .^ n .* W)
    r = [x, y, z]
    return real(r)
end

function moment(rgrid::RGrid, W::AbstractArray{T,4}, n::U) where {T<:Complex,U<:Integer}
    n_wann = size(W, 4)
    r = Matrix{real(T)}(undef, 3, n_wann)
    for i in 1:n_wann
        r[:, i] = moment(rgrid, W[:, :, :, i], n)
    end
    return r
end

"""
    center(rgrid::RGrid, W::AbstractArray)

Compute WF center in real space.

Returned value in Cartesian coordinates.

See also [`moment`](@ref moment).
"""
center(rgrid::RGrid, W::AbstractArray) = moment(rgrid, W, 1)

"""
    omega(rgrid::RGrid, W::AbstractArray)

Compute WF spread in real space.

Returned value in Å^2 unit.

See also [`moment`](@ref moment).
"""
function omega(rgrid::RGrid, W::AbstractArray)
    return sum(moment(rgrid, W, 2) - center(rgrid, W) .^ 2; dims=1)
end

"""
    position_op(rgrid::RGrid, W::AbstractArray{T,4})

Compute position operator matrices in real space.

Returned value in Cartesian coordinates.

See also [`center`](@ref center).
"""
function position_op(rgrid::RGrid, W::AbstractArray{T,4}) where {T<:Complex}
    Xᶜ, Yᶜ, Zᶜ = cartesianize_xyz(rgrid)
    n_wann = size(W, 4)
    # last index is x,y,z
    r = zeros(T, n_wann, n_wann, 3)
    for i in 1:n_wann
        for j in 1:n_wann
            Wᵢ = W[:, :, :, i]
            Wⱼ = W[:, :, :, j]
            r[i, j, 1] = sum(conj(Wᵢ) .* Xᶜ .* Wⱼ)
            r[i, j, 2] = sum(conj(Wᵢ) .* Yᶜ .* Wⱼ)
            r[i, j, 3] = sum(conj(Wᵢ) .* Zᶜ .* Wⱼ)
        end
    end
    return r
end

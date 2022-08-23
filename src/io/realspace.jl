using Printf: @sprintf
using LazyGrids: ndgrid

struct RGrid{T<:Real,XT<:AbstractArray3,YT<:AbstractArray3,ZT<:AbstractArray3}
    # spanning vectors, 3 * 3, each column is a spanning vector
    basis::Mat3{T}

    # usually these are LazyGrids
    # x fractional coordinate, nx * ny * nz
    X::XT
    # y fractional coordinate, nx * ny * nz
    Y::YT
    # z fractional coordinate, nx * ny * nz
    Z::ZT
end

function RGrid(basis::AbstractMatrix, X, Y, Z)
    size(X) == size(Y) == size(Z) || error("X, Y, Z must have the same size")
    return RGrid(Mat3(basis), X, Y, Z)
end

@doc """
Get origin of the RGrid.
"""
function origin(rgrid::RGrid)
    O = [rgrid.X[1, 1, 1], rgrid.Y[1, 1, 1], rgrid.Z[1, 1, 1]]
    origin = rgrid.lattice * O
    return origin
end

@doc """
Get the span vectors of the RGrid.

Assumptions:
1. the grid is uniformly spaced
2. the last point along each direction is NOT the periodic point of the 1st point
"""
function span_vectors(rgrid::RGrid)
    O = [rgrid.X[1, 1, 1], rgrid.Y[1, 1, 1], rgrid.Z[1, 1, 1]]
    v1 = [rgrid.X[end, 1, 1], rgrid.Y[end, 1, 1], rgrid.Z[end, 1, 1]] - O
    v2 = [rgrid.X[1, end, 1], rgrid.Y[1, end, 1], rgrid.Z[1, end, 1]] - O
    v3 = [rgrid.X[1, 1, end], rgrid.Y[1, 1, end], rgrid.Z[1, 1, end]] - O
    # the last point is NOT the periodic repetation of the 1st point,
    # so the spanning vector includes another displacement
    nx, ny, nz = size(rgrid.X)
    v1 .*= (nx + 1) / nx
    v2 .*= (ny + 1) / ny
    v3 .*= (nz + 1) / nz
    # to cartesian
    v1 = rgrid.basis * v1
    v2 = rgrid.basis * v2
    v3 = rgrid.basis * v3
    # each column is a vector
    spanvec = hcat(v1, v2, v3)
    return spanvec
end

@doc """
Get the extension name of UNK files, e.g. `NC` from `UNK00001.NC`.
e.g. `1` for no-spin calcluation, `NC` for non-collinear calculation.
"""
function _get_unk_ext(unkdir::AbstractString)
    files = readdir(unkdir)
    i = findfirst(x -> startswith(x, "UNK"), files)
    isnothing(i) && error("No UNK files found in $unkdir")
    suffix = split(files[i], ".")[2]  # "1" or "NC"
    return suffix
end

"""
Read unk files, and generate realspace Wannier functions.

A: n_bands x n_wann x n_kpts, rotation matrix
kpoints: 3 x n_kpts, fractional coordinates
R: fractional coordinates w.r.t lattice (actually integers), the cell for WF |nR>,
    default generate WF at home unit cell |n0>
"""
function read_realspace_wf(
    A::AbstractArray{Complex{T},3},
    kpoints::AbstractMatrix{T},
    n_supercells::AbstractVector{Int},
    unkdir::String=".";
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    n_bands, n_wann, n_kpts = size(A)
    length(n_supercells) != 3 && error("n_supercells must be 3-vector")

    supercells = Vector{UnitRange}()
    for i in 1:3
        n_s = n_supercells[i]
        d, r = divrem(n_s, 2)
        if r == 0
            push!(supercells, (-d):(d - 1))
        else
            push!(supercells, (-d):d)
        end
    end
    n_sx, n_sy, n_sz = n_supercells

    ext = _get_unk_ext(unkdir)
    unkname(ik) = joinpath(unkdir, @sprintf("UNK%05d.%s", ik, ext))

    # read 1st kpoint to allocate matrices
    ik = 1
    _, Ψₖ = read_unk(unkname(ik))

    n_gx, n_gy, n_gz, _ = size(Ψₖ)
    size(Ψₖ, 4) != n_bands && error("incompatible n_bands")

    # WF in realspace
    W = zeros(eltype(A), n_gx * n_sx, n_gy * n_sy, n_gz * n_sz, n_wann)

    # generate X, Y, Z fractional coordinates relative to lattice
    X = zeros(T, n_gx * n_sx)
    Y = zeros(T, n_gy * n_sy)
    Z = zeros(T, n_gz * n_sz)
    X .= (n_gx * supercells[1][1]):(n_gx * (supercells[1][end] + 1) - 1)
    Y .= (n_gy * supercells[2][1]):(n_gy * (supercells[2][end] + 1) - 1)
    Z .= (n_gz * supercells[3][1]):(n_gz * (supercells[3][end] + 1) - 1)
    X ./= n_gx
    Y ./= n_gy
    Z ./= n_gz
    # fractional w.r.t. lattice (here unknown)
    Xg, Yg, Zg = ndgrid(X, Y, Z)

    """Modify W"""
    function add_k!(ik, Ψₖ)
        k = kpoints[:, ik]
        # rotate Ψ
        # * does not support high dimensional matrix multiplication,
        # I need to reshape it to 2D matrix
        ΨAₖ = reshape(reshape(Ψₖ, :, n_bands) * A[:, :, ik], n_gx, n_gy, n_gz, n_wann)
        for ix in axes(X, 1)
            for iy in axes(Y, 1)
                for iz in axes(Z, 1)
                    # r grid in fractional coordinates
                    r = [X[ix], Y[iy], Z[iz]]
                    # UNK file (and ΨAₖ) is the periodic part of Bloch wavefunction,
                    # need a factor exp(ikr)
                    f = exp(2π * im * k' * (r - R))
                    # make sure Ψ is normalized to 1 in unit cell
                    f /= (n_gx * n_gy * n_gz)^(1 / 2)
                    # ΨAₖ is only defined in the home unit cell, find corresponding indexes
                    # e.g. 1 -> 1, n_gx -> n_gx, n_gx + 1 -> 1
                    jx = (ix - 1) % n_gx + 1
                    jy = (iy - 1) % n_gy + 1
                    jz = (iz - 1) % n_gz + 1
                    W[ix, iy, iz, :] += f * ΨAₖ[jx, jy, jz, :]
                end
            end
        end
        return nothing
    end

    # 1st kpoint
    add_k!(ik, Ψₖ)

    for ik in 2:n_kpts
        unk = unkname(ik)
        _, Ψₖ = read_unk(unk)
        add_k!(ik, Ψₖ)
    end

    W ./= n_kpts

    return Xg, Yg, Zg, W
end

function read_realspace_wf(
    A::AbstractArray{Complex{T},3},
    kpoints::AbstractMatrix{T},
    n_supercells::Int=2,
    unkdir::String=".";
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    return read_realspace_wf(
        A, kpoints, [n_supercells, n_supercells, n_supercells], unkdir; R=R
    )
end

function read_realspace_wf(
    lattice::AbstractMatrix{T},
    A::AbstractArray{Complex{T},3},
    kpoints::AbstractMatrix{T},
    n_supercells::Int=2,
    unkdir::String=".";
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    X, Y, Z, W = read_realspace_wf(A, kpoints, n_supercells, unkdir; R=R)
    rgrid = RGrid(lattice, X, Y, Z)
    return rgrid, W
end

function read_realspace_wf(
    model::Model{T},
    n_supercells::Int=2,
    unkdir::String=".";
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    return read_realspace_wf(
        model.lattice, model.A, model.kpoints, n_supercells, unkdir; R=R
    )
end

"""
Write realspace Wannier functions to cube files.

seedname: the name prefix for cube files, e.g., `seedname_00001.cube`
A: gauge rotation matrix
atom_positions: fractional coordinates w.r.t. lattice
part: which part to plot? pass a Function, e.g. real, imag, abs2
"""
function write_realspace_wf(
    seedname::String,
    A::AbstractArray,
    kpoints::AbstractMatrix,
    lattice::AbstractMatrix,
    atom_positions::AbstractMatrix,
    atom_labels::AbstractVector{String};
    n_supercells::Int=2,
    unkdir::String=".",
    part::Function=real,
    format::Symbol=:xsf,
)
    format ∈ [:xsf, :cube] || error("format must be :xsf or :cube")

    rgrid, W = read_realspace_wf(lattice, A, kpoints, n_supercells, unkdir)
    n_wann = size(W, 4)

    for i in 1:n_wann
        f = fix_global_phase(W[:, :, :, i])
        W[:, :, :, i] .*= f
        r = compute_imre_ratio(W[:, :, :, i])
        @printf("Im/Re ratio of WF %4d = %10.4f\n", i, r)
    end

    # seems W90 always write the real part, so I use real as default
    W2 = part.(W)

    atom_numbers = get_atom_number(atom_labels)
    origin = origin(rgrid)
    spanvec = span_vectors(rgrid)

    if format == :xsf
        for i in 1:n_wann
            filename = @sprintf("%s_%05d.%s", seedname, i, "xsf")
            write_xsf(
                filename,
                lattice,
                atom_positions,
                atom_numbers,
                origin,
                spanvec,
                W2[:, :, :, i],
            )
        end
    elseif format == :cube
        positions_cart = lattice * atom_positions
        for i in 1:n_wann
            filename = @sprintf("%s_%05d.%s", seedname, i, "cube")
            write_cube(
                filename, positions_cart, atom_numbers, origin, spanvec, W2[:, :, :, i]
            )
        end
    end

    return nothing
end

function write_realspace_wf(
    seedname::String,
    model::Model;
    n_supercells::Int=2,
    unkdir::String=".",
    part::Function=real,
    format::Symbol=:xsf,
)
    return write_realspace_wf(
        seedname,
        model.A,
        model.kpoints,
        model.lattice,
        model.atom_positions,
        model.atom_labels;
        n_supercells,
        unkdir,
        part,
        format,
    )
end

@doc """Return X, Y, Z in cartesian coordinates"""
function cartesianize_xyz(rgrid::RGrid)
    XYZᶜ =
        rgrid.basis * vcat(reshape(rgrid.X, :)', reshape(rgrid.Y, :)', reshape(rgrid.Z, :)')
    Xᶜ = reshape(XYZᶜ[1, :], size(rgrid.X)...)
    Yᶜ = reshape(XYZᶜ[2, :], size(rgrid.X)...)
    Zᶜ = reshape(XYZᶜ[3, :], size(rgrid.X)...)
    return Xᶜ, Yᶜ, Zᶜ
end

@doc """
Compute WF moment (mean, variance, ...) in realspace.

Note WFs are defined in a supercell that is n_kpts times unit cell,
however, usuall we only calculate realspace WFs in a smaller supercell
that is 2^3 or 3^3 times unit cell (as defined by the `n_supercells` of
`read_realspace_wf`). Some times this is not sufficient if the WFs are
truncated by the boundries of the smaller supercell, thus the center
calculated by this function is inexact. In principle, we should calculate
centers in the n_kpts supercell, however, this is memory-consuming.

rgrid: realspace grid on which W is defined
W: Wannier functions
n: order of moment, e.g., 1 for WF center, 2 for variance, etc.
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

center(rgrid::RGrid, W::AbstractArray) = moment(rgrid, W, 1)
omega(rgrid::RGrid, W::AbstractArray) = moment(rgrid, W, 2) - center(rgrid, W) .^ 2

@doc """Position operator matrices computed with realspace WFs"""
function position(rgrid::RGrid, W::AbstractArray{T,4}) where {T<:Complex}
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

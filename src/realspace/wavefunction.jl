using Printf: @sprintf
using LazyGrids: ndgrid

export read_realspace_wf, write_realspace_wf

"""
    _get_unk_ext(unkdir::AbstractString)

Get the extension name of `UNK` files.

For example:
- `1` of `UNK00001.1`, for no-spin calcluation
- `NC` of `UNK00001.NC`, for non-collinear calculation
"""
function _get_unk_ext(unkdir::AbstractString)
    files = readdir(unkdir)
    i = findfirst(x -> startswith(x, "UNK"), files)
    isnothing(i) && error("No UNK files found in $unkdir")
    suffix = split(files[i], ".")[2]  # "1" or "NC"
    return suffix
end

"""
    read_realspace_wf(U, kpoints, n_supercells, unkdir="."; R=[0, 0, 0])

Read `UNK` files, rotate gauge, and generate real space WFs.

# Arguments
- `U`: `n_bands * n_wann * n_kpts`, gauge rotation matrix
- `kpoints`: `3 * n_kpts`, each column is a vector for fractional coordinates
- `n_supercells`: an integer or a vector of 3 integers for 3 directions along lattice vectors,
    defines the number of super cells where the real space WF lives
- `unkdir`: folder containing `UNK` files

# Keyword arguments
- `R`: fractional coordinates w.r.t lattice (actually integers), the cell for WF
    ``|n \\bm{R} \\rangle``, default is generating WF at home unit cell ``|n \\bm{0} \\rangle``

!!! tip

    See also the section [Normalization convention of WFs](@ref) for further explanation.
"""
function read_realspace_wf(
    U::AbstractArray{Complex{T},3},
    kpoints::AbstractMatrix{T},
    n_supercells::AbstractVector{Int},
    unkdir::AbstractString;
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    n_bands, n_wann, n_kpts = size(U)
    length(n_supercells) == 3 || error("n_supercells must be 3-vector")

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
    size(Ψₖ, 4) == n_bands || error("incompatible n_bands")

    # WF in realspace
    W = zeros(eltype(U), n_gx * n_sx, n_gy * n_sy, n_gz * n_sz, n_wann)

    # generate X, Y, Z fractional coordinates relative to lattice (here unknown)
    # actually X./n_gx, Y./n_gy, Z./n_gz are the fractional coordinates w.r.t lattice
    X = (supercells[1][1] * n_gx):((supercells[1][end] + 1) * n_gx - 1)
    Y = (supercells[2][1] * n_gy):((supercells[2][end] + 1) * n_gy - 1)
    Z = (supercells[3][1] * n_gz):((supercells[3][end] + 1) * n_gz - 1)
    # we start from origin, however W90 starts from the left of origin
    # uncomment these lines to exactly reproduce W90 output
    # seems subtracting 1 leads to much worse performance
    #  11.015472 seconds (136.09 M allocations: 5.570 GiB, 10.77% gc time)
    #   3.513259 seconds (25.60 M allocations: 2.728 GiB, 24.68% gc time)
    # X = X .- 1
    # Y = Y .- 1
    # Z = Z .- 1

    """Modify W"""
    function add_k!(ik, Ψₖ)
        k = kpoints[:, ik]
        # rotate Ψ
        # * does not support high dimensional matrix multiplication,
        # I need to reshape it to 2D matrix
        ΨUₖ = reshape(reshape(Ψₖ, :, n_bands) * U[:, :, ik], n_gx, n_gy, n_gz, n_wann)
        # make sure Ψ is normalized to 1 in unit cell
        # the QE output UNK needs this factor to be normalized
        # Note in QE, the r->G FFT has a factor 1/N,
        # while the G->r IFFT has factor 1,
        # but FFT/IFFT is unitary only if they have factor sqrt(1/N)
        f0 = 1 / sqrt(n_gx * n_gy * n_gz)
        for (iz, z) in enumerate(Z)
            # ΨUₖ is only defined in the home unit cell, find corresponding indexes
            # e.g. x=0 -> W[1,1,1], x=1 -> W[2,1,1], x=n_gx -> W[1,1,1], x=-1 -> W[end,1,1]
            jz = mod(z, n_gz) + 1
            for (iy, y) in enumerate(Y)
                jy = mod(y, n_gy) + 1
                for (ix, x) in enumerate(X)
                    jx = mod(x, n_gx) + 1
                    # r grid in fractional coordinates
                    r = [x / n_gx, y / n_gy, z / n_gz]
                    # UNK file (and ΨUₖ) is the periodic part of Bloch wavefunction,
                    # need a factor exp(ikr)
                    f = f0 * exp(2π * im * k' * (r - R))
                    W[ix, iy, iz, :] += f * ΨUₖ[jx, jy, jz, :]
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
    # to fractional coordinates
    Xg, Yg, Zg = ndgrid(X ./ n_gx, Y ./ n_gy, Z ./ n_gz)
    return Xg, Yg, Zg, W
end

function read_realspace_wf(
    U::AbstractArray{Complex{T},3},
    kpoints::AbstractMatrix{T},
    n_supercells::Int,
    unkdir::AbstractString;
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    return read_realspace_wf(
        U, kpoints, [n_supercells, n_supercells, n_supercells], unkdir; R=R
    )
end

"""
    read_realspace_wf(lattice, U, kpoints, n_supercells=2, unkdir="."; R=[0, 0, 0])

Read `UNK` files, rotate gauge, and generate real space WFs.

This is a more user-friendly version, which returns a tuple of `(RGrid, W)`,
where `RGrid` is the grid on which `W` is defined, and `W` is volumetric data for WFs.

# Arguments
- `lattice`: each column is a lattice vector
- `U`: `n_bands * n_wann * n_kpts`, gauge rotation matrix
"""
function read_realspace_wf(
    lattice::AbstractMatrix{T},
    U::AbstractArray{Complex{T},3},
    kpoints::AbstractMatrix{T},
    n_supercells::Union{AbstractArray{Int},Int}=2,
    unkdir::AbstractString=".";
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    X, Y, Z, W = read_realspace_wf(U, kpoints, n_supercells, unkdir; R=R)
    rgrid = RGrid(lattice, X, Y, Z)
    return rgrid, W
end

"""
    read_realspace_wf(model, U, n_supercells=2, unkdir="."; R=[0, 0, 0])

Read `UNK` files, rotate gauge, and generate real space WFs.

This is a most user-friendly version, use `lattice`, `U` and `kpoints` from `model`
and returns a tuple of `(RGrid, W)`,
where `RGrid` is the grid on which `W` is defined, and `W` is volumetric data for WFs.

# Arguments
- `model`: a `Model`
- `U`: `n_bands * n_wann * n_kpts`, gauge rotation matrix
"""
function read_realspace_wf(
    model::Model{T},
    U::AbstractArray{Complex{T},3},
    n_supercells::Union{AbstractArray{Int},Int}=2,
    unkdir::AbstractString=".";
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    return read_realspace_wf(model.lattice, U, model.kpoints, n_supercells, unkdir; R=R)
end

"""
    read_realspace_wf(model, n_supercells=2, unkdir="."; R=[0, 0, 0])

Read `UNK` files, rotate gauge, and generate real space WFs.

This is a most user-friendly version, use `lattice`, `U` and `kpoints` from `model`
and returns a tuple of `(RGrid, W)`,
where `RGrid` is the grid on which `W` is defined, and `W` is volumetric data for WFs.

# Arguments
- `model`: a `Model`
"""
function read_realspace_wf(
    model::Model{T},
    n_supercells::Union{AbstractArray{Int},Int}=2,
    unkdir::AbstractString=".";
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    return read_realspace_wf(model, model.U, n_supercells, unkdir; R=R)
end

"""
    write_realspace_wf(seedname, U, kpoints, lattice, atom_positions, atom_labels;
        n_supercells=2, unkdir=".", part=real, format=:xsf, wf_center=nothing)

Write real space WFs to `xsf` or `cube` files.

# Arguments
- `seedname`: the name prefix for `cube` files, e.g., `seedname_00001.cube`
- `U`: gauge rotation matrix
- `kpoints`: each column is a kpoint, fractional coordinates
- `lattice`: each column is a lattice vector
- `atom_positions`: each column is an atom position, fractional coordinates w.r.t. lattice
- `atom_labels`: each element is an atom label

# Keyword arguments
- `n_supercells`: number of supercells in each direction,
    equivalent to `wannier_plot_supercell` of `Wannier90`
- `unkdir`: directory of `UNK` files
- `part`: which part to plot? pass a `Function`, e.g. `real`, `imag`, `abs2`
- `format`: `:xsf` or `:cube`
- `wf_center`: WF centers in fractional coordinates w.r.t. lattice.
    Only used for `cube` format, add additional atoms around WF centers.

!!! note

    `Wannier90` only plot the real part of WFs, so `part=real` is the default.

!!! tip

    See also the section [Normalization convention of WFs](@ref) for further explanation.
"""
function write_realspace_wf(
    seedname::AbstractString,
    U::AbstractArray,
    kpoints::AbstractMatrix,
    lattice::AbstractMatrix,
    atom_positions::AbstractMatrix,
    atom_labels::AbstractVector{String};
    n_supercells::Union{AbstractArray{Int},Int}=2,
    unkdir::AbstractString=".",
    part::Function=real,
    format::Symbol=:xsf,
    wf_center::Union{Nothing,AbstractMatrix}=nothing,
)
    format ∈ [:xsf, :cube] || error("format must be :xsf or :cube")

    rgrid, W = read_realspace_wf(lattice, U, kpoints, n_supercells, unkdir)
    n_wann = size(W, 4)
    if format == :cube && wf_center === nothing
        # if not given, compute in realspace, lower accuracy. in fractional coordinates
        wf_center = inv(lattice) * center(rgrid, W)
    end

    # In principle, since we normalize Bloch wavefunction inside one unit cell,
    # thus the WF is normalized inside the n_kpts times unit cell.
    # However, W90 does not check the normalization of UNK files, it just plainly
    # mulitply the unitary matrices. But in read_readspace_wf, we normalized
    # the unk files, so here we need to multiply a factor to reproduce the
    # W90 output XSF or cube files.
    if isa(n_supercells, Integer)
        n_supcells = n_supercells^3
    else
        n_supcells = prod(n_supercells)
    end
    W .*= sqrt(length(W) / n_wann / n_supcells)

    for i in 1:n_wann
        f = fix_global_phase(W[:, :, :, i])
        W[:, :, :, i] .*= f
        r = compute_imre_ratio(W[:, :, :, i])
        @printf("Max Im/Re ratio of WF %4d = %10.4f\n", i, r)
    end

    # seems W90 always write the real part, so I use real as default
    W2 = part.(W)

    atom_numbers = get_atom_number(atom_labels)

    if format == :xsf
        for i in 1:n_wann
            filename = @sprintf("%s_%05d.%s", seedname, i, "xsf")
            write_xsf(
                filename, lattice, atom_positions, atom_numbers, rgrid, W2[:, :, :, i]
            )
        end
    elseif format == :cube
        for i in 1:n_wann
            filename = @sprintf("%s_%05d.%s", seedname, i, "cube")
            write_cube(
                filename,
                lattice,
                atom_positions,
                atom_numbers,
                wf_center[:, i],
                rgrid,
                W2[:, :, :, i],
            )
        end
    end

    return nothing
end

"""
    write_realspace_wf(seedname, model; n_supercells=2, unkdir=".", part=real, format=:xsf)

Write real space WFs to `xsf` or `cube` files.

This is a user-friendly version that use `model` to fill the arguments of `write_realspace_wf`.
"""
function write_realspace_wf(
    seedname::String,
    model::Model;
    n_supercells::Union{AbstractArray{Int},Int}=2,
    unkdir::AbstractString=".",
    part::Function=real,
    format::Symbol=:xsf,
)
    wf_center = nothing
    if format == :cube
        # compute in recip space, more accurate than realspace
        wf_center = inv(model.lattice) * center(model)
    end
    return write_realspace_wf(
        seedname,
        model.U,
        model.kpoints,
        model.lattice,
        model.atom_positions,
        model.atom_labels;
        n_supercells,
        unkdir,
        part,
        format,
        wf_center,
    )
end

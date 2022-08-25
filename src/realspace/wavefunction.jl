using Printf: @sprintf
using LazyGrids: ndgrid

export read_realspace_wf, write_realspace_wf

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
    W = zeros(eltype(A), n_gx * n_sx, n_gy * n_sy, n_gz * n_sz, n_wann)

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
        ΨAₖ = reshape(reshape(Ψₖ, :, n_bands) * A[:, :, ik], n_gx, n_gy, n_gz, n_wann)
        # make sure Ψ is normalized to 1 in unit cell
        # the QE output UNK needs this factor to be normalized
        # Note in QE, the r->G FFT has a factor 1/N,
        # while the G->r IFFT has factor 1,
        # but FFT/IFFT is unitary only if they have factor sqrt(1/N)
        f0 = 1 / sqrt(n_gx * n_gy * n_gz)
        for (iz, z) in enumerate(Z)
            # ΨAₖ is only defined in the home unit cell, find corresponding indexes
            # e.g. x=0 -> W[1,1,1], x=1 -> W[2,1,1], x=n_gx -> W[1,1,1], x=-1 -> W[end,1,1]
            jz = mod(z, n_gz) + 1
            for (iy, y) in enumerate(Y)
                jy = mod(y, n_gy) + 1
                for (ix, x) in enumerate(X)
                    jx = mod(x, n_gx) + 1
                    # r grid in fractional coordinates
                    r = [x / n_gx, y / n_gy, z / n_gz]
                    # UNK file (and ΨAₖ) is the periodic part of Bloch wavefunction,
                    # need a factor exp(ikr)
                    f = f0 * exp(2π * im * k' * (r - R))
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
    # to fractional coordinates
    Xg, Yg, Zg = ndgrid(X ./ n_gx, Y ./ n_gy, Z ./ n_gz)
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

    # In principle, since we normalize Bloch wavefunction inside one unit cell,
    # thus the WF is normalized inside the n_kpts times unit cell.
    # However, W90 does not check the normalization of UNK files, it just plainly
    # mulitply the unitary matrices. But in read_readspace_wf, we normalized
    # the unk files, so here we need to multiply a factor to reproduce the
    # W90 output XSF or cube files.
    W .*= sqrt(length(W) / n_wann / n_supercells^3)

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
        positions_cart = lattice * atom_positions
        for i in 1:n_wann
            filename = @sprintf("%s_%05d.%s", seedname, i, "cube")
            write_cube(filename, positions_cart, atom_numbers, rgrid, W2[:, :, :, i])
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

using Printf: @sprintf

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
            push!(supercells, (-d + 1):d)
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

    # generate X, Y, Z fractional coordinates
    X = zeros(T, n_gx * n_sx)
    Y = zeros(T, n_gy * n_sy)
    Z = zeros(T, n_gz * n_sz)
    X .= (n_gx * supercells[1][1]):(n_gx * (supercells[1][end] + 1) - 1)
    Y .= (n_gy * supercells[2][1]):(n_gy * (supercells[2][end] + 1) - 1)
    Z .= (n_gz * supercells[3][1]):(n_gz * (supercells[3][end] + 1) - 1)
    X ./= n_gx
    Y ./= n_gy
    Z ./= n_gz

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
                    # ΨAₖ is only defined in the home unit cell, find corresponding indexes
                    # e.g. 1 -> 1, n_gx -> n_gx, n_gx + 1 -> 1
                    iix = (ix - 1) % n_gx + 1
                    iiy = (iy - 1) % n_gy + 1
                    iiz = (iz - 1) % n_gz + 1
                    W[ix, iy, iz, :] += f * ΨAₖ[iix, iiy, iiz, :]
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

    # X, Y, Z are in fractional coordinates
    return X, Y, Z, W
end

function read_realspace_wf(
    A::AbstractArray{Complex{T},3},
    kpoints::AbstractMatrix{T},
    n_supercells::Int=2,
    unkdir::String=".",
) where {T<:Real}
    return read_realspace_wf(A, kpoints, [n_supercells, n_supercells, n_supercells], unkdir)
end

"""
Write realspace Wannier functions to cube files.

seedname: the name prefix for cube files, e.g., `seedname_00001.cube`
A: gauge rotation matrix
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

    X, Y, Z, W = read_realspace_wf(A, kpoints, n_supercells, unkdir)
    n_wann = size(W, 4)

    for i in 1:n_wann
        fix_global_phase!(W[:, :, :, i])
        r = compute_imre_ratio(W[:, :, :, i])
        @printf("Im/Re ratio of WF %4d = %10.4f\n", i, r)
    end

    # seems W90 always write the real part, so I use real as default
    W2 = part.(W)

    atom_numbers = get_atom_number(atom_labels)

    if format == :xsf
        ext = "xsf"
        func = write_xsf
    elseif format == :cube
        ext = "cube"
        func = write_cube
    end
    for i in 1:n_wann
        filename = @sprintf("%s_%05d.%s", seedname, i, ext)
        func(filename, atom_positions, atom_numbers, lattice, X, Y, Z, W2[:, :, :, i])
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

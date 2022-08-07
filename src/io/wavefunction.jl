using Printf: @sprintf

"""
Read unk files, and generate realspace Wannier functions.

A: n_bands x n_wann x n_kpts, rotation matrix
kpoints: 3 x n_kpts, fractional coordinates
"""
function read_realspace_wf(
    A::AbstractArray{Complex{T},3},
    kpoints::AbstractMatrix{T},
    n_supercells::AbstractVector{Int},
    unkdir::String=".",
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

    unkname(ik) = joinpath(unkdir, "UNK" * @sprintf("%05d", ik) * ".1")

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

    # preallocate a reshaped matrix
    n_g = n_gx * n_gy * n_gz
    ΨAₖ = zeros(eltype(A), n_g, n_wann)

    """Modify W, X, Y, Z"""
    function add_k!(ik, Ψₖ)
        k = kpoints[:, ik]
        # unk is the periodic part of Bloch wavefunction,
        # need a factor exp(ikr)
        for z in 1:n_gz
            for y in 1:n_gy
                for x in 1:n_gx
                    # r grid in fractional coordinates
                    r = [(x - 1) / n_gx, (y - 1) / n_gy, (z - 1) / n_gz]
                    Ψₖ[x, y, z, :] *= exp(2π * im * k' * r)
                end
            end
        end
        # rotate Ψ
        # * does not support high dimensional matrix multiplication,
        # I need to reshape it to 2D matrix
        ΨAₖ .= reshape(Ψₖ, n_g, n_bands) * A[:, :, ik]

        for sx in supercells[1]
            for sy in supercells[2]
                for sz in supercells[3]
                    # fractional coordinates
                    R = [sx, sy, sz]
                    f = exp(-2π * im * k' * R)
                    # starting index in W
                    ix = n_gx * (sx - supercells[1][1]) + 1
                    iy = n_gy * (sy - supercells[2][1]) + 1
                    iz = n_gz * (sz - supercells[3][1]) + 1
                    # ending index in W
                    jx = ix + n_gx - 1
                    jy = iy + n_gy - 1
                    jz = iz + n_gz - 1
                    #
                    X[ix:jx] = (n_gx * sx):(n_gx * (sx + 1) - 1)
                    Y[iy:jy] = (n_gy * sy):(n_gy * (sy + 1) - 1)
                    Z[iz:jz] = (n_gz * sz):(n_gz * (sz + 1) - 1)
                    #
                    W[ix:jx, iy:jy, iz:jz, :] += f * reshape(ΨAₖ, n_gx, n_gy, n_gz, n_wann)
                end
            end
        end
        X ./= n_gx
        Y ./= n_gy
        Z ./= n_gz
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
function write_realspace_wf_cube(
    seedname::String,
    A::AbstractArray,
    kpoints::AbstractMatrix,
    lattice::AbstractMatrix,
    atom_positions::AbstractMatrix,
    atom_labels::AbstractVector{String};
    n_supercells::Int=2,
    unkdir::String=".",
    part::Function=real,
)
    X, Y, Z, W = read_realspace_wf(A, kpoints, n_supercells, unkdir)
    n_wann = size(W, 4)

    for i in 1:n_wann
        fix_global_phase!(W[:, :, :, i])
        r = compute_imre_ratio(W[:, :, :, i])
        @printf("Im/Re ratio of WF %4d = %10.4f\n", i, r)
    end

    # seems W90 always write the real part, so I use real as default
    W2 = part(W)

    atom_numbers = get_atom_number(atom_labels)

    for i in 1:n_wann
        filename = seedname * "_" * @sprintf("%05d", i) * ".cube"
        write_cube(filename, atom_positions, atom_numbers, lattice, X, Y, Z, W2[:, :, :, i])
    end

    return nothing
end

function write_realspace_wf_cube(
    seedname::String,
    model::Model;
    n_supercells::Int=2,
    unkdir::String=".",
    part::Function=real,
)
    return write_realspace_wf_cube(
        seedname,
        model.A,
        model.kpoints,
        model.lattice,
        model.atom_positions,
        model.atom_labels;
        n_supercells,
        unkdir,
        part,
    )
end

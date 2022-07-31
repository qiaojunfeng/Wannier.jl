using Printf: @printf

"""Struct for storing matrices in seedname.chk file."""
struct Chk{T<:Real}
    header::String

    # These commented variables are written inside chk file,
    # but I move them to the end of the struct and use a constructor
    # to set them according to the shape of the corresponding arrays.

    # n_bands:: Int
    # n_exclude_bands:: Int

    # 1D array, size: num_exclude_bands
    exclude_bands::Vector{Int}

    # 2D array, size: 3 x 3, each column is a lattice vector
    lattice::Mat3{T}

    # 2D array, size: 3 x 3, each column is a lattice vector
    recip_lattice::Mat3{T}

    # n_kpts: Int

    # List of 3 int
    kgrid::Vec3{Int}

    # 2D array, size: 3 x n_kpts, each column is a vector
    kpoints::Matrix{T}

    # n_bvecs:: Int

    # n_wann: int

    checkpoint::String

    have_disentangled::Bool

    # omega invariant
    ΩI::T

    # Bands taking part in disentanglement, not frozen bands!
    # This is needed since W90 puts all the disentanglement bands
    # in the first several rows of Uᵈ,
    # (and the first few columns of Uᵈ are the frozen bands)
    # so directly multiplying eigenvalues e.g.
    # (Uᵈ * U)' * diag(eigenvalues) * (Uᵈ * U) is wrong!
    # 2D bool array, size: n_bands x n_kpts
    dis_bands::BitMatrix

    # 1D int array, size: n_kpts
    # n_dimfrozen:: Vector{Int}

    # u_matrix_opt, 3D array, size: num_bands x num_wann x num_kpts
    Uᵈ::Array{Complex{T},3}

    # u_matrix, 3D array, size: num_wann x num_wann x num_kpts
    U::Array{Complex{T},3}

    # m_matrix, 4D array, size:num_wann x num_wann x n_bvecs x num_kpts
    M::Array{Complex{T},4}

    # wannier_centres, 2D array, size: 3 x num_wann
    r::Matrix{T}

    # wannier_spreads, 1D array, size: num_wann
    ω::Vector{T}

    # these variables are auto-set in constructor
    n_bands::Int
    n_exclude_bands::Int
    n_kpts::Int
    n_bvecs::Int
    n_wann::Int
    n_dis::Vector{Int}
end

function Chk(
    header::String,
    exclude_bands::Vector{Int},
    lattice::Mat3{T},
    recip_lattice::Mat3{T},
    kgrid::Vec3{Int},
    kpoints::Matrix{T},
    checkpoint::String,
    have_disentangled::Bool,
    ΩI::T,
    dis_bands::BitMatrix,
    Uᵈ::Array{Complex{T},3},
    U::Array{Complex{T},3},
    M::Array{Complex{T},4},
    r::Matrix{T},
    ω::Vector{T},
) where {T<:Real}
    if have_disentangled
        n_bands = size(Uᵈ, 1)
    else
        n_bands = size(U, 1)
    end

    n_exclude_bands = length(exclude_bands)

    n_kpts = size(M, 4)

    n_bvecs = size(M, 3)

    n_wann = size(U, 1)

    if have_disentangled
        n_dis = zeros(Int, n_kpts)
        for ik in 1:n_kpts
            n_dis[ik] = count(dis_bands[:, ik])
        end
    else
        n_dis = zeros(Int, 0)
    end

    return Chk(
        header,
        exclude_bands,
        lattice,
        recip_lattice,
        kgrid,
        kpoints,
        checkpoint,
        have_disentangled,
        ΩI,
        dis_bands,
        Uᵈ,
        U,
        M,
        r,
        ω,
        n_bands,
        n_exclude_bands,
        n_kpts,
        n_bvecs,
        n_wann,
        n_dis,
    )
end

@doc raw"""
Read formatted (not Fortran binary) CHK file.
"""
function read_chk(filename::String)
    if isbinary_file(filename)
        error("$filename is a binary file? Consider using `w90chk2chk.x`?")
    end

    @info "Reading chk file:" filename
    println()

    io = open(filename)

    # strip and read line
    srline() = strip(readline(io))

    # Read formatted chk file
    header = String(srline())

    n_bands = parse(Int, srline())

    n_exclude_bands = parse(Int, srline())

    exclude_bands = zeros(Int, n_exclude_bands)

    if n_exclude_bands > 0
        for i in 1:n_exclude_bands
            exclude_bands[i] = parse(Int, srline())
        end
    end

    # Each column is a lattice vector
    # but W90 writes x components first, then y, z. Not a1 first, then a2, a3.
    line = parse.(Float64, split(srline()))
    lattice = Mat3{Float64}(reshape(line, (3, 3))')

    # Each column is a lattice vector
    line = parse.(Float64, split(srline()))
    recip_lattice = Mat3{Float64}(reshape(line, (3, 3))')

    n_kpts = parse(Int, srline())

    kgrid = Vec3{Int}(parse.(Int, split(srline())))

    kpoints = zeros(Float64, 3, n_kpts)
    for ik in 1:n_kpts
        kpoints[:, ik] = parse.(Float64, split(srline()))
    end

    n_bvecs = parse(Int, srline())

    n_wann = parse(Int, srline())

    checkpoint = String(srline())

    # 1 -> True, 0 -> False
    have_disentangled = Bool(parse(Int, srline()))

    if have_disentangled
        # omega_invariant
        ΩI = parse(Float64, srline())

        dis_bands = falses(n_bands, n_kpts)
        for ik in 1:n_kpts
            for ib in 1:n_bands
                # 1 -> True, 0 -> False
                dis_bands[ib, ik] = Bool(parse(Int, srline()))
            end
        end

        n_dis = zeros(Int, n_kpts)
        for ik in 1:n_kpts
            n_dis[ik] = parse(Int, srline())
            @assert n_dis[ik] == count(dis_bands[:, ik])
        end

        # u_matrix_opt
        Uᵈ = zeros(ComplexF64, n_bands, n_wann, n_kpts)
        for ik in 1:n_kpts
            for iw in 1:n_wann
                for ib in 1:n_bands
                    vals = parse.(Float64, split(srline()))
                    Uᵈ[ib, iw, ik] = vals[1] + im * vals[2]
                end
            end
        end

    else
        ΩI = -1.0
        dis_bands = falses(0, 0)
        n_dis = zeros(Int, 0)
        Uᵈ = zeros(ComplexF64, 0, 0, 0)
    end

    # u_matrix
    U = zeros(ComplexF64, n_wann, n_wann, n_kpts)
    for ik in 1:n_kpts
        for iw in 1:n_wann
            for ib in 1:n_wann
                vals = parse.(Float64, split(srline()))
                U[ib, iw, ik] = vals[1] + im * vals[2]
            end
        end
    end

    #  m_matrix
    M = zeros(ComplexF64, n_wann, n_wann, n_bvecs, n_kpts)
    for ik in 1:n_kpts
        for inn in 1:n_bvecs
            for iw in 1:n_wann
                for ib in 1:n_wann
                    vals = parse.(Float64, split(srline()))
                    M[ib, iw, inn, ik] = vals[1] + im * vals[2]
                end
            end
        end
    end

    # wannier_centres
    r = zeros(Float64, 3, n_wann)
    for iw in 1:n_wann
        r[:, iw] = parse.(Float64, split(srline()))
    end

    # wannier_spreads
    ω = zeros(Float64, n_wann)
    for iw in 1:n_wann
        ω[iw] = parse(Float64, srline())
    end

    close(io)

    return Chk(
        header,
        exclude_bands,
        lattice,
        recip_lattice,
        kgrid,
        kpoints,
        checkpoint,
        have_disentangled,
        ΩI,
        dis_bands,
        Uᵈ,
        U,
        M,
        r,
        ω,
    )
end

@doc raw"""
Write formatted (not Fortran binary) CHK file.
"""
function write_chk(filename::String, chk::Chk)
    @info "Writing chk file:" filename
    io = open(filename, "w")

    n_bands = chk.n_bands
    n_wann = chk.n_wann
    n_kpts = chk.n_kpts
    n_bvecs = chk.n_bvecs

    # Read formatted chk file
    @printf(io, "%33s\n", chk.header)

    @printf(io, "%d\n", n_bands)

    @printf(io, "%d\n", chk.n_exclude_bands)

    if chk.n_exclude_bands > 0
        for i in 1:(chk.n_exclude_bands)
            @printf(io, "%d\n", chk.exclude_bands[i])
        end
    end

    # Each column is a lattice vector
    # but W90 writes x components first, then y, z. Not a1 first, then a2, a3.
    for v in reshape(chk.lattice', 9)
        @printf(io, "%25.17f", v)
    end
    @printf(io, "\n")

    # Each column is a lattice vector
    for v in reshape(chk.recip_lattice', 9)
        @printf(io, "%25.17f", v)
    end
    @printf(io, "\n")

    @printf(io, "%d\n", n_kpts)

    @printf(io, "%d %d %d\n", chk.kgrid...)

    for ik in 1:n_kpts
        @printf(io, "%25.17f %25.17f %25.17f\n", chk.kpoints[:, ik]...)
    end

    @printf(io, "%d\n", n_bvecs)

    @printf(io, "%d\n", n_wann)

    # left-justified
    @printf(io, "%-20s\n", chk.checkpoint)

    # 1 -> True, 0 -> False
    # v = chk.have_disentangled ? 1 : 0
    @printf(io, "%d\n", chk.have_disentangled)

    if chk.have_disentangled
        # omega_invariant
        @printf(io, "%25.17f\n", chk.ΩI)

        for ik in 1:n_kpts
            for ib in 1:n_bands
                # 1 -> True, 0 -> False
                @printf(io, "%d\n", chk.dis_bands[ib, ik])
            end
        end

        for ik in 1:n_kpts
            @printf(io, "%d\n", chk.n_dis[ik])
        end

        # u_matrix_opt
        for ik in 1:n_kpts
            for iw in 1:n_wann
                for ib in 1:n_bands
                    v = chk.Uᵈ[ib, iw, ik]
                    @printf(io, "%25.17f %25.17f\n", real(v), imag(v))
                end
            end
        end
    end

    # u_matrix
    for ik in 1:n_kpts
        for iw in 1:n_wann
            for ib in 1:n_wann
                v = chk.U[ib, iw, ik]
                @printf(io, "%25.17f %25.17f\n", real(v), imag(v))
            end
        end
    end

    #  m_matrix
    for ik in 1:n_kpts
        for inn in 1:n_bvecs
            for iw in 1:n_wann
                for ib in 1:n_wann
                    v = chk.M[ib, iw, inn, ik]
                    @printf(io, "%25.17f %25.17f\n", real(v), imag(v))
                end
            end
        end
    end

    # wannier_centres
    for iw in 1:n_wann
        @printf(io, "%25.17f %25.17f %25.17f\n", chk.r[:, iw]...)
    end

    # wannier_spreads
    for iw in 1:n_wann
        @printf(io, "%25.17f\n", chk.ω[iw])
    end

    close(io)

    @info "Written to file: $(filename)"
    println()
    return nothing
end

@doc raw"""
Write formatted (not Fortran binary) CHK file.
"""
function write_chk(
    filename::String,
    model::Model;
    exclude_bands::Union{AbstractVector{Int},Nothing}=nothing,
)
    header = @sprintf "Created by Wannier.jl %s" string(now())

    if isnothing(exclude_bands)
        exclude_bands = Vector{Int}()
    end

    checkpoint = "postwann"
    have_disentangled = true
    Ω = omega(model)
    dis_bands = trues(model.n_bands, model.n_kpts)
    Uᵈ = model.A
    U = eyes_A(eltype(Uᵈ), model.n_wann, model.n_kpts)
    M = rotate_M(model.M, model.bvectors.kpb_k, model.A)

    chk = Chk(
        header,
        exclude_bands,
        model.lattice,
        model.recip_lattice,
        model.kgrid,
        model.kpoints,
        checkpoint,
        have_disentangled,
        Ω.ΩI,
        dis_bands,
        Uᵈ,
        U,
        M,
        Ω.r,
        Ω.ω,
    )

    write_chk(filename, chk)
    return nothing
end

"""
Construct a model from CHK file
"""
function get_model(chk::Chk)
    atom_positions = zeros(Float64, 3, 0)
    atom_labels = Vector{String}()

    recip_lattice = get_recip_lattice(chk.lattice)
    # I try to generate bvectors, but it might happen that the generated bvectors
    # are different from the calculation corresponding to the chk file,
    # e.g. kmesh_tol is different
    bvectors = get_bvectors(chk.kpoints, recip_lattice)
    if bvectors.n_bvecs != chk.n_bvecs
        error("Number of bvectors is different from the number in the chk file")
    end

    frozen_bands = falses(chk.n_bands, chk.n_kpts)

    # the M in chk is already rotated by the A matrix
    M = chk.M
    # so I set A matrix as identity
    A = eyes_A(eltype(M), chk.n_wann, chk.n_kpts)

    # no eig in chk file
    E = zeros(Float64, chk.n_wann, chk.n_kpts)

    model = Model(
        chk.lattice,
        atom_positions,
        atom_labels,
        chk.kgrid,
        chk.kpoints,
        bvectors,
        frozen_bands,
        M,
        A,
        E,
    )
    return model
end

"""
Extract AMN matrices from chk.
"""
function get_A(chk::Chk)
    n_kpts = chk.n_kpts
    n_bands = chk.n_bands
    n_wann = chk.n_wann

    U = similar(chk.U, n_bands, n_wann, n_kpts)

    if !chk.have_disentangled
        U .= chk.U
        return U
    end

    # need to permute wavefunctions since Uᵈ is stored in a way that
    # the bands taking part in disentanglement are in the first few rows
    Iᵏ = Matrix{eltype(U)}(I, n_bands, n_bands)

    for ik in 1:n_kpts
        # sortperm is stable, and
        # need descending order (dis bands at the front)
        p = sortperm(chk.dis_bands[:, ik]; order=Base.Order.Reverse)
        # usually we don't need this permutation, but if
        # 1. the dis_win_min > minimum(E), then these below
        #    dis_win_min bands are shifted to the last rows of Uᵈ
        # 2. use projectability disentanglement, then
        #    there might be cases that the lower (bonding) and
        #    higher (anti-bonding) bands participate in disentanglement,
        #    but some low-projectability bands are excluded from
        #    disentanglement, then these low-proj bands are shifted to
        #    the last rows of Uᵈ
        # so we need to permute the Bloch states before multiplying Uᵈ
        # Uᵈ: semi-unitary matrices from disentanglement
        # U: unitary matrices from maximal localization
        U[:, :, ik] = Iᵏ[:, p] * chk.Uᵈ[:, :, ik] * chk.U[:, :, ik]
    end

    return U
end

using Printf: @printf
using Dates: now
import DelimitedFiles as Dlm
import LinearAlgebra as LA


function read_win(filename::String)::Dict
    @info "Reading win file: $filename"
    io = open(filename)

    num_wann = missing
    num_bands = missing
    mp_grid = missing
    unit_cell = missing
    kpoints = missing
    kpoint_path = missing

    parse_array(line::String) = map(x -> parse(Float64, x), split(line))
    read_array(f::IOStream) = parse_array(readline(f))

    while !eof(io)
        line = readline(io)
        # handle case insensitive win files (relic of Fortran)
        line = strip(lowercase(line))

        line = replace(line, "=" => " ", ":" => " ", "," => " ")

        i = findfirst(r"!|#", line)
        if i != nothing
            line = strip(line[1:i.start-1])
        end

        if isempty(line)
            continue
        end

        if occursin("mp_grid", line)
            mp_grid = map(x -> parse(Int, x), split(line)[2:4])
            n_kpts = prod(mp_grid)
        elseif occursin("num_bands", line)
            num_bands = parse(Int, split(line)[2])
        elseif occursin("num_wann", line)
            num_wann = parse(Int, split(line)[2])
        elseif occursin(r"begin\s+unit_cell_cart", line)
            unit_cell = zeros(Float64, 3, 3)
            line = readline(io)
            unit = strip(lowercase(line))
            if !startswith(unit, "b") && !startswith(unit, "a")
                unit = "ang"
            else
                line = readline(io)
            end
            for i = 1:3
                # in win file, each line is a lattice vector, here it is stored as column vec
                unit_cell[:, i] = parse_array(line)
                line = readline(io)
            end
            if startswith(unit, "b")
                # convert to angstrom
                unit_cell .*= Bohr
            end
            unit_cell = Mat3{Float64}(unit_cell)
        elseif occursin(r"begin\s+kpoints", line)
            line = strip(readline(io))
            # kpoints block might be defined before mp_grid!
            # I need to read all lines and get n_kpts
            lines = Vector{String}()
            while !occursin(r"end\s+kpoints", line)
                push!(lines, line)
                line = strip(readline(io))
            end

            n_kpts = length(lines)
            kpoints = Matrix{Float64}(undef, 3, n_kpts)
            for i = 1:n_kpts
                # There might be weight at 4th column, but we don't use it.
                kpoints[:, i] = parse_array(lines[i])[1:3]
            end
        elseif occursin(r"begin\skpoint_path", line)
            StrVec3 = Pair{String,Vec3{Float64}}
            kpoint_path = Vector{Tuple{StrVec3,StrVec3}}()

            line = strip(readline(io))
            while !occursin(r"end\s+kpoint_path", line)
                l = split(line)
                length(l) != 8 && error("Invalid kpoint_path line: $line")

                start_label = String(l[1])
                start_kpt = Vec3{Float64}(parse.(Float64, l[2:4]))
                end_label = String(l[5])
                end_kpt = Vec3{Float64}(parse.(Float64, l[6:8]))
                push!(kpoint_path, (start_label => start_kpt, end_label => end_kpt))

                line = strip(readline(io))
            end
        end
    end
    close(io)

    ismissing(num_wann) && error("num_wann not found in win file")
    num_wann <= 0 && error("num_wann must be positive")

    ismissing(num_bands) && (num_bands = num_wann)
    num_bands <= 0 && error("num_bands must be positive")

    any(i -> i <= 0, mp_grid) && error("mp_grid must be positive")

    any(x -> ismissing(x), unit_cell) && error("unit_cell not found in win file")

    @info "$filename OK" num_wann num_bands mp_grid[1] mp_grid[2] mp_grid[3]

    Dict(
        "num_wann" => num_wann,
        "num_bands" => num_bands,
        "mp_grid" => mp_grid,
        "kpoints" => kpoints,
        "unit_cell" => unit_cell,
        "kpoint_path" => kpoint_path,
    )
end


function read_mmn(filename::String)
    @info "Reading mmn file: $filename"
    io = open(filename)

    # skip header
    header = readline(io)

    line = readline(io)
    n_bands, n_kpts, n_bvecs = map(x -> parse(Int64, x), split(line))

    # overlap matrix
    M = zeros(ComplexF64, n_bands, n_bands, n_bvecs, n_kpts)
    # for each point, list of neighbors, (K) representation
    kpb_k = zeros(Int64, n_bvecs, n_kpts)
    kpb_b = zeros(Int64, 3, n_bvecs, n_kpts)

    while !eof(io)
        for ib = 1:n_bvecs
            line = readline(io)
            arr = split(line)
            k = parse(Int64, arr[1])
            kpb_k[ib, k] = parse(Int64, arr[2])
            kpb_b[:, ib, k] = map(x -> parse(Int64, x), arr[3:5])
            for n = 1:n_bands
                for m = 1:n_bands
                    line = readline(io)
                    arr = split(line)
                    o = parse(Float64, arr[1]) + im * parse(Float64, arr[2])
                    M[m, n, ib, k] = o
                end
            end
        end
    end

    @assert !any(isnan.(M))

    close(io)

    @info "$filename OK" header n_bands n_bvecs n_kpts

    M, kpb_k, kpb_b
end


"""
Output amn file
"""
function write_mmn(
    filename::String,
    M::Array{ComplexF64,4},
    kpb_k::Matrix{Int},
    kpb_b::Array{Int,3},
)

    n_bands, _, n_bvecs, n_kpts = size(M)
    n_bands != size(M, 2) && error("M must be n_bands x n_bands x n_bvecs x n_kpts")

    size(kpb_k) != (n_bvecs, n_kpts) && error("kpb_k has wrong size")
    size(kpb_b) != (3, n_bvecs, n_kpts) && error("kpb_b has wrong size")

    open(filename, "w") do io
        write(io, "Created by Wannier.jl ", string(now()), "\n")

        @printf(io, "    %d   %d    %d \n", n_bands, n_kpts, n_bvecs)

        for ik = 1:n_kpts
            for ib = 1:n_bvecs
                @printf(io, "%d %d %d %d %d\n", ik, kpb_k[ib, ik], kpb_b[:, ib, ik]...)

                for n = 1:n_bands
                    for m = 1:n_bands
                        o = M[m, n, ib, ik]
                        @printf(io, "  %16.12f  %16.12f\n", real(o), imag(o))
                    end
                end
            end
        end

    end

    @info "Written to file: $(filename)"

    nothing
end


function read_amn(filename::String)
    @info "Reading $filename"

    io = open("$filename")

    header = readline(io)

    arr = split(readline(io))
    n_bands, n_kpts, n_wann = parse.(Int64, arr[1:3])

    A = zeros(ComplexF64, n_bands, n_wann, n_kpts)

    while !eof(io)
        line = readline(io)
        arr = split(line)

        m, n, k = parse.(Int64, arr[1:3])

        a = parse(Float64, arr[4]) + im * parse(Float64, arr[5])
        A[m, n, k] = a
    end

    close(io)

    @info "$filename OK" header n_bands n_wann n_kpts

    A
end


"""
Output amn file
"""
function write_amn(filename::String, A::Array{ComplexF64,3})
    n_bands, n_wann, n_kpts = size(A)

    io = open(filename, "w")

    write(io, "Created by Wannier.jl ", string(now()), "\n")

    @printf(io, "%3d %4d %4d\n", n_bands, n_kpts, n_wann)

    for ik = 1:n_kpts
        for iw = 1:n_wann
            for ib = 1:n_bands
                a = A[ib, iw, ik]
                @printf(io, "%5d %4d %4d  %16.12f  %16.12f\n", ib, iw, ik, real(a), imag(a))
            end
        end
    end

    close(io)

    @info "Written to file: $(filename)"

    nothing
end


function read_eig(filename::String)
    @info "Reading $filename"

    lines = open(filename) do io
        readlines(io)
    end

    n_lines = length(lines)
    idx_b = zeros(Int, n_lines)
    idx_k = zeros(Int, n_lines)
    eig = zeros(Float64, n_lines)

    for i = 1:n_lines
        arr = split(lines[i])
        idx_b[i] = parse(Int, arr[1])
        idx_k[i] = parse(Int, arr[2])
        eig[i] = parse(Float64, arr[3])
    end

    # find unique elements
    n_bands = length(Set(idx_b))
    n_kpts = length(Set(idx_k))
    E = reshape(eig, (n_bands, n_kpts))

    # check eigenvalues should be in order
    # some times there are small noise
    round_digits(x) = round(x, digits = 9)
    for ik = 1:n_kpts
        @assert issorted(E[:, ik], by = round_digits) display(ik) display(E[:, ik])
    end

    @info "$filename OK" n_bands n_kpts

    E
end


"""
Write eig file
"""
function write_eig(filename::String, E::Matrix{T}) where {T<:Real}
    n_bands, n_kpts = size(E)

    open(filename, "w") do io

        for ik = 1:n_kpts
            for ib = 1:n_bands
                @printf(io, "%5d%5d%18.12f\n", ib, ik, E[ib, ik])
            end
        end

    end

    @info "Written to file: $(filename)"

    nothing
end


@doc raw"""
read win, and optionally amn, mmn, eig
"""
function read_seedname(
    seedname::String;
    amn::Bool = true,
    mmn::Bool = true,
    eig::Bool = true,
)
    win = read_win("$seedname.win")

    n_bands = win["num_bands"]
    n_wann = win["num_wann"]
    kpoints = win["kpoints"]
    kgrid = Vec3{Int}(win["mp_grid"])
    n_kpts = prod(kgrid)

    lattice = win["unit_cell"]
    recip_lattice = get_recip_lattice(lattice)

    bvectors = get_bvectors(kpoints, recip_lattice)
    n_bvecs = bvectors.n_bvecs

    if mmn
        M, kpb_k_mmn, kpb_b_mmn = read_mmn("$seedname.mmn")

        # check consistency for mmn
        n_bands != size(M)[1] && error("n_bands != size(M)[1]")
        n_bvecs != size(M)[3] && error("n_bvecs != size(M)[3]")
        n_kpts != size(M)[4] && error("n_kpts != size(M)[4]")
        bvectors.kpb_k != kpb_k_mmn && error("kpb_k != kpb_k from mmn file")
        bvectors.kpb_b != kpb_b_mmn && error("kpb_b != kpb_b from mmn file")
    else
        M = zeros(ComplexF64, n_bands, n_bands, n_bvecs, n_kpts)
    end

    if amn
        A = read_amn("$seedname.amn")
        n_bands != size(A)[1] && error("n_bands != size(A)[1]")
        n_wann != size(A)[2] && error("n_wann != size(A)[2]")
        n_kpts != size(A)[3] && error("n_kpts != size(A)[3]")
    else
        A = zeros(ComplexF64, n_bands, n_wann, n_kpts)
        for n = 1:n_wann
            A[n, n, :] .= 1
        end
    end

    if eig
        E = read_eig("$seedname.eig")
        n_bands != size(E)[1] && error("n_bands != size(E)[1]")
        n_kpts != size(E)[2] && error("n_kpts != size(E)[2]")
    else
        E = zeros(Float64, n_bands, n_kpts)
    end

    frozen_bands = falses(n_bands, n_kpts)

    Model(lattice, kgrid, kpoints, bvectors, frozen_bands, M, A, E)
end


function read_nnkp(filename::String)
    @info "Reading nnkp file: $filename"

    io = open(filename)

    read_array(type::Type) = map(x -> parse(type, x), split(readline(io)))

    real_lattice = zeros(Float64, 3, 3)
    recip_lattice = similar(real_lattice)

    n_kpts = nothing
    n_bvecs = nothing
    kpoints = nothing
    kpb_k = nothing
    kpb_b = nothing

    while !eof(io)
        line = readline(io)

        if occursin("begin real_lattice", line)
            for i = 1:3
                real_lattice[:, i] = read_array(Float64)
            end

            line = strip(readline(io))
            line != "end real_lattice" && error("expected end real_lattice")

        elseif occursin("begin recip_lattice", line)
            for i = 1:3
                recip_lattice[:, i] = read_array(Float64)
            end

            line = strip(readline(io))
            line != "end recip_lattice" && error("expected end recip_lattice")

        elseif occursin("begin kpoints", line)
            n_kpts = parse(Int, readline(io))
            kpoints = zeros(Float64, 3, n_kpts)

            for i = 1:n_kpts
                kpoints[:, i] = read_array(Float64)
            end

            line = strip(readline(io))
            line != "end kpoints" && error("expected end kpoints")

        elseif occursin("begin nnkpts", line)
            n_kpts === nothing && error("no kpoints block before nnkpts block?")

            n_bvecs = parse(Int, readline(io))
            kpb_k = zeros(Int, n_bvecs, n_kpts)
            kpb_b = zeros(Int, 3, n_bvecs, n_kpts)

            for ik = 1:n_kpts
                for ib = 1:n_bvecs
                    arr = read_array(Int)
                    ik != arr[1] && error("expected ik = $ik, got $(arr[1])")
                    kpb_k[ib, ik] = arr[2]
                    kpb_b[:, ib, ik] = arr[3:end]
                end
            end

            line = strip(readline(io))
            line != "end nnkpts" && error("expected end nnkpts")
        end
    end
    close(io)

    @info "$filename OK" n_kpts n_bvecs

    bvectors = Matrix{Float64}(undef, 3, n_bvecs)
    fill!(bvectors, NaN)

    weights = zeros(Float64, n_bvecs)
    fill!(weights, NaN)

    BVectors(Mat3{Float64}(recip_lattice), kpoints, bvectors, weights, kpb_k, kpb_b)
end


"""
read si_band.dat  si_band.kpt  si_band.labelinfo.dat
"""
function read_w90_bands(seed_name::String)
    kpaths_coord = Dlm.readdlm("$(seed_name)_band.kpt", Float64; skipstart = 1)
    num_kpts = size(kpaths_coord, 1)

    dat = Dlm.readdlm("$(seed_name)_band.dat", Float64)
    kpaths = reshape(dat[:, 1], num_kpts, :)[:, 1]
    energies = reshape(dat[:, 2], num_kpts, :)
    num_bands = size(energies, 2)

    flabel = open("$(seed_name)_band.labelinfo.dat")
    labels = readlines(flabel)
    close(flabel)

    num_symm_points = length(labels)
    symm_points = Vector{Int}(undef, num_symm_points)
    symm_points_label = Vector{String}(undef, num_symm_points)
    for (i, line) in enumerate(labels)
        lab, idx = split(line)[1:2]
        symm_points[i] = parse(Int, idx)
        symm_points_label[i] = lab
    end

    return Bands(
        num_kpts = num_kpts,
        num_bands = num_bands,
        kpaths = kpaths,
        kpaths_coord = kpaths_coord,
        energies = energies,
        num_symm_points = num_symm_points,
        symm_points = symm_points,
        symm_points_label = symm_points_label,
    )
end

# """
# read si_band.dat  si_band.kpt  si_band.labelinfo.dat
# """
# function write_w90_bands(bands::Bands, seed_name::String)
#     open("$(seed_name)_band.kpt", "w") do io
#         write(io, "$(bands.num_kpts)\n")
#         for i = 1:bands.num_kpts
#             write(io, join([bands.kpaths_coord[:, i]..., 1.0], " "))
#             write(io, "\n")
#         end
#     end
#     println("Written to $(seed_name)_band.kpt")

#     open("$(seed_name)_band.dat", "w") do io
#         for i = 1:bands.num_bands
#             for j = 1:bands.num_kpts
#                 write(io, "$(bands.kpaths[j]) $(bands.energies[j,i])\n")
#             end
#             write(io, "\n")
#         end
#     end
#     println("Written to $(seed_name)_band.dat")

#     open("$(seed_name)_band.labelinfo.dat", "w") do io
#         for i = 1:bands.num_symm_points
#             line = "$(bands.symm_points_label[i]) "
#             idx = bands.symm_points[i]
#             line *= "$(idx) "
#             line *= "$(bands.kpaths[idx]) "
#             coord = join(bands.kpaths_coord[:, idx], " ")
#             line *= "$(coord)\n"
#             write(io, line)
#         end
#     end
#     println("Written to $(seed_name)_band.labelinfo.dat")
# end


function read_wout(filename::String)
    fwout = open(filename)

    ret = Dict()

    start_cell = "Lattice Vectors ("
    start_atom = "|   Site       Fractional Coordinate          Cartesian Coordinate"
    end_atom = "*----------------------------------------------------------------------------*"
    start_finalstate = "Final State"
    end_finalstate = "Sum of centres and spreads"

    unit_is_ang = false
    while !eof(fwout)
        line = strip(readline(fwout))
        if occursin("|  Length Unit", line)
            @assert occursin("Ang", line)
            unit_is_ang = true
            continue
        end
        if occursin(start_cell, line)
            @assert occursin("Ang", line)
            unit_cell = zeros(Float64, 3, 3)
            line = split(strip(readline(fwout)))
            @assert line[1] == "a_1"
            unit_cell[:, 1] = parse.(Float64, line[2:end])
            line = split(strip(readline(fwout)))
            @assert line[1] == "a_2"
            unit_cell[:, 2] = parse.(Float64, line[2:end])
            line = split(strip(readline(fwout)))
            @assert line[1] == "a_3"
            unit_cell[:, 3] = parse.(Float64, line[2:end])
            ret["unit_cell"] = unit_cell
            continue
        end
        if occursin(start_atom, line)
            @assert occursin("Ang", line)
            readline(fwout)
            lines = Vector{String}()
            line = strip(readline(fwout))
            while line != end_atom
                push!(lines, line)
                line = strip(readline(fwout))
            end
            num_atoms = length(lines)
            atoms = zeros(Float64, 3, num_atoms)
            for (i, line) in enumerate(lines)
                line = split(line)
                @assert line[1] == "|" line
                atoms[:, i] = parse.(Float64, line[8:10])
            end
            ret["atoms"] = atoms
            continue
        end
        if occursin(start_finalstate, line)
            lines = Vector{String}()
            line = strip(readline(fwout))
            while !occursin(end_finalstate, line)
                push!(lines, line)
                line = strip(readline(fwout))
            end
            num_wann = length(lines)
            wf_centers = zeros(Float64, 3, num_wann)
            wf_spreads = zeros(Float64, num_wann)
            for (i, line) in enumerate(lines)
                line = split(line)
                @assert join(line[1:4], " ") == "WF centre and spread"
                @assert i == parse(Int, line[5])
                x = parse(Float64, replace(line[7], "," => ""))
                y = parse(Float64, replace(line[8], "," => ""))
                z = parse(Float64, replace(line[9], "," => ""))
                s = parse(Float64, line[11])
                wf_centers[:, i] = [x, y, z]
                wf_spreads[i] = s
            end
            ret["wf_centers"] = wf_centers
            ret["wf_spreads"] = wf_spreads
            continue
        end
    end
    close(fwout)
    @assert unit_is_ang

    return ret
end

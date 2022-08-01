using Printf: @printf

function read_win(filename::String)
    @info "Reading win file: $filename"
    io = open(filename)

    # The win file uses "num_wann", so I keep it as is, and not using "n_wann".
    num_wann = missing
    num_bands = missing
    mp_grid = missing
    unit_cell = missing
    kpoints = missing
    kpoint_path = missing
    kmesh_tol = missing
    dis_froz_min = missing
    dis_froz_max = missing
    dis_win_min = missing
    dis_win_max = missing
    atoms_frac = missing
    atoms_cart = missing
    atom_labels = missing

    # handle case insensitive win files (relic of Fortran)
    read_lowercase_line() = strip(lowercase(readline(io)))
    parse_array(line::AbstractString) = map(x -> parse(Float64, x), split(line))
    read_array(f::IOStream) = parse_array(readline(f))

    while !eof(io)
        line = read_lowercase_line()
        line = replace(line, "=" => " ", ":" => " ", "," => " ")

        i = findfirst(r"!|#", line)
        if i != nothing
            line = strip(line[1:(i.start - 1)])
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
        elseif occursin("kmesh_tol", line)
            kmesh_tol = parse_float(split(line)[2])
        elseif occursin("dis_froz_min", line)
            dis_froz_min = parse_float(split(line)[2])
        elseif occursin("dis_froz_max", line)
            dis_froz_max = parse_float(split(line)[2])
        elseif occursin("dis_win_min", line)
            dis_win_min = parse_float(split(line)[2])
        elseif occursin("dis_win_max", line)
            dis_win_max = parse_float(split(line)[2])
        elseif occursin(r"begin\s+unit_cell_cart", line)
            unit_cell = zeros(Float64, 3, 3)
            line = read_lowercase_line()
            unit = line
            if !startswith(unit, "b") && !startswith(unit, "a")
                unit = "ang"
            else
                line = read_lowercase_line()
            end
            for i in 1:3
                # in win file, each line is a lattice vector, here it is stored as column vec
                unit_cell[:, i] = parse_array(line)
                line = read_lowercase_line()
            end
            if startswith(unit, "b")
                # convert to angstrom
                unit_cell .*= Bohr
            end
            unit_cell = Mat3{Float64}(unit_cell)
        elseif occursin(r"begin\s+atoms_frac", line)
            line = strip(readline(io))
            # I need to read all lines and get n_atoms
            lines = Vector{String}()
            while !occursin(r"end\s+atoms_frac", lowercase(line))
                push!(lines, line)
                # do not lowercase due to atomic label
                line = strip(readline(io))
            end
            n_atoms = length(lines)
            atoms_frac = zeros(Float64, 3, n_atoms)
            atom_labels = Vector{String}()
            for i in 1:n_atoms
                l = split(lines[i])
                push!(atom_labels, l[1])
                atoms_frac[:, i] = parse_float.(l[2:end])
            end
        elseif occursin(r"begin\s+atoms_cart", line)
            line = strip(readline(io))
            unit = lowercase(line)
            if !startswith(unit, "b") && !startswith(unit, "a")
                unit = "ang"
            else
                # do not lowercase due to atomic label
                line = strip(readline(io))
            end
            # I need to read all lines and get n_atoms
            lines = Vector{String}()
            while !occursin(r"end\s+atoms_cart", lowercase(line))
                push!(lines, line)
                line = strip(readline(io))
            end
            n_atoms = length(lines)
            atoms_cart = zeros(Float64, 3, n_atoms)
            atom_labels = Vector{String}()
            for i in 1:n_atoms
                l = split(lines[i])
                push!(atom_labels, l[1])
                atoms_cart[:, i] = parse_float.(l[2:end])
            end
            if startswith(unit, "b")
                # convert to angstrom
                atoms_cart .*= Bohr
            end
        elseif occursin(r"begin\s+kpoints", line)
            line = read_lowercase_line()
            # kpoints block might be defined before mp_grid!
            # I need to read all lines and get n_kpts
            lines = Vector{String}()
            while !occursin(r"end\s+kpoints", line)
                push!(lines, line)
                line = read_lowercase_line()
            end

            n_kpts = length(lines)
            kpoints = Matrix{Float64}(undef, 3, n_kpts)
            for i in 1:n_kpts
                # There might be weight at 4th column, but we don't use it.
                kpoints[:, i] = parse_array(lines[i])[1:3]
            end
        elseif occursin(r"begin\skpoint_path", line)
            kpoint_path = Kpath{Float64}()

            # allow uppercase
            line = strip(readline(io))
            while !occursin(r"end\s+kpoint_path", lowercase(line))
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

    # if atoms_cart, convert to fractional
    if ismissing(atoms_frac)
        ismissing(atoms_cart) && error("both atoms_frac and atoms_cart are missing")
        atoms_frac = inv(unit_cell) * atoms_cart
    end

    println("  num_wann  = ", num_wann)
    println("  num_bands = ", num_bands)
    @printf("  mp_grid   = %d %d %d\n", mp_grid...)
    println()

    return (
        num_wann=num_wann,
        num_bands=num_bands,
        mp_grid=mp_grid,
        kpoints=kpoints,
        kmesh_tol=kmesh_tol,
        unit_cell=unit_cell,
        atoms_frac=atoms_frac,
        atom_labels=atom_labels,
        kpoint_path=kpoint_path,
        dis_froz_min=dis_froz_min,
        dis_froz_max=dis_froz_max,
        dis_win_min=dis_win_min,
        dis_win_max=dis_win_max,
    )
end

"""
Parse wout file.

Return lattice in angstrom (each column is a lattice vector),
atom positions in fractional coordinates (each column is a coordinate),
WF centers in angstrom, and WF spreads in angstrom^2.
"""
function read_wout(filename::String)
    io = open(filename)

    start_cell = "Lattice Vectors ("
    start_atom = "|   Site       Fractional Coordinate          Cartesian Coordinate"
    end_atom = "*----------------------------------------------------------------------------*"
    start_finalstate = "Final State"
    end_finalstate = "Sum of centres and spreads"

    ang_unit = false
    lattice = nothing
    atom_labels = nothing
    atom_positions = nothing
    centers = nothing
    spreads = nothing

    while !eof(io)
        line = strip(readline(io))

        if occursin("|  Length Unit", line)
            @assert occursin("Ang", line)
            ang_unit = true
            continue
        end

        if occursin(start_cell, line)
            @assert occursin("Ang", line)
            lattice = zeros(Float64, 3, 3)

            line = split(strip(readline(io)))
            @assert line[1] == "a_1"
            lattice[:, 1] = parse.(Float64, line[2:end])

            line = split(strip(readline(io)))
            @assert line[1] == "a_2"
            lattice[:, 2] = parse.(Float64, line[2:end])

            line = split(strip(readline(io)))
            @assert line[1] == "a_3"
            lattice[:, 3] = parse.(Float64, line[2:end])

            continue
        end

        if occursin(start_atom, line)
            @assert occursin("Ang", line)
            readline(io)

            lines = Vector{String}()
            line = strip(readline(io))
            while line != end_atom
                push!(lines, line)
                line = strip(readline(io))
            end

            n_atom = length(lines)
            atom_labels = Vector{String}()
            atom_positions = zeros(Float64, 3, n_atom)
            for (i, line) in enumerate(lines)
                line = split(line)
                @assert line[1] == "|" line
                push!(atom_labels, line[2])
                # cartesian
                # atom_positions[:, i] = parse.(Float64, line[8:10])
                # fractional
                atom_positions[:, i] = parse.(Float64, line[4:6])
            end

            continue
        end

        if occursin(start_finalstate, line)
            lines = Vector{String}()
            line = strip(readline(io))
            while !occursin(end_finalstate, line)
                push!(lines, line)
                line = strip(readline(io))
            end

            n_wann = length(lines)
            centers = zeros(Float64, 3, n_wann)
            spreads = zeros(Float64, n_wann)
            for (i, line) in enumerate(lines)
                line = split(line)
                @assert join(line[1:4], " ") == "WF centre and spread"
                @assert i == parse(Int, line[5])

                x = parse(Float64, replace(line[7], "," => ""))
                y = parse(Float64, replace(line[8], "," => ""))
                z = parse(Float64, replace(line[9], "," => ""))
                s = parse(Float64, line[11])

                centers[:, i] = [x, y, z]
                spreads[i] = s
            end

            continue
        end
    end

    close(io)

    @assert ang_unit

    return (
        lattice=lattice,
        atom_labels=atom_labels,
        atom_positions=atom_positions,
        centers=centers,
        spreads=spreads,
    )
end

module InOut

using Base: String
using Dates: now
import DelimitedFiles as Dlm
import LinearAlgebra as LA

using ..Constants: Bohr
using ..Parameters: CoreData, InputParams, Bands, AtomicWavefunction, Projectabilities
using ..Utilities: get_recipcell
using ..BVectors: generate_bvectors

function read_win(filename::String)::Dict
    println("Reading $filename")
    fwin = open(filename)

    tol = 1e-6
    num_wann = missing
    num_bands = missing
    num_kpts = missing
    kpts_size = zeros(3)
    unit_cell = zeros(3, 3)
    kpts = zeros(3, 1)
    # 1st index: coordinates, 2nd index: start & end, 3rd index: number of paths. Need to resize the 3rd index
    kpath = Array{Float64, 3}(undef, 3, 2, 1)
    kpath_label = Array{String, 2}(undef, 2, 1)

    read_array(f) = map(x -> parse(Float64, x), split(readline(f)))

    while !eof(fwin)
        line = readline(fwin)
        # handle case insensitive win files (relic of Fortran)
        line = strip(lowercase(line))
        line = replace(line, "=" => " ")
        line = replace(line, ":" => " ")
        line = replace(line, "," => " ")
        i = findfirst(isequal('#'), line)
        if i !== nothing
            line = strip(line[1:i-1])
        end
        if startswith(line, r"!|#") || isempty(line)
            continue
        elseif occursin("mp_grid", line)
            kpts_size = map(x -> parse(Int, x), split(line)[2:4])
            num_kpts = prod(kpts_size)
            kpts = zeros(3, num_kpts)
        elseif occursin("num_bands", line)
            num_bands = parse(Int, split(line)[2])
        elseif occursin("num_wann", line)
            num_wann = parse(Int, split(line)[2])
        elseif occursin("begin unit_cell_cart", line)
            unit = strip(lowercase(readline(fwin)))
            for i = 1:3
                # in win file, each line is a lattice vector, here it is stored as column vec
                unit_cell[:, i] = read_array(fwin)
            end
            if startswith(unit, r"b")
                # convert to angstrom
                unit_cell .*= Bohr
            end
        elseif occursin("begin kpoints", line)
            for i = 1:num_kpts
                kpts[:, i] = read_array(fwin)[1:3]
            end
        elseif occursin("begin kpoint_path", line)
            lines = Vector{String}()
            line = strip(readline(fwin))
            while !occursin("end kpoint_path", line)
                push!(lines, line)
                line = strip(readline(fwin))
            end
            num_kpath = length(lines)
            kpath = Array{Float64, 3}(undef, 3, 2, num_kpath)
            kpath_label = Array{String, 2}(undef, 2, num_kpath)
            for i = 1:num_kpath
                l = split(lines[i])
                @assert length(l) == 8
                kpath_label[1,i] = l[1]
                kpath_label[2,i] = l[5]
                kpath[:,1,i] = parse.(Float64, l[2:4])
                @assert all(-1.0-tol .<= kpath[:,1,i] .<= 1.0+tol)
                kpath[:,2,i] = parse.(Float64, l[6:8])
                @assert all(-1.0-tol .<= kpath[:,2,i] .<= 1.0+tol)
            end
        end
    end
    close(fwin)

    @assert !ismissing(num_wann)
    @assert num_wann > 0

    if ismissing(num_bands)
        num_bands = num_wann
    end
    @assert num_bands > 0

    @assert all(i -> i > 0, num_kpts)

    @assert all(x -> !ismissing(x), unit_cell)

    recip_cell = get_recipcell(unit_cell)

    println("$filename OK, num_wann = $num_wann, num_bands = $num_bands, num_kpts = $num_kpts")

    return Dict(
        "num_wann" => num_wann,
        "num_bands" => num_bands,
        "num_kpts" => num_kpts,
        "kpts_size" => kpts_size,
        "kpts" => kpts,
        "unit_cell" => unit_cell,
        "recip_cell" => recip_cell,
        "kpath" => kpath,
        "kpath_label" => kpath_label
    )
end

function read_mmn(filename::String)::Tuple
    println("Reading $filename")
    fmmn = open(filename)

    # skip header
    readline(fmmn)
    line = readline(fmmn)
    num_bands, num_kpts, num_bvecs = map(x -> parse(Int64, x), split(line))

    # overlap matrix
    mmn = zeros(ComplexF64, num_bands, num_bands, num_bvecs, num_kpts)
    # for each point, list of neighbors, (K) representation
    bvecs = zeros(Int64, num_bvecs, num_kpts)
    bvecs_disp = zeros(Int64, 3, num_bvecs, num_kpts)

    while !eof(fmmn)
        for ib = 1:num_bvecs
            line = readline(fmmn)
            arr = split(line)
            k = parse(Int64, arr[1])
            kpb = parse(Int64, arr[2])
            bvecs_disp[:, ib, k] = map(x -> parse(Int64, x), arr[3:5])
            bvecs[ib, k] = kpb
            for n = 1:num_bands
                for m = 1:num_bands
                    line = readline(fmmn)
                    arr = split(line)
                    o = parse(Float64, arr[1]) + im * parse(Float64, arr[2])
                    mmn[m, n, ib, k] = o
                    @assert !isnan(o)
                end
            end
        end
    end
    close(fmmn)
    println("$filename OK, size = ", size(mmn))

    return mmn, bvecs, bvecs_disp
end

function read_amn(filename::String)
    println("Reading $filename")

    famn = open("$filename")
    readline(famn)
    arr = split(readline(famn))
    num_bands = parse(Int64, arr[1])
    num_kpts = parse(Int64, arr[2])
    num_wann = parse(Int64, arr[3])
    amn = zeros(ComplexF64, num_bands, num_wann, num_kpts)

    while !eof(famn)
        line = readline(famn)
        arr = split(line)
        m = parse(Int64, arr[1])
        n = parse(Int64, arr[2])
        k = parse(Int64, arr[3])
        a = parse(Float64, arr[4]) + im * parse(Float64, arr[5])
        amn[m, n, k] = a
    end
    close(famn)

    println("$filename OK, size = ", size(amn))
    return amn
end

function read_eig(filename::String)
    println("Reading $filename")

    feig = open(filename)
    lines = readlines(filename)
    close(feig)

    len = length(lines)
    indexb = zeros(len)
    indexk = zeros(len)
    eig = zeros(len)

    for i = 1:len
        arr = split(lines[i])
        indexb[i] = parse(Int, arr[1])
        indexk[i] = parse(Int, arr[2])
        eig[i] = parse(Float64, arr[3])
    end
    
    # find unique elements
    num_bands = length(Set(indexb))
    num_kpts = length(Set(indexk))
    eig = reshape(eig, (num_bands, num_kpts))

    # check eigenvalues should be in order
    # some times there are small noise
    round_digits(x) = round(x, digits=9)
    for ik = 1:num_kpts
        # @assert issorted(eig[:, ik]) display(ik) display(eig[:, ik])
        @assert issorted(eig[:, ik], by=round_digits) display(ik) display(eig[:, ik])
    end
    
    println("$filename OK, size = ", size(eig))
    return eig
end

function read_seedname(seedname::String; amn::Bool=true, mmn::Bool=true, eig::Bool=true)
    # read win, mmn and optionally amn
    win = read_win("$seedname.win")
    num_bands = win["num_bands"]
    num_wann = win["num_wann"]
    num_kpts = win["num_kpts"]

    kpbs, kpbs_disp, kpbs_weight = generate_bvectors(win["kpts"], win["recip_cell"])
    num_bvecs = size(kpbs_weight, 1)

    if mmn
        m_mat, kpbs2, kpbs_disp2 = read_mmn("$seedname.mmn")
        # check consistency for mmn
        @assert num_bands == size(m_mat)[1]
        @assert num_bvecs == size(m_mat)[3]
        @assert num_kpts == size(m_mat)[4]
        @assert kpbs == kpbs2
        @assert kpbs_disp == kpbs_disp2
    else
        m_mat = zeros(ComplexF64, num_bands, num_bands, num_bvecs, num_kpts)
    end

    if amn
        a_mat = read_amn("$seedname.amn")
        @assert num_bands == size(a_mat)[1]
        @assert num_wann == size(a_mat)[2]
        @assert num_kpts == size(a_mat)[3]
    else
        # TODO: not tested
        a_mat = zeros(ComplexF64, num_bands, num_wann, num_kpts)
        for n = 1:num_wann
            a_mat[n, n, :] .= 1
        end
    end

    if eig
        eig_mat = read_eig("$seedname.eig")
        @assert num_bands == size(eig_mat)[1]
        @assert num_kpts == size(eig_mat)[2]
    else
        eig_mat = zeros(num_bands, num_kpts)
    end

    # FIX: not tested
    # map = true
    # logMethod = true

    println("num_bands = $num_bands, num_wann = $num_wann, kpt = ", 
    win["kpts_size"], " num_bvecs = $num_bvecs")

    frozen = falses(num_bands, num_kpts)

    return CoreData(
        unit_cell = win["unit_cell"],
        recip_cell = win["recip_cell"],
        num_bands = win["num_bands"],
        num_wann = win["num_wann"],
        num_kpts = win["num_kpts"],
        kpts_size = win["kpts_size"],
        kpts = win["kpts"],
        num_bvecs = num_bvecs,
        kpbs = kpbs,
        kpbs_weight = kpbs_weight,
        kpbs_disp = kpbs_disp,
        frozen = frozen,
        mmn = m_mat,
        amn = a_mat,
        eig = eig_mat
    )
end

function read_nnkp(filename::String)
    println("Reading $filename")

    fnnkp = open(filename)

    read_array(f, type) = map(x -> parse(type, x), split(readline(f)))

    real_lattice = zeros(3, 3)
    recip_lattice = zeros(3, 3)
    num_kpts = 0
    num_bvecs = 0
    kpoints = zeros(3, 1)
    nnkpts = zeros(Int, 4, 1, 1)

    while !eof(fnnkp)
        line = readline(fnnkp)
        if occursin("begin real_lattice", line)
            for i = 1:3
                real_lattice[:, i] = read_array(fnnkp, Float64)
            end
            line = strip(readline(fnnkp))
            @assert line == "end real_lattice"
        elseif occursin("begin recip_lattice", line)
            for i = 1:3
                recip_lattice[:, i] = read_array(fnnkp, Float64)
            end
            line = strip(readline(fnnkp))
            @assert line == "end recip_lattice"
        elseif occursin("begin kpoints", line)
            num_kpts = parse(Int, readline(fnnkp))
            kpoints = zeros(3, num_kpts)
            for i = 1:num_kpts
                kpoints[:, i] = read_array(fnnkp, Float64)
            end
            line = strip(readline(fnnkp))
            @assert line == "end kpoints"
        elseif occursin("begin nnkpts", line)
            num_bvecs = parse(Int, readline(fnnkp))
            nnkpts = zeros(Int, 4, num_bvecs, num_kpts)
            for _ = 1:num_kpts
                for b = 1:num_bvecs
                    arr = read_array(fnnkp, Int)
                    k = arr[1]
                    nnkpts[:, b, k] = arr[2:end]
                end
            end
            line = strip(readline(fnnkp))
            @assert line == "end nnkpts"
        end
    end
    close(fnnkp)

    println("$filename OK, num_kpts = $num_kpts, num_bvecs = $num_bvecs")
    return Wannier90Nnkp(
        real_lattice,
        recip_lattice,
        num_kpts,
        kpoints,
        num_bvecs,
        nnkpts
    )
end

"""
Output amn file
"""
function write_amn(filename::String, A::Array{ComplexF64,3})
    famn = open(filename, "w")
    write(famn, "Created by Wannier.jl ", string(now()), "\n")
    num_bands, num_wann, num_kpts = size(A)
    write(famn, "$num_bands $num_kpts $num_wann\n")
    for ik = 1:num_kpts
        for iw = 1:num_wann
            for ib = 1:num_bands
                coeff = A[ib,iw,ik]
                write(famn, "$ib $iw $ik $(real(coeff)) $(imag(coeff))\n")
            end
        end
    end
    close(famn)
    println("Written to file $(filename)")
end

"""
read si_band.dat  si_band.kpt  si_band.labelinfo.dat
"""
function read_wannier90_bands(seed_name::String)
    kpaths_coord = Dlm.readdlm("$(seed_name)_band.kpt", Float64; skipstart=1)
    num_kpts = size(kpaths_coord, 1)

    dat = Dlm.readdlm("$(seed_name)_band.dat", Float64)
    kpaths = reshape(dat[:,1], num_kpts, :)[:,1]
    energies = reshape(dat[:,2], num_kpts, :)
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
        symm_points_label = symm_points_label
    )
end

"""
read si_band.dat  si_band.kpt  si_band.labelinfo.dat
"""
function write_wannier90_bands(bands::Bands, seed_name::String)
    open("$(seed_name)_band.kpt", "w") do io
        write(io, "$(bands.num_kpts)\n")
        for i = 1:bands.num_kpts
            write(io, join([bands.kpaths_coord[:,i]..., 1.0], " "))
            write(io, "\n")
        end
    end
    println("Written to $(seed_name)_band.kpt")

    open("$(seed_name)_band.dat", "w") do io
        for i = 1:bands.num_bands
            for j = 1:bands.num_kpts
                write(io, "$(bands.kpaths[j]) $(bands.energies[j,i])\n")
            end
            write(io, "\n")
        end
    end
    println("Written to $(seed_name)_band.dat")

    open("$(seed_name)_band.labelinfo.dat", "w") do io
        for i = 1:bands.num_symm_points
            line = "$(bands.symm_points_label[i]) "
            idx = bands.symm_points[i]
            line *= "$(idx) "
            line *= "$(bands.kpaths[idx]) "
            coord = join(bands.kpaths_coord[:,idx], " ")
            line *= "$(coord)\n"
            write(io, line)
        end
    end
    println("Written to $(seed_name)_band.labelinfo.dat")
end

"""
read QE bands.x output dat file

```
 &plot nbnd=  20, nks=   380 /
           -0.500000  0.500000  0.500000
   -3.320   -0.666    5.173    5.173    7.994    9.725    9.725   14.147   16.993   16.993
   17.841   17.841   17.902   19.666   25.961   26.563   28.186   28.186   28.368   28.368
           -0.495000  0.495000  0.495000
   -3.322   -0.664    5.173    5.173    7.994    9.725    9.725   14.148   16.980   16.980
```
"""
function read_qe_bands(filename::String)
    fdat = open(filename)

    line = readline(fdat)
    regex = r"&plot nbnd=\s*(\d+), nks=\s*(\d+) /"
    m = match(regex, line)
    if m !== nothing
        nbnd, nks = parse.(Int, m.captures)
    else
        regex = r"&plot nbnd=\s*(\d+), nks=\s*(\d+) alat=\s*([+-]?([0-9]*[.])?[0-9]+) /"
        m = match(regex, line)
        nbnd, nks = parse.(Int, m.captures[1:2])
        alat = parse.(Float64, m.captures[3])
    end

    kpaths_coord = Matrix{Float64}(undef, 3, nks)
    energies = Matrix{Float64}(undef, nks, nbnd)

    for ik = 1:nks
        # QE kpt are in absolute coordinates, but is scaled by `alat`
        kpaths_coord[:,ik] = parse.(Float64, split(readline(fdat)))
        ib = 0
        while ib < nbnd
            e = parse.(Float64, split(readline(fdat)))
            len_e = length(e)
            energies[ik, ib+1:ib+len_e] = e
            ib += len_e
        end
        @assert ib == nbnd
    end

    @assert eof(fdat)

    # detect high symmetry points - by the criteria that 
    # there are angles between two consecutive kpath
    symm_points = Vector{Int}()
    symm_points_label = Vector{String}()
    # push the first kpt
    push!(symm_points, 1)
    push!(symm_points_label, "")
    for ik = 2:nks-1
        vec0 = kpaths_coord[:,ik] - kpaths_coord[:,ik-1]
        vec1 = kpaths_coord[:,ik+1] - kpaths_coord[:,ik]
        if !all(isapprox.(LA.cross(vec0, vec1), 0; atol=2e-6))
            push!(symm_points, ik)
            push!(symm_points_label, "")
        end
    end
    # push the last kpt
    push!(symm_points, nks)
    push!(symm_points_label, "")
    num_symm_points = length(symm_points)

    # generate kpath
    kpaths = zeros(Float64, nks)
    ik_prev = 1
    for ik = 1:nks
        dk = LA.norm(kpaths_coord[:,ik] - kpaths_coord[:,ik_prev])
        if ik != 2 && ik_prev in symm_points
            dk = 0
        end
        kpaths[ik] = kpaths[ik_prev] + dk
        ik_prev = ik
    end
    
    return Bands(nks, nbnd, kpaths, kpaths_coord, energies, 
    num_symm_points, symm_points, symm_points_label)
end

function read_qe_projwfcup(filename::String)
    fdat = open(filename)

    splitline() = split(strip(readline(fdat)))

    # header
    title = strip(readline(fdat), '\n')
    nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp = parse.(Int, splitline())
    line = splitline()
    ibrav = parse(Int, line[1])
    celldm = parse.(Float64, line[2:end])
    line = splitline()
    # some version of projwfc.x output the unit_cell
    if length(line) == 3
        readline(fdat)
        readline(fdat)
        line = splitline()
    end
    gcutm, dual, ecutwfc = parse.(Float64, line[1:end-1])
    magicnum = parse(Int, line[end])
    @assert magicnum == 9
    atm = Vector{String}(undef, ntyp)
    zv = zeros(Float64, ntyp)
    for i = 1:ntyp
        line = splitline()
        nt = parse(Int, line[1])
        @assert nt == i
        atm[i] = line[2]
        zv[i] = parse(Float64, line[3])
    end
    tau = zeros(Float64, 3, nat)
    ityp = zeros(Int, nat)
    for i = 1:nat
        line = splitline()
        na = parse(Int, line[1])
        @assert na == i
        tau[:, i] = parse.(Float64, line[2:4])
        ityp[i] = parse(Int, line[5])
    end
    natomwfc, nkstot, nbnd = parse.(Int, splitline())
    parsebool(s::Union{String,SubString}) = lowercase(s) == "t" ? true : false
    noncolin, lspinorb = parsebool.(splitline())
    @assert !noncolin && !lspinorb

    # projection data
    nlmchi = Vector{Dict}()
    proj = zeros(Float64, nkstot, nbnd, natomwfc)
    for iw = 1:natomwfc
        line = splitline()
        nwfc = parse(Int, line[1])
        @assert nwfc == iw
        na = parse(Int, line[2])
        atm_name = line[3]
        @assert atm_name == atm[ityp[na]]
        els = line[4]
        n, l, m = parse.(Int, line[5:end])
        push!(nlmchi, Dict("na"=>na, "els"=>els, "n"=>n, "l"=>l, "m"=>m))
        for ik = 1:nkstot
            for ib = 1:nbnd
                line = splitline()
                k, b = parse.(Int, line[1:2])
                @assert k == ik && b == ib
                p = parse(Float64, line[3])
                proj[ik,ib,iw] = p
            end
        end
    end

    wfcs_type = Vector{AtomicWavefunction}(undef, natomwfc)
    for i=1:natomwfc
        atom_index = nlmchi[i]["na"]
        atom_label = atm[ityp[atom_index]]
        wfc_label = nlmchi[i]["els"]
        n = nlmchi[i]["n"]
        l = nlmchi[i]["l"]
        m = nlmchi[i]["m"]
        wfc = AtomicWavefunction(atom_index, atom_label, wfc_label, n, l, m)
        wfcs_type[i] = wfc
    end

    return Projectabilities(
        nkstot,
        nbnd,
        natomwfc,
        wfcs_type,
        proj
    )
end

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
            unit_cell[:,1] = parse.(Float64, line[2:end])
            line = split(strip(readline(fwout)))
            @assert line[1] == "a_2"
            unit_cell[:,2] = parse.(Float64, line[2:end])
            line = split(strip(readline(fwout)))
            @assert line[1] == "a_3"
            unit_cell[:,3] = parse.(Float64, line[2:end])
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
            for (i,line) in enumerate(lines)
                line = split(line)
                @assert line[1] == "|" line
                atoms[:,i] = parse.(Float64, line[8:10])
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
            for (i,line) in enumerate(lines)
                line = split(line)
                @assert join(line[1:4], " ") == "WF centre and spread"
                @assert i == parse(Int, line[5])
                x = parse(Float64, replace(line[7], "," => ""))
                y = parse(Float64, replace(line[8], "," => ""))
                z = parse(Float64, replace(line[9], "," => ""))
                s = parse(Float64, line[11])
                wf_centers[:,i] = [x,y,z]
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

end

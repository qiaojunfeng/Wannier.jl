module InOut

using Dates: now

using ..Constants: Bohr
using ..Parameters: InputParams
using ..Utilities: get_recipcell
using ..BVectors: generate_bvectors

function read_win(filename::String)::Dict
    println("Reading $filename")
    fwin = open(filename)

    num_wann = missing
    num_bands = missing
    num_kpts = missing
    kpts_size = zeros(3)
    unit_cell = zeros(3, 3)
    kpts = zeros(3, 1)

    read_array(f) = map(x -> parse(Float64, x), split(readline(f)))

    while !eof(fwin)
        line = readline(fwin)
        # handle case insensitive win files (relic of Fortran)
        line = strip(lowercase(line))
        line = replace(line, "=" => " ")
        line = replace(line, ":" => " ")
        line = replace(line, "," => " ")
        if startswith(line, r"!|#")
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
        "recip_cell" => recip_cell
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

    # FIX: normalization should be done later
    # for k = 1:num_kpts
    #     amn[:, :, k] = orthonormalize_lowdin(amn[:, :, k])
    # end

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

function read_seedname(seedname::String, amn::Bool=true, eig::Bool=true)
    # read win, mmn and optionally amn
    win = read_win("$seedname.win")
    num_bands = win["num_bands"]
    num_wann = win["num_wann"]
    num_kpts = win["num_kpts"]

    kpbs, kpbs_disp, kpbs_weight = generate_bvectors(win["kpts"], win["recip_cell"])
    num_bvecs = size(kpbs_weight, 1)

    m_mat, kpbs2, kpbs_disp2 = read_mmn("$seedname.mmn")
    # check consistency for mmn
    @assert num_bands == size(m_mat)[1]
    @assert num_bvecs == size(m_mat)[3]
    @assert num_kpts == size(m_mat)[4]
    @assert kpbs == kpbs2
    @assert kpbs_disp == kpbs_disp2

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

    return InputParams(
        seedname,
        win["unit_cell"],
        win["recip_cell"],
        win["num_bands"],
        win["num_wann"],
        win["num_kpts"],
        win["kpts_size"],
        win["kpts"],
        num_bvecs,
        kpbs,
        kpbs_weight,
        kpbs_disp,
        m_mat,
        a_mat,
        eig_mat
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

end

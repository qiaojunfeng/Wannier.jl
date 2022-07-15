using LinearAlgebra

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

    for ik in 1:nks
        # QE kpt are in absolute coordinates, but is scaled by `alat`
        kpaths_coord[:, ik] = parse.(Float64, split(readline(fdat)))
        ib = 0
        while ib < nbnd
            e = parse.(Float64, split(readline(fdat)))
            len_e = length(e)
            energies[ik, (ib + 1):(ib + len_e)] = e
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
    for ik in 2:(nks - 1)
        vec0 = kpaths_coord[:, ik] - kpaths_coord[:, ik - 1]
        vec1 = kpaths_coord[:, ik + 1] - kpaths_coord[:, ik]
        if !all(isapprox.(cross(vec0, vec1), 0; atol=2e-6))
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
    for ik in 1:nks
        dk = norm(kpaths_coord[:, ik] - kpaths_coord[:, ik_prev])
        if ik != 2 && ik_prev in symm_points
            dk = 0
        end
        kpaths[ik] = kpaths[ik_prev] + dk
        ik_prev = ik
    end

    return Bands(
        nks,
        nbnd,
        kpaths,
        kpaths_coord,
        energies,
        num_symm_points,
        symm_points,
        symm_points_label,
    )
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
    gcutm, dual, ecutwfc = parse.(Float64, line[1:(end - 1)])
    magicnum = parse(Int, line[end])
    @assert magicnum == 9
    atm = Vector{String}(undef, ntyp)
    zv = zeros(Float64, ntyp)
    for i in 1:ntyp
        line = splitline()
        nt = parse(Int, line[1])
        @assert nt == i
        atm[i] = line[2]
        zv[i] = parse(Float64, line[3])
    end
    tau = zeros(Float64, 3, nat)
    ityp = zeros(Int, nat)
    for i in 1:nat
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
    for iw in 1:natomwfc
        line = splitline()
        nwfc = parse(Int, line[1])
        @assert nwfc == iw
        na = parse(Int, line[2])
        atm_name = line[3]
        @assert atm_name == atm[ityp[na]]
        els = line[4]
        n, l, m = parse.(Int, line[5:end])
        push!(nlmchi, Dict("na" => na, "els" => els, "n" => n, "l" => l, "m" => m))
        for ik in 1:nkstot
            for ib in 1:nbnd
                line = splitline()
                k, b = parse.(Int, line[1:2])
                @assert k == ik && b == ib
                p = parse(Float64, line[3])
                proj[ik, ib, iw] = p
            end
        end
    end

    wfcs_type = Vector{AtomicWavefunction}(undef, natomwfc)
    for i in 1:natomwfc
        atom_index = nlmchi[i]["na"]
        atom_label = atm[ityp[atom_index]]
        wfc_label = nlmchi[i]["els"]
        n = nlmchi[i]["n"]
        l = nlmchi[i]["l"]
        m = nlmchi[i]["m"]
        wfc = AtomicWavefunction(atom_index, atom_label, wfc_label, n, l, m)
        wfcs_type[i] = wfc
    end

    return Projectabilities(nkstot, nbnd, natomwfc, wfcs_type, proj)
end

@testitem "generate kpath from win" begin
    using Wannier: SymbolVec3
    using Brillouin: LATTICE

    lattice = [
        -2.6988 0.0 -2.6988
        0.0 2.6988 2.6988
        2.6988 2.6988 0.0
    ]
    kpoint_path = [
        [:L => [0.5, 0.5, 0.5], :G => [0.0, 0.0, 0.0]],
        [:G => [0.0, 0.0, 0.0], :X => [0.5, 0.0, 0.5]],
        [:X => [0.5, -0.5, 0.0], :K => [0.375, -0.375, 0.0]],
        [:K => [0.375, -0.375, 0.0], :G => [0.0, 0.0, 0.0]],
    ]
    win_kpath = generate_kpath(lattice, kpoint_path)

    ref_kpath = Dict(
        "basis" => [
            [-1.164069, -1.164069, 1.164069],
            [1.164069, 1.164069, 1.164069],
            [-1.164069, 1.164069, -1.164069],
        ],
        "paths" => [[:L, :G, :X], [:X_1, :K, :G]],
        "points" => Dict(
            :X_1 => [0.5, -0.5, 0.0],
            :G => [0.0, 0.0, 0.0],
            :K => [0.375, -0.375, 0.0],
            :L => [0.5, 0.5, 0.5],
            :X => [0.5, 0.0, 0.5],
        ),
        "setting" => Ref(LATTICE),
    )
    @test ref_kpath["points"] == win_kpath.points
    @test ref_kpath["paths"] == win_kpath.paths
    @test isapprox(ref_kpath["basis"], win_kpath.basis; atol=1e-5)
    @test Symbol(ref_kpath["setting"]) == Symbol(win_kpath.setting)
end

@testitem "generate w90 bands kpi" begin
    using Wannier.Datasets
    win = read_win(dataset"Si2/Si2.win")
    lattice = win.unit_cell_cart
    recip_lattice = reciprocal_lattice(lattice)
    ref_kpi, _ = read_w90_band(dataset"Si2/reference/Si2", recip_lattice)

    # number of points along the 1st segment
    n_points = 100
    kpath = generate_kpath(lattice, win.kpoint_path)
    kpi = generate_w90_kpoint_path(kpath, n_points)

    @test all(isapprox.(kpi.kpaths, ref_kpi.kpaths; atol=1e-5))
    # If in the kpath block of win file, there are two kpoints with same label but
    # different coordinates, I will append a number to the repeated label in read_win,
    # so I only compare label without number.
    labels = map(kpi.labels) do line
        new_line = typeof(line)()
        for (idx, lab) in line
            slab = string(lab)
            delim = "_"
            if occursin(delim, slab)
                # only remove the last number suffix
                slab = join(split(slab, delim)[1:(end - 1)], delim)
            end
            push!(new_line, idx => Symbol(slab))
        end
        new_line
    end
    @test labels == ref_kpi.labels
    @test kpi.basis ≈ ref_kpi.basis
    @test Symbol(kpi.setting) == Symbol(ref_kpi.setting)
end

@testitem "get_linear_path" begin
    using WannierIO
    using Wannier.Datasets
    win = read_win(dataset"Si2/Si2.win")
    lattice = win.unit_cell_cart
    banddat = WannierIO.read_w90_band_dat(dataset"Si2/reference/Si2_band.dat")

    kpoint_path = generate_kpath(lattice, win.kpoint_path)
    kpi = Wannier.generate_w90_kpoint_path(kpoint_path)
    @test all(isapprox.(banddat.x, Wannier.get_linear_path(kpi); atol=1e-5))
end

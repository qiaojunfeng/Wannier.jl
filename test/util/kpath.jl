
@testset "get_kpath" begin
    test_data = YAML.load_file(String(@__DIR__) * "/test_data/kpath.yaml")

    win = read_win(joinpath(FIXTURE_PATH, "silicon/silicon.win"))

    # yaml_dict = Dict(
    #     "kpoint_path" => Dict(
    #         "basis" => win.kpoint_path.basis,
    #         "paths" => win.kpoint_path.paths,
    #         "points" => win.kpoint_path.points,
    #         "setting" => win.kpoint_path.setting,
    #     ),
    # )
    # YAML.write_file(String(@__DIR__) * "/test_data/kpath.yaml", yaml_dict)

    test_kpath = test_data["kpoint_path"]
    win_kpath = Wannier.get_kpath(win.unit_cell, win.kpoint_path)
    t_points = Dict(Symbol(k) => v for (k, v) in test_kpath["points"])
    @test t_points == win_kpath.points
    t_paths = [[Symbol(i) for i in l] for l in test_kpath["paths"]]
    @test t_paths == win_kpath.paths
    @test test_kpath["basis"] == win_kpath.basis
    @test Symbol(test_kpath["setting"]) == Symbol(win_kpath.setting)
end

@testset "interpolate w90 kpath" begin
    win = read_win(joinpath(FIXTURE_PATH, "valence/band/silicon.win"))
    lattice = win.unit_cell
    recip_lattice = Wannier.get_recip_lattice(lattice)
    kpi, E = read_w90_band(
        joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon"), recip_lattice
    )

    # num points of 1st segment
    n_points = 100
    kpath = Wannier.get_kpath(lattice, win.kpoint_path)
    test_kpi = Wannier.interpolate_w90(kpath, n_points)

    @test all(isapprox.(test_kpi.kpaths, kpi.kpaths; atol=1e-5))
    # If in the kpath block of win file, there are two kpoints with same label but
    # different coordinates, I will append a number to the repeated label in read_win,
    # so I only compare label without number.
    test_labels = test_kpi.labels
    for (i, lab) in enumerate(test_labels)
        new_lab = deepcopy(lab)
        for (k, v) in lab
            sv = String(v)
            delim = "_"
            if occursin(delim, sv)
                # only remove the last number suffix
                sv = join(split(sv, delim)[1:(end - 1)], delim)
                pop!(new_lab, k)
                push!(new_lab, k => Symbol(sv))
            end
        end
        test_labels[i] = new_lab
    end
    @test test_labels == kpi.labels
    @test test_kpi.basis ≈ kpi.basis
    @test Symbol(test_kpi.setting) == Symbol(kpi.setting)
end

@testset "get_x kpath" begin
    win = read_win(joinpath(FIXTURE_PATH, "valence/band/silicon.win"))
    lattice = win.unit_cell

    band = WannierIO.read_w90_band(joinpath(FIXTURE_PATH, "valence/band/mdrs/silicon"))

    kpoint_path = Wannier.get_kpath(lattice, win.kpoint_path)
    kpi = Wannier.interpolate_w90(kpoint_path, 100)

    @test all(isapprox.(band.x, Wannier.get_x(kpi); atol=1e-5))
end

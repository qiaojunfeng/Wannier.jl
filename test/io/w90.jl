using YAML
using Brillouin

mat2vec(A::AbstractMatrix) = [v for v in eachcol(A)]

@testset "read win" begin
    test_data = YAML.load_file(String(@__DIR__) * "/test_data/win.yaml")

    win = read_win(joinpath(FIXTURE_PATH, "silicon/silicon.win"))

    # Convert type so YAML can write it.
    unit_cell = mat2vec(win.unit_cell)
    atoms_frac = mat2vec(win.atoms_frac)
    kpoints = mat2vec(win.kpoints)

    # yaml_dict = Dict(
    #     "num_bands" => win.num_bands,
    #     "num_wann" => win.num_wann,
    #     "atoms_frac" => atoms_frac,
    #     "atom_labels" => win.atom_labels,
    #     "unit_cell" => unit_cell,
    #     "mp_grid" => win.mp_grid,
    #     "kpoints" => kpoints,
    #     "dis_froz_max" => win.dis_froz_max,
    #     "dis_win_max" => win.dis_win_max,
    #     "kpoint_path" => Dict(
    #         "basis" => win.kpoint_path.basis,
    #         "paths" => win.kpoint_path.paths,
    #         "points" => win.kpoint_path.points,
    #         "setting" => win.kpoint_path.setting,
    #     ),
    # )
    # YAML.write_file(String(@__DIR__) * "/test_data/win.yaml", yaml_dict)

    test_kpath = test_data["kpoint_path"]
    win_kpath = win.kpoint_path
    t_points = Dict(Symbol(k) => v for (k, v) in test_kpath["points"])
    @test t_points == win_kpath.points
    t_paths = [[Symbol(i) for i in l] for l in test_kpath["paths"]]
    @test t_paths == win_kpath.paths
    @test test_kpath["basis"] == win_kpath.basis
    @test Symbol(test_kpath["setting"]) == Symbol(win_kpath.setting)

    @test test_data["num_wann"] == win.num_wann
    @test test_data["num_bands"] == win.num_bands
    @test test_data["unit_cell"] ≈ unit_cell
    @test test_data["atoms_frac"] ≈ atoms_frac
    @test test_data["atom_labels"] == win.atom_labels
    @test test_data["mp_grid"] == win.mp_grid
    @test test_data["kpoints"] ≈ kpoints
    @test ismissing(win.dis_froz_min)
    @test ismissing(win.dis_win_min)
    @test test_data["dis_froz_max"] == win.dis_froz_max
    @test test_data["dis_win_max"] == win.dis_win_max
end

@testset "read/write mmn" begin
    M, kpb_k, kpb_b = read_mmn(joinpath(FIXTURE_PATH, "silicon/silicon.mmn"))

    tmpfile = tempname(; cleanup=true)

    write_mmn(tmpfile, M, kpb_k, kpb_b)

    M2, kpb_k2, kpb_b2 = read_mmn(tmpfile)

    @test M ≈ M2
    @test kpb_k ≈ kpb_k2
    @test kpb_b ≈ kpb_b2
end

@testset "read/write eig" begin
    E = read_eig(joinpath(FIXTURE_PATH, "silicon/silicon.eig"))

    tmpfile = tempname(; cleanup=true)
    write_eig(tmpfile, E)
    E2 = read_eig(tmpfile)

    @test E ≈ E2
end

@testset "read_w90" begin
    model = read_w90(joinpath(FIXTURE_PATH, "silicon/silicon"))

    @test model.n_bands ≈ 12
    @test model.n_wann ≈ 8
    @test model.n_kpts ≈ 64
end

@testset "read nnkp" begin
    test_data = YAML.load_file(String(@__DIR__) * "/test_data/nnkp.yaml")

    bvectors = read_nnkp(joinpath(FIXTURE_PATH, "silicon/silicon.nnkp"))

    # Convert type so YAML can write it.
    kpb_b = bvectors.kpb_b
    dict = Dict(
        "recip_lattice" => mat2vec(bvectors.recip_lattice),
        "kpoints" => mat2vec(bvectors.kpoints),
        "bvectors" => mat2vec(bvectors.bvectors),
        "kpb_k" => mat2vec(bvectors.kpb_k),
        "kpb_b" => [mat2vec(kpb_b[:, :, ik]) for ik in axes(kpb_b, 3)],
    )

    # YAML.write_file(String(@__DIR__) * "/test_data/nnkp.yaml", dict)

    for (key, value) in dict
        @test value ≈ test_data[key]
    end
end

@testset "read/write nnkp" begin
    bvectors = read_nnkp(joinpath(FIXTURE_PATH, "silicon/silicon.nnkp"))
    tmpfile = tempname(; cleanup=true)
    n_wann = 8
    write_nnkp(tmpfile, bvectors, n_wann)

    bvectors2 = read_nnkp(tmpfile)
    @test bvectors.recip_lattice ≈ bvectors2.recip_lattice
    @test bvectors.kpoints ≈ bvectors2.kpoints
    @test bvectors.bvectors ≈ bvectors2.bvectors
    @test bvectors.kpb_k ≈ bvectors2.kpb_k
    @test bvectors.kpb_b ≈ bvectors2.kpb_b
end

@testset "read/write unk" begin
    ik, Ψ = read_unk(joinpath(FIXTURE_PATH, "silicon/UNK00001.1"))

    tmpfile = tempname(; cleanup=true)
    write_unk(tmpfile, ik, Ψ)
    ik2, Ψ2 = read_unk(tmpfile)

    @test ik ≈ ik2
    @test Ψ ≈ Ψ2
end

@testset "read chk" begin
    chk = read_chk(joinpath(FIXTURE_PATH, "silicon/silicon.chk.fmt"))

    @test chk.n_wann == 8
    @test chk.n_bands == 12
end

@testset "read/write w90 band dat" begin
    band = read_w90_band(joinpath(FIXTURE_PATH, "valence/band/silicon"))

    outdir = mktempdir(; cleanup=true)
    outseedname = joinpath(outdir, "silicon")

    write_w90_band(
        outseedname, band.kpoints, band.E, band.x, band.symm_idx, band.symm_label
    )

    band2 = read_w90_band(outseedname)

    @test band.kpoints ≈ band2.kpoints
    @test band.E ≈ band2.E
    @test band.x ≈ band2.x
    @test band.symm_idx == band2.symm_idx
    @test band.symm_label == band2.symm_label
end

@testset "read/write w90 band dat kpi" begin
    win = read_win(joinpath(FIXTURE_PATH, "valence/band/silicon.win"))
    recip_lattice = Wannier.get_recip_lattice(win.unit_cell)
    kpi, E = read_w90_band(joinpath(FIXTURE_PATH, "valence/band/silicon"), recip_lattice)

    outdir = mktempdir(; cleanup=false)#true)
    outseedname = joinpath(outdir, "silicon")

    write_w90_band(outseedname, kpi, E)

    kpi2, E2 = read_w90_band(outseedname, recip_lattice)

    @test kpi.kpaths ≈ kpi2.kpaths
    @test kpi.labels == kpi2.labels
    @test kpi.basis ≈ kpi2.basis
    @test Symbol(kpi.setting) == Symbol(kpi2.setting)
    @test E ≈ E2
end

@testset "read/write chk" begin
    chk = read_chk(joinpath(FIXTURE_PATH, "silicon/silicon.chk.fmt"))

    win = read_win(joinpath(FIXTURE_PATH, "silicon/silicon.win"))
    @test chk.lattice ≈ win.unit_cell

    tmpfile = tempname(; cleanup=true)
    write_chk(tmpfile, chk)
    chk2 = read_chk(tmpfile)

    @test chk.header == chk2.header
    @test chk.n_bands == chk2.n_bands
    @test chk.n_exclude_bands == chk2.n_exclude_bands
    @test chk.exclude_bands == chk2.exclude_bands
    @test chk.lattice ≈ chk2.lattice
    @test chk.recip_lattice ≈ chk2.recip_lattice
    @test chk.n_kpts ≈ chk2.n_kpts
    @test chk.kgrid == chk2.kgrid
    @test chk.kpoints ≈ chk2.kpoints
    @test chk.n_bvecs == chk2.n_bvecs
    @test chk.n_wann == chk2.n_wann
    @test chk.checkpoint == chk2.checkpoint
    @test chk.have_disentangled == chk2.have_disentangled
    @test chk.ΩI ≈ chk2.ΩI
    @test chk.dis_bands == chk2.dis_bands
    @test chk.Uᵈ ≈ chk2.Uᵈ
    @test chk.U ≈ chk2.U
    @test chk.M ≈ chk2.M
    @test chk.r ≈ chk2.r
    @test chk.ω ≈ chk2.ω
    @test chk.n_dis == chk2.n_dis
end

@testset "read wout" begin
    wout = read_wout(joinpath(FIXTURE_PATH, "valence/band/silicon.wout"))

    ref_lattice =
        [
            -2.698804 0.000000 2.698804
            0.000000 2.698804 2.698804
            -2.698804 2.698804 0.000000
        ]'
    ref_atom_labels = ["Si", "Si"]
    ref_atom_positions = [
        -0.25000 0.75000 -0.25000
        0.00000 0.00000 0.00000
    ]'
    ref_centers =
        [
            -0.659352 0.658238 -0.680969
            0.669283 0.695828 0.666806
            0.682490 -0.683846 -0.683726
            -0.701673 -0.656575 0.703751
        ]'
    ref_spreads = [
        2.39492617
        2.19372718
        1.83863803
        1.88512458
    ]

    @test wout.lattice ≈ ref_lattice
    @test wout.atom_labels == ref_atom_labels
    @test wout.atom_positions ≈ ref_atom_positions
    @test wout.centers ≈ ref_centers
    @test wout.spreads ≈ ref_spreads
end


@testset "get_kpoint_mappings" begin
    win = read_win(joinpath(FIXTURE_PATH, "silicon/silicon.win"))
    kpoints = win.kpoints
    kgrid = win.mp_grid

    k_xyz, xyz_k = Wannier.get_kpoint_mappings(kpoints, kgrid)

    n_kx, n_ky, n_kz = kgrid
    n_kpts = prod(kgrid)

    k_xyz_ref = Vector{Wannier.Vec3{Int}}(undef, n_kpts)
    xyz_k_ref = Array{Int,3}(undef, n_kx, n_ky, n_kz)

    # z increases the fastest
    for ikx in 1:n_kx
        for iky in 1:n_ky
            for ikz in 1:n_kz
                ik = (ikx - 1) * n_kx * n_ky + (iky - 1) * n_ky + ikz
                k_xyz_ref[ik] = [ikx, iky, ikz]
                xyz_k_ref[k_xyz_ref[ik]...] = ik
            end
        end
    end

    @test k_xyz == k_xyz_ref
    @test xyz_k == xyz_k_ref
end

@testset "get_kpoints" begin
    kpoints = Wannier.get_kpoints([2, 2, 2])

    ref_kpoints = Vec3[
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5],
        [0.0, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.0],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.5],
    ]
    @test kpoints == ref_kpoints
end

@testset "get_kpoints endpoint" begin
    kpoints = Wannier.get_kpoints([2, 2, 2]; endpoint=true)

    ref_kpoints = Vec3[
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0],
        [0.0, 1.0, 1.0],
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 0.0],
        [1.0, 1.0, 1.0],
    ]
    @test kpoints == ref_kpoints
end

@testset "sort_kpoints" begin
    ref_kpoints = Vec3[
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.5],
        [0.0, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.0],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.5],
    ]
    unsorted_kpoints = deepcopy(ref_kpoints)
    unsorted_kpoints[1] = Vec3(0.5, 0.5, 0.0)
    unsorted_kpoints[end - 1] = Vec3(0, 0, 0)
    sorted_kpoints = Wannier.sort_kpoints(unsorted_kpoints)
    @test sorted_kpoints ≈ ref_kpoints
end

@testset "get_kgrid" begin
    kgrid = [2, 2, 2]
    kpoints = Wannier.get_kpoints(kgrid)
    kpoints[1] = Vec3(0.5, 0.5, 0.0)
    kpoints[end - 1] = Vec3(0, 0, 0)
    @test Wannier.get_kgrid(kpoints) == kgrid
end

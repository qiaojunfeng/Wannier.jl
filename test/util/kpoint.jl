
@testset "get_kpoint_mappings" begin
    win = read_win(joinpath(FIXTURE_PATH, "silicon/silicon.win"))
    kpoints = win.kpoints
    kgrid = win.mp_grid

    k_xyz, xyz_k = get_kpoint_mappings(kpoints, kgrid)

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

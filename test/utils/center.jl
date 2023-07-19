@testitem "find_nearest_neighbor" begin
    using Wannier: Vec3

    lattice = [
        -2.698804 0.000000 -2.698804
        0.000000 2.698804 2.698804
        2.698804 2.698804 0.000000
    ]
    # atom_positions in fractional
    atom_positions = [[-0.25000, 0.75000, -0.25000], [0.00000, 0.00000, 0.00000]]
    # centers in Cartesian
    centers = [
        [-0.659352, 0.658238, -0.680969],
        [0.669283, 0.695828, 0.666806],
        [0.682490, -0.683846, -0.683726],
        [-0.701673, -0.656575, 0.703751],
    ]
    centers = map(c -> inv(lattice) * c, centers)  # to fractional
    distances, indices, translations = find_nearest_neighbor(
        centers, lattice, atom_positions
    )

    ref_distances = [1.15401, 1.16433, 1.15363, 1.14736]
    ref_indices = [2, 1, 1, 1]
    ref_translations = Vec3[[0, 0, 0], [0, 0, 0], [0, -1, 0], [1, -1, 0]]

    @test all(isapprox.(distances, ref_distances; atol=2e-5))
    @test indices == ref_indices
    @test translations == ref_translations
end

@testitem "wrap_centers" begin
    using Wannier: Vec3

    lattice = [
        -2.698804 0.000000 -2.698804
        0.000000 2.698804 2.698804
        2.698804 2.698804 0.000000
    ]
    # centers in Cartesian
    centers = [
        [-0.659352, 0.658238, -0.680969],
        [0.669283, 0.695828, 0.666806],
        [0.682490, -0.683846, -0.683726],
        [-0.701673, -0.656575, 0.703751],
    ]
    centers = Wannier.wrap_centers(centers, lattice)

    ref_centers = Vec3[
        [-3.35816, 3.35704, 4.71664],
        [-4.72832, 3.39463, 3.36561],
        [-4.71512, 4.71376, 4.71388],
        [-3.40048, 4.74103, 3.40255],
    ]
    @test all(isapprox.(centers, ref_centers; atol=2e-5))
end

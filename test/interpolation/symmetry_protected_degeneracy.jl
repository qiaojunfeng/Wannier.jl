@testitem "symmetry-protected interpolation degeneracies" begin
    using WannierIO
    using LinearAlgebra

    function synthetic_model(lattice, grid_size)
        number_1, number_2, number_3 = grid_size
        kpoints = [
            Vec3(i / number_1, j / number_2, k / number_3)
                for k in 0:(number_3 - 1)
                for j in 0:(number_2 - 1)
                for i in 0:(number_1 - 1)
        ]
        stencil = generate_kspace_stencil(
            reciprocal_lattice(lattice),
            Vec3(grid_size),
            kpoints,
            Wannier.CubicNearestKSpaceStencil(),
        )
        number_bands = 2
        number_kpoints = length(kpoints)
        overlaps = zeros(
            ComplexF64,
            number_bands,
            number_bands,
            Wannier.n_bvectors(stencil),
            number_kpoints,
        )
        gauges = zeros(ComplexF64, number_bands, number_bands, number_kpoints)
        identity_2 = Matrix{ComplexF64}(I, number_bands, number_bands)
        for kpoint_index in 1:number_kpoints
            gauges[:, :, kpoint_index] .= identity_2
            for neighbor_index in 1:Wannier.n_bvectors(stencil)
                overlaps[:, :, neighbor_index, kpoint_index] .= identity_2
            end
        end

        eigenvalues = Matrix{Float64}(undef, number_bands, number_kpoints)
        for (kpoint_index, kpoint) in enumerate(kpoints)
            kpoint_1, kpoint_2, kpoint_3 = kpoint
            eigenvalues[1, kpoint_index] =
                cospi(2kpoint_1) +
                0.37sinpi(2kpoint_2) +
                0.11cospi(2kpoint_3)
            eigenvalues[2, kpoint_index] =
                -0.43cospi(2kpoint_1) +
                0.29sinpi(2(kpoint_1 + kpoint_2)) -
                0.07cospi(2kpoint_3)
        end
        return Wannier.Model(
            lattice,
            Vec3{Float64}[],
            String[],
            stencil,
            overlaps,
            gauges,
            eigenvalues,
            falses(number_bands, number_kpoints),
        )
    end

    function make_symmetry(
            real_space_matrices,
            representations;
            antiunitary = falses(length(real_space_matrices)),
            inverses,
        )
        identity_spin = Matrix{ComplexF64}(I, 2, 2)
        operations = map(eachindex(real_space_matrices)) do index
            real_space_matrix = Wannier.Mat3{Int}(real_space_matrices[index])
            reciprocal_matrix = Wannier.Mat3{Int}(transpose(inv(real_space_matrix)))
            WannierIO.SymOp(
                "synthetic symmetry $index",
                real_space_matrix,
                Vec3(0.0, 0.0, 0.0),
                reciprocal_matrix,
                antiunitary[index],
                identity_spin,
                index,
                inverses[index],
            )
        end
        orbital_representations = map(eachindex(representations)) do index
            WannierIO.OrbitalRep{2}(index, representations[index])
        end
        centers = [Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0)]
        return WannierSymmetry(operations, orbital_representations, centers)
    end

    identity_2 = [1 0; 0 1]
    threefold = [-1 -1; 1 0]
    mirror = [1 1; 0 -1]
    planar_operations = [
        identity_2,
        threefold,
        threefold^2,
        mirror,
        mirror * threefold,
        mirror * threefold^2,
    ]
    spatial_operations = map(planar_operations) do operation
        [
            operation[1, 1] operation[1, 2] 0
            operation[2, 1] operation[2, 2] 0
            0 0 1
        ]
    end
    angle = 2π / 3
    rotation = [cos(angle) -sin(angle); sin(angle) cos(angle)]
    reflection = [1.0 0.0; 0.0 -1.0]
    representations = [
        Matrix{Float64}(I, 2, 2),
        rotation,
        rotation^2,
        reflection,
        reflection * rotation,
        reflection * rotation^2,
    ]
    d3 = make_symmetry(
        spatial_operations, representations; inverses = [1, 3, 2, 4, 5, 6]
    )

    hexagonal_lattice = Wannier.Mat3{Float64}(
        [1.0 0.5 0.0; 0.0 sqrt(3) / 2 0.0; 0.0 0.0 3.0]
    )
    graphene_like_model = synthetic_model(hexagonal_lattice, (3, 3, 1))
    graphene_like_interpolation = InterpolationModel(
        graphene_like_model; real_space = WignerSeitz(), symmetry = d3
    )
    k_corner = Vec3(2 / 3, 1 / 3, 0.0)
    @test all(
        Wannier.isequivalent(Wannier.rotate_kpoint(k_corner, operation), k_corner)
            for operation in d3.symops
    )
    graphene_like_energies = interpolate(
        graphene_like_interpolation,
        [k_corner, k_corner + Vec3(1.0e-3, -2.0e-3, 0.0)],
        BandEnergy(),
    ).band_energy
    @test abs(graphene_like_energies[2, 1] - graphene_like_energies[1, 1]) < 1.0e-13
    @test abs(graphene_like_energies[2, 2] - graphene_like_energies[1, 2]) > 1.0e-3
    @test Wannier._symmetry_covariance_residual(
        graphene_like_interpolation, :hamiltonian, [Vec3(0.123, 0.234, 0.0)]
    ).maximum < 1.0e-13

    identity_3 = Matrix{Int}(I, 3, 3)
    spinful_time_reversal = [0.0 1.0; -1.0 0.0]
    @test spinful_time_reversal * conj(spinful_time_reversal) ≈ -I
    kramers_symmetry = make_symmetry(
        [identity_3, identity_3],
        [Matrix{Float64}(I, 2, 2), spinful_time_reversal];
        antiunitary = BitVector([false, true]),
        inverses = [1, 2],
    )
    spinful_model = synthetic_model(Wannier.Mat3{Float64}(I), (3, 3, 3))
    spinful_interpolation = InterpolationModel(
        spinful_model; real_space = WignerSeitz(), symmetry = kramers_symmetry
    )
    time_reversal_invariant_momenta = [
        Vec3(i / 2, j / 2, k / 2) for k in 0:1 for j in 0:1 for i in 0:1
    ]
    kramers_energies = interpolate(
        spinful_interpolation, time_reversal_invariant_momenta, BandEnergy()
    ).band_energy
    @test maximum(abs, diff(kramers_energies; dims = 1)) < 2.0e-13

    generic_energies = interpolate(
        spinful_interpolation, [Vec3(0.13, 0.21, 0.17)], BandEnergy()
    ).band_energy
    @test abs(generic_energies[2, 1] - generic_energies[1, 1]) > 0.05
    @test Wannier._symmetry_covariance_residual(
        spinful_interpolation, :hamiltonian, [Vec3(0.13, 0.21, 0.17)]
    ).maximum < 1.0e-13
end

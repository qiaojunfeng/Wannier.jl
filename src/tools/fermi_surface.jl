"""
Interpolate Fermi surface.

# Arguments

- `model`: interpolation model

# Options

- `--nk`: number of interpolated kpoints along each reciprocal lattice vector.
    Default is `50`. Note in output bxsf, there will be 51 points, since the
    bxsf format requires general grid, i.e., the last kpoint is the periodic
    image of the first one. This also restores the behavior of wannier90.
- `--dk`: spacing of interpolated kpoints along each reciprocal lattice vector,
    in unit of `Å^-1`. Default is a negative number, meaning I will use `--nk`
    instead; if `dk` > 0, `dk` will take precedence.
    The number of interpolated kpoints is the nearest integer to `b/dk`,
    where `b` is the norm of each reciprocal lattice vector. Also plus 1 due to
    the same reason in `--nk`.
- `--ef`: Fermi energy. Default is `0.0`.
- `--outprefix`: output prefix for bxsf. Default is `wjl`
"""
function fermisurf(
        model::InterpolationModel;
        nk::Int = 50,
        dk::Float64 = -0.04,
        ef::Float64 = 0.0,
        outprefix::String = "wjl",
    )
    @printf("* Fermi surface interpolation started on %s\n", string(now()))

    recip_latt = reciprocal_lattice(model)

    if dk > 0
        nk = map(eachcol(recip_latt)) do b
            round(Int, norm(b) / dk)
        end
    end

    nks = _get_bxsf_nks(nk)
    @printf("Interpolation grid: %d %d %d\n", nks...)

    # kpoints are in fractional coordinates
    kpoints = Wannier.get_kpoints(nks; endpoint = true)

    eigenvalues = interpolate(model, kpoints, BandEnergy()).band_energy

    # kpoints = reshape(kpoints, 3, n_kz, n_ky, n_kx)
    # E = reshape(E, nwann, n_kz, n_ky, n_kx)
    # # and permutedims
    # kpoints = permutedims(kpoints, (1, 4, 3, 2))
    # E = permutedims(E, (1, 4, 3, 2))

    # reshape eigenvalues to a nwann * nkx * nky * nkz array
    E = zeros(n_wannier(model), nks...)
    # The kz increase the fastest in kpoints
    counter = 1
    for i in 1:nks[1]
        for j in 1:nks[2]
            for k in 1:nks[3]
                E[:, i, j, k] .= eigenvalues[:, counter]
                counter += 1
            end
        end
    end

    coordinates = map(size -> collect(range(0.0, 1.0; length = size)), nks)
    bxsf = WannierIO.Bxsf(
        ef,
        Vec3(0.0, 0.0, 0.0),
        recip_latt,
        coordinates...,
        E,
    )
    WannierIO.write_bxsf("$outprefix.bxsf", bxsf)

    @printf("* Fermi surface interpolation finished on %s\n", string(now()))
    return nothing
end

function _get_bxsf_nks(nks::AbstractVector)
    @assert length(nks) == 3 "nks must be a 3-element vector"
    # bxsf need general grid, i.e., the last kpoint is the periodic image of the
    # first one, so I increase the n_k by 1
    return map(n -> n + 1, nks)
end

function _get_bxsf_nks(nk::Integer)
    return _get_bxsf_nks([nk, nk, nk])
end

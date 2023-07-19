using Printf: @printf
using Dates: now

"""
Interpolate Fermi surface.

# Args

- `seedname`: seedname for `amn`/`mmn`/`eig`/`chk` files or `tb/wsvec.dat` file

# Options

- `--nk`: number of interpolated kpoints along each reciprocal lattice vector. Default is `50`,
    on output bxsf, there will be 51 points (the last point is the periodic image of the first point).
- `--dk`: spacing of interpolated kpoints along each reciprocal lattice vector, in unit of `Å^-1`.
    Default is a negative number, meaning I will use `--nk` instead; if dk > 0, dk will take precedence.
    The number of interpolated kpoints is the nearest integer to `a/dk`,
    where `a` is the norm of each reciprocal lattice vector. Also plus 1 due to
    the same reason in `--nk`.
- `--ef`: Fermi energy. Default is `0.0`.
- `--out`: output seedname for bxsf. Default is `wjl`
- `--amn`: `amn` filename. If not given, default is read `chk.fmt` file

# Flags

- `--tb`: read `tb/wsvec.dat` file instead of `amn`/`mmn`/`eig`/`chk` files. Default is `false`.
"""
@cast function fermisurf(
    seedname::String;
    nk::Int=50,
    dk::Float64=-0.04,
    ef::Float64=0.0,
    out::String="wjl",
    amn::String="",
    tb::Bool=false,
)
    @printf("* Fermi surface interpolation started on %s\n", string(now()))

    if tb
        model = read_w90_tb(seedname)
    else
        if isempty(amn)
            model = read_w90_with_chk(seedname)
        else
            model = read_w90_with_chk(seedname; chk=false, amn=amn)
        end
    end
    Rvecs = model.kRvectors.Rvectors
    _print_type(Rvecs)

    recip_latt = reciprocal_lattice(Rvecs.lattice)

    if dk > 0
        nk = [round(Int, norm(b) / dk) for b in eachcol(recip_latt)]
    end
    kpoints, E = fermi_surface(Rvecs, model.H; n_k=nk)

    origin = zeros(Float64, 3)
    WannierIO.write_bxsf("$out.bxsf", ef, origin, recip_latt, E)

    @printf("* Fermi surface interpolation finished on %s\n", string(now()))
    return nothing
end

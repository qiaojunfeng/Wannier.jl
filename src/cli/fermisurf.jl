"""
Interpolate Fermi surface.

# Args

- `seedname`: seedname for `amn`/`mmn`/`eig`/`chk` files or `tb/wsvec.dat` file

# Options

- `--nk`: number of interpolated kpoints along each reciprocal lattice vector. Default is `50`,
    on output bxsf, there will be 51 points (the last point is the periodic image of the first point).
- `--ef`: Fermi energy. Default is `0.0`.
- `--out`: output seedname for bxsf. Default is `wjl`
- `--amn`: `amn` filename. If not given, default is read `chk.fmt` file

# Flags

- `--tb`: read `tb/wsvec.dat` file instead of `amn`/`mmn`/`eig`/`chk` files. Default is `false`.
- `--ws`: force using WS interpolation to reproduce Wannier90 results. Default is `false`.
"""
@cast function fermisurf(
    seedname::String;
    nk::Int=50,
    ef::Float64=0.0,
    out::String="wjl",
    amn::String="",
    tb::Bool=false,
    ws::Bool=false,
)
    if tb
        Rvecs, H, r = read_w90_tb(seedname)
    else
        if isempty(amn)
            model = read_w90_post(seedname)
        else
            model = read_w90_post(seedname; chk=false, amn=amn)
        end
        Rvecs = model.kRvectors.Rvectors
    end
    _print_type(Rvecs)

    # If you really want, you can use WS interpolation.
    # Note W90 only use WS interpolation when computing the Fermi surface,
    # however we use MDRS here (if the tb.dat is in MDRS), so the results are (usually) different.
    if ws && (Rvecs isa RVectorsMDRS)
        Rvecs = Rvecs.Rvectors
    end
    kpoints, E = fermi_surface(Rvecs, H; n_k=nk)

    origin = zeros(Float64, 3)
    recip_latt = get_recip_lattice(Rvecs.lattice)
    write_bxsf("$out.bxsf", ef, origin, recip_latt, E)

    return nothing
end

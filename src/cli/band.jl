"""
Interpolate band structure.

# Args

- `seedname`: seedname for `amn`/`mmn`/`eig`/`chk` files

# Options

- `--out`: output seedname for bands.dat. Default is `wjl`
- `--amn`: `amn` filename. If not given, default is read `chk.fmt` file

"""
@cast function band(seedname::String; out::String="wjl", amn::String="")
    if isempty(amn)
        model = read_w90_post(seedname)
    else
        model = read_w90_post(seedname; chk=false, amn=amn)
    end

    type_Rvectors = typeof(model.kRvectors.Rvectors)
    if type_Rvectors <: RVectorsMDRS
        interp_type = "MDRS"
    elseif type_Rvectors <: RVectors
        interp_type = "WS"
    else
        interp_type = "???"
    end
    @info "Using $interp_type interpolation"
    println()

    kpi, E = interpolate(model)

    write_w90_band(out, kpi, E)
    return nothing
end

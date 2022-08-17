"""
Interpolate band structure.

# Args

- `seedname`: seedname for AMN/MMN/EIG/CHK files

# Options

- `--out`: output seedname for bands.dat. Default is `wjl`
- `--amn`: AMN filename. If not given, default is read chk.fmt file

"""
@cast function band(seedname::String; out::String="wjl", amn::String="")
    if isempty(amn)
        model = read_w90_post(seedname)
    else
        model = read_w90_post(seedname; chk=false, amn=amn)
    end

    kpi, E = interpolate(model)

    write_w90_band(out, kpi, E)
    return nothing
end

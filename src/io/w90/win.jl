"""Check `win.use_ws_distance`."""
ismdrs(win::NamedTuple) = get(win, :use_ws_distance, true)

"""
Construct a `KPath` from `win` input `unit_cell_cart` and `kpoint_path`.

If `kpoint_path` is empty, an empty `KPath` will be returned.
"""
function get_kpath(win::NamedTuple)
    if haskey(win, :kpoint_path)
        kpath = get_kpath(win.unit_cell_cart, win.kpoint_path)
    else
        # an empty kpath
        points = Dict{Symbol,Vec3{Float64}}()
        paths = Vector{Vector{Symbol}}()
        basis = ReciprocalBasis([v for v in eachcol(model.recip_lattice)])
        setting = Ref(Brillouin.LATTICE)
        kpath = KPath{3}(points, paths, basis, setting)
    end
    return kpath
end

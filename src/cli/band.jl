"""
Interpolate band structure.

Read tight-binding Hamiltonian from `tb.dat`/`wsvec.dat` files, kpoint path
from `prefix.win` file, and interpolate band structure.

# Args

- `prefix`: prefix of `prefix.win`/`prefix_tb.dat`/`prefix_wsvec.dat` files

# Options

- `--outprefix`: prefix for output `_band.dat`/`_band.kpt`/`_band.labelinfo.dat`
    files. Default is `wjl`.

"""
@cast function band(prefix::String; outprefix::String="wjl")
    hamiltonian, _ = read_w90_tb(prefix)
    win = read_win("$prefix.win")
    kpi = generate_w90_kpoint_path(win.unit_cell_cart, win.kpoint_path)
    interp = HamiltonianInterpolator(hamiltonian)
    eigenvalues, _ = interp(kpi)
    write_w90_band(outprefix, kpi, eigenvalues)
    return nothing
end

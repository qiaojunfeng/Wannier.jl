# Must prepend a dot due to Requires.jl
using .PlotlyJS

export plot_band_diff

"""
    plot_band_diff(kpi::KPathInterpolant, E1::AbstractArray, E22::AbstractArray; kwargs...)

Compare two band structures.

`E1` in black, while `E2` in dashed orange.

# Arguments
- `kpi`: KPathInterpolant
- `E1`: band energies, 1D of length `n_kpts`, or 2D array of size `n_bands * n_kpts`
- `E2`: band energies, 1D of length `n_kpts`, or 2D array of size `n_bands * n_kpts`

# Keyword Arguments
See also the keyword arguments of [`_get_band_plot`](@ref).
"""
function plot_band_diff(
    kpi::KPathInterpolant, E1::AbstractArray, E2::AbstractArray; kwargs...
)
    x = get_x(kpi)
    symm_idx, symm_label = _get_symm_idx_label(kpi)

    P1 = _get_band_plot(
        x, E1; color="black", symm_idx=symm_idx, symm_label=symm_label, kwargs...
    )
    # orange and slightly thinner
    P2 = _get_band_plot(
        x,
        E2;
        color="orange",
        dash="dash",
        width=0.9,
        symm_idx=symm_idx,
        symm_label=symm_label,
        kwargs...,
    )

    addtraces!(P1, P2.data...)

    return plot(P1)
end

"""
Plot density of states.
"""
function Wannier.plot_dos(x::AbstractVector, y::AbstractVector; kwargs...)
    P = Wannier.get_dos_plot(x, y; kwargs...)
    return PlotlyJS.plot(P)
end

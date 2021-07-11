#!/usr/bin/env julia
import ArgParse
import Wannier as Wan

function parse_commandline()
    s = ArgParse.ArgParseSettings()
    @ArgParse.add_arg_table s begin
        "wout"
            help = "Filename of Wannier90 win output file"
            required = true
    end
    return ArgParse.parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    
    f_wout = parsed_args["wout"]

    wout = Wan.InOut.read_wout(f_wout)

    num_wann = length(wout["wf_spreads"])
    dists = zeros(Float64, num_wann)
    for iw = 1:num_wann
        d, _, _ = Wan.find_nearests(wout["unit_cell"], wout["atoms"], wout["wf_centers"][:,iw];
        reduced_coord=false, search_neighbors=1)
        dists[iw] = d[1]
    end

    @info "unit_cell" wout["unit_cell"]
    @info "atoms" wout["atoms"]
    @info "wf_spreads" wout["wf_spreads"]'
    @info "wf_centers" wout["wf_centers"]
    @info "distance to nearest atom" dists'
    @info "sum(dists), average(dists)" sum(dists) sum(dists)/num_wann
end

main()

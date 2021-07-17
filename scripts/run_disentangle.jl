#!/usr/bin/env julia
import Random
import LinearAlgebra as LA
import ArgParse
import TOML
import Configurations
import Wannier as Wan
import Serialization

function parse_commandline()
    s = ArgParse.ArgParseSettings()
    @ArgParse.add_arg_table s begin
        # "--opt1"
        #     help = "an option with an argument"
        # "--opt2", "-o"
        #     help = "another option with an argument"
        #     arg_type = Int
        #     default = 0
        # "--flag1"
        #     help = "an option without argument, i.e. a flag"
        #     action = :store_true
        "seedname"
            help = "Name of the input TOML file: seedname.toml"
            required = true
    end
    return ArgParse.parse_args(s)
end

function check_inputs(params::Wan.Parameters.InputParams)
    # check parameters
    if params.do_randomize_gauge && params.read_amn
        error("do not set do_randomize_gauge and read_amn")
    end
end

"""
check consistancy between input parameters and readed amn, mmn matrices
"""
function check_inputs2(data::Wan.Parameters.CoreData, params::Wan.Parameters.InputParams)
    if params.dis_froz_pao
        if params.dis_froz_num != data.num_wann
            error("if using dis_froz_pao, dis_froz_num should equal to number of WFs")
        end
    end
end

function wannierize(params::Wan.Parameters.InputParams)
    Random.seed!(0)

    data = Wan.InOut.read_seedname(
        params.seed_name, amn=params.read_amn, eig=params.read_eig)

    # check_inputs2(data, params)

    if params.read_amn
        amn0 = copy(data.amn)
    else
        amn0 = randn(size(data.amn)) + im * randn(size(data.amn))
    end

    # Orthonormalize Amn, make it semiunitary
    for ik = 1:data.num_kpts
        amn0[:,:,ik] = Wan.Utilities.orthonormalize_lowdin(amn0[:,:,ik])
    end

    Wan.Disentangle.set_frozen_bands!(data, params)

    if params.dis_froz_pao
        V = zeros(ComplexF64, data.num_bands, data.num_bands, data.num_kpts)
        for ik = 1:data.num_kpts
            # amn0[:,:,ik], V[:,:,ik] = Wan.Disentangle.max_projectability(amn0[:,:,ik])
            amn0[:,:,ik], V[:,:,ik] = Wan.Disentangle.maxproj_froz(amn0[:,:,ik], data.frozen[:,ik])
        end
        for ik = 1:data.num_kpts
            for ib = 1:data.num_bvecs
                data.mmn[:,:,ib,ik] = V[:,:,ik]' * data.mmn[:,:,ib,ik] * V[:,:,data.kpbs[ib,ik]]
                # for band interpolation we need to rotate data.eig as well!
            end
        end
    end
    
    for ik = 1:data.num_kpts
        l_frozen = data.frozen[:,ik]
        l_non_frozen = .!l_frozen
        # initial guess for U matrix
        amn0[:,:,ik] = Wan.Disentangle.orthonormalize_and_freeze(amn0[:,:,ik], l_frozen, l_non_frozen)
    end

    # for ik = 1:data.num_kpts
    #     for ib = 1:data.num_bvecs
    #         A1 = amn0[:,:,ik]
    #         A2 = amn0[:,:,data.kpbs[ib,ik]]
    #         @info "projectability @ $ik" real(LA.diag(proj)') real(LA.diag(A_new * A_new')')
    #         amn0[:,:,ik] = A_new
    #     end
    # end

    # 
    A = Wan.Disentangle.minimize(data, params, amn0)

    # fix global phase
    if params.do_normalize_phase
        for i = 1:data.num_wann
            imax = indmax(abs.(A[:,i,1,1,1]))
            @assert abs(A[imax,i,1,1,1]) > 1e-2
            A[:,i,:,:,:] *= conj(A[imax,i,1,1,1] / abs(A[imax,i,1,1,1]))
        end
    end

    if params.dis_froz_pao
        for ik = 1:data.num_kpts
            A[:,:,ik] = V[:,:,ik] * A[:,:,ik]
        end
    end

    if params.write_optimized_amn
        Wan.InOut.write_amn("$(params.seed_name).amn.optimized", A)
    end
end

function parse_inputs(seed_name::String)
    toml_suffix = ".toml"
    if endswith(seed_name, toml_suffix)
        seed_name = seed_name[1:length(seed_name)-length(toml_suffix)]
    end

    toml = TOML.parsefile("$(seed_name)$(toml_suffix)")
    toml["seed_name"] = seed_name
    params = Configurations.from_dict(Wan.Parameters.InputParams, toml)

    # non-core stuffs in win are also stored in params
    win = Wan.InOut.read_win("$(params.seed_name).win")
    params.kpath = win["kpath"]
    params.kpath_label = win["kpath_label"]

    check_inputs(params)

    # println(params)
    # exit()

    return params
end

function main()
    parsed_args = parse_commandline()

    params = parse_inputs(parsed_args["seedname"])

    if params.restart == "plot"
        Wan.interpolate(params)
    else
        wannierize(params)
    end
end

main()

# Notes on band interpolation
# After seedname.amn.optimized is generated, when using w90 to 
# interpolate bands, w90 does
#    1. Lowdin orthogonalization of AMN matrix, in disentangle.F90:dis_project, line 1418
#       This should do no harm, since the optmized amn is already semi-unitary, 
#       a SVD of it should not change the optmized amn (apart from numerical noise)
#    2. Generate a new amn according to frozen window, in disentangle.F90:dis_proj_froz, line 1830
#       This will DESTROY the optmized amn matrix, if we restart w90 from the optmized Amn
#       with dis_num_iter = 0, the wout spreads is very different from the output of wannier.jl,
#       we must skip this step by commenting out ALL the dis_froz_min/max in the win file and
#       use w90 to interpolate bands, remember also set num_iter and dis_num_iter = 0.
# In the future I might use FFT to directly interpolate bands in wannier.jl

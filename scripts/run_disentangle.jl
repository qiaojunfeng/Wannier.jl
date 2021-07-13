import Random
import LinearAlgebra as LA
import ArgParse
import TOML
import Configurations
import Wannier as Wan

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

function check_inputs(inputs::Wan.Parameters.InputParams)
    # check parameters
    if inputs.do_randomize_gauge && inputs.read_amn
        error("do not set do_randomize_gauge and read_amn")
    end
end

function wannierize(params::Wan.Parameters.InputParams)
    Random.seed!(0)

    data = Wan.InOut.read_seedname(
        params.seed_name, params.read_amn, params.read_eig)

    if params.read_amn
        amn0 = copy(data.amn)
    else
        amn0 = randn(size(data.amn)) + im * randn(size(data.amn))
    end

    # Orthonormalize Amn, make it semiunitary
    for ik = 1:data.num_kpts
        amn0[:,:,ik] = Wan.Utilities.orthonormalize_lowdin(amn0[:,:,ik])
    end

    Wan.Disentangle.set_frozen_bands(data, params)

    # initial guess for U matrix
    # Vt = zeros(ComplexF64, data.num_bands, data.num_bands, data.num_kpts)
        # amn0[:,:,ik] = Wan.Disentangle.max_projectability(amn0[:,:,ik])
    
    for ik = 1:data.num_kpts
        l_frozen = data.frozen[:,ik]
        l_non_frozen = .!l_frozen
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

    if params.write_optimized_amn
        Wan.InOut.write_amn("$(params.seed_name).amn.optimized", A)
    end
end

function main()
    parsed_args = parse_commandline()

    seedname = parsed_args["seedname"]
    toml_suffix = ".toml"
    if endswith(seedname, toml_suffix)
        seedname = seedname[1:length(seedname)-length(toml_suffix)]
    end

    toml = TOML.parsefile("$(seedname)$(toml_suffix)")
    toml["seed_name"] = seedname
    params = Configurations.from_dict(Wan.Parameters.InputParams, toml)

    check_inputs(params)
    wannierize(params)
end

main()

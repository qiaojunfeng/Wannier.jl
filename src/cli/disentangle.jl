#!/usr/bin/env julia


"""
Maximally localize a group of entangled bands.

# Args

- `seedname`: seedname for WIN/AMN/MMN/EIG files

# Options

- `-o, --output`: filename for output AMN. Default is `seedname.dis.amn`
"""
@cast function dis(seedname::String; output::Union{String,Nothing} = nothing)

    # seedname = "/home/jqiao/git/Wannier.jl/test/fixtures/silicon"
    if output === nothing
        output = basename(seedname) * ".dis.amn"
    end

    model = read_seedname(seedname)

    # Orthonormalize Amn, make it semiunitary
    model.A .= orthonorm_lowdin(model.A)

    # eV
    dis_froz_max = 10.0

    Wannier.set_frozen_win!(model, dis_froz_max)

    # if params.dis_froz_pao
    #     V = zeros(ComplexF64, data.num_bands, data.num_bands, data.num_kpts)
    #     if params.dis_froz_win || params.dis_froz_proj || params.dis_froz_num > 0
    #         # freeze both PAO space and user specified space
    #         for ik = 1:data.num_kpts
    #             amn0[:,:,ik], V[:,:,ik] = Wan.Disentangle.maxproj_froz(amn0[:,:,ik], data.frozen[:,ik])
    #         end
    #     else
    #         # freeze only PAO space
    #         for ik = 1:data.num_kpts
    #             amn0[:,:,ik], V[:,:,ik] = Wan.Disentangle.max_projectability(amn0[:,:,ik])
    #         end
    #     end
    #     for ik = 1:data.num_kpts
    #         for ib = 1:data.num_bvecs
    #             data.mmn[:,:,ib,ik] = V[:,:,ik]' * data.mmn[:,:,ib,ik] * V[:,:,data.kpbs[ib,ik]]
    #             # for band interpolation we need to rotate data.eig as well!
    #         end
    #     end
    #     # set it to num_wann, so during optimization the first num_wann rotated Bloch wfcs are always frozen
    #     params.dis_froz_num = data.num_wann
    #     params.dis_froz_win = false
    #     params.dis_froz_proj = false
    #     params.dis_froz_degen = false
    #     Wan.Disentangle.set_frozen_bands!(data, params)
    # end

    Amin = disentangle(model)

    # if dis_froz_pao
    #     for ik = 1:data.num_kpts
    #         A[:,:,ik] = V[:,:,ik] * A[:,:,ik]
    #     end
    # end

    write_amn(output, Amin)

    nothing
end


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

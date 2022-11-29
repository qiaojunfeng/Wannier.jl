"""
Maximally localize a group of entangled bands.

# Args

- `seedname`: seedname for `win`/`amn`/`mmn`/`eig` files

# Options

- `-o, --output=<str>`: filename for output `amn`. Default is `seedname.dis.amn`
- `-m, --maxiter=<int>`: max number of iterations. Default is `50`
"""
@cast function dis(seedname::String; output::String="", maxiter::Int=50)
    if output == ""
        output = basename(seedname) * ".dis.amn"
    end

    model = read_w90(seedname)

    # the dis_froz_max was set during reading win file, you can still change it here
    # dis_froz_max = 10.0  # eV
    # Wannier.set_frozen_win!(model, dis_froz_max)

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

    Umin = disentangle(model; max_iter=maxiter)

    # if dis_froz_pao
    #     for ik = 1:data.num_kpts
    #         U[:,:,ik] = V[:,:,ik] * U[:,:,ik]
    #     end
    # end

    write_amn(output, Umin)

    return nothing
end

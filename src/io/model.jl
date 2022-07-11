function write_model(seedname::String, model::Model)
    # Note I need to use basename! If seedname is absolute path the joinpath
    # will just return seedname.
    outname(suffix::String) = "$seedname.$suffix"

    write_eig(outname("eig"), model.E)

    kpb_k = model.bvectors.kpb_k
    kpb_b = model.bvectors.kpb_b
    write_mmn(outname("mmn"), model.M, kpb_k, kpb_b)

    write_amn(outname("amn"), model.A)

    return nothing
end

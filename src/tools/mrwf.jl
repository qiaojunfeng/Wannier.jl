using OrderedCollections
using TOML: parsefile

"""
    $(SIGNATURES)

Split valence and conduction Wannier functions into subgroups.

Usually start from a Wannierization of valence+conduction bands,
then this command split MLWFs into several subgroups.

# Arguments
- `model`: valence+conduction band `Model`, should be already disentangled
    and max localized.
- `indices`: a vector, where each element is a vector of Wannier band indices
    (of the input valence+conduction MLWFs) that belong to a target (isolated
    in energy) subgroup.
- `model_cubic`: for some system, the ``b``-vectors might not contain the 6
    cubic nearest neighbors. Therefore, we need a special `Model` that contains
    the overlaps for such cubic neighbors. If not provided, defaults to the
    1st argument `model`.

# Keyword Arguments
- `run_dis`: run disentangle on the `model`
- `run_optrot`: max localize w.r.t. single unitary matrix after parallel transport.
    Should further reduce the spread and much closer to the true max localization.
- `run_maxloc`: run a final max localize w.r.t. all kpoints.
    Should reach the true max localization.
"""
function mrwf(
        model::Model,
        indices::AbstractVector{T},
        model_cubic::Model = model;
        run_dis::Bool = false,
        run_optrot::Bool = false,
        run_maxloc::Bool = false,
    ) where {T <: AbstractVector}
    # Check input
    length(indices) > 0 || @error "indices must not be empty"

    # The valcond MLWFs
    run_dis && (model.gauges .= disentangle(model))

    # @info "Spread of input model"
    # show(omega(model))
    # println("\n")

    models_Us = split_wannierize(model_cubic, indices)
    if model != model_cubic
        models_Us = map(models_Us) do (m_cubic, U)
            m = transform_gauge(model, U)
            m.gauges .= m_cubic.gauges
            (m, U)
        end
    end

    for (i, (m, _)) in enumerate(models_Us)
        @info "Group $i after parallel transport:" omega(m)
        println("\n")

        if run_optrot
            @info "Run optimal rotation"
            println()
            W = opt_rotate(m)
            m.gauges .= merge_gauge(m.gauges, W)
        end

        if run_maxloc
            @info "Run max localization"
            println()
            m.gauges .= max_localize(m)
        end
    end

    return models_Us
end

"""
    $(SIGNATURES)

# Arguments
- `prefix`: prefix for W90 `win`/`amn`/`mmn`/`eig`/`chk` files. Will read
    `chk` file if exists.
- `outdirs`: a vector of output directories, where each element is the
    directory name for the corresponding subgroup. The files will be written
    in the form of `prefix.amn`, `prefix.mmn`, `prefix.eig`, and `prefix.chk`
    in the corresponding subgroup directory.
- `mmn_cubic`: filename for the `mmn` file that contains the 6 cubic nearest
    neighbors for [`parallel_transport`](@ref). If not provided, the overlap
    matrices from the `prefix.mmn` will be used instead.

# Keyword Arguments
- `chk`: path of chk file to get the unitary matrices for the `prefix.mmn`
    calculation, if `nothing`, read the initial gauges from `prefix.amn` file.
- `rot_unk`: generate `unk` files for subgroups, for plotting WFs
"""
function mrwf(
        prefix::AbstractString,
        indices::AbstractVector{T},
        outdirs::AbstractVector{S},
        mmn_cubic::Union{AbstractString, Nothing} = nothing;
        chk::Union{AbstractString, Nothing} = "$prefix.chk",
        rot_unk::Bool = false,
        kwargs...,
    ) where {T <: AbstractVector, S <: AbstractString}
    length(indices) == length(outdirs) ||
        error("length of indices and outdirs must be equal")

    if isnothing(chk)
        # Read initial gauges from amn file, do not use chk
        @warn "chk file not found, using gauges from amn file instead"
        model = read_w90(prefix)
    else
        # Read MLWF gauge from chk file
        model = read_w90_with_chk(prefix, chk)
    end

    println("Model will be split into $(length(indices)) groups")
    for (i, (idxs, od)) in enumerate(zip(indices, outdirs))
        println("  Group $i:")
        println("    indices: $idxs")
        println("    outdir : $od")
    end

    if isnothing(mmn_cubic)
        models_Us = mrwf(model, indices; kwargs...)
    else
        # Read mmn file for cubic neighbors
        # I assume the kpoints,
        mmn = read_mmn(mmn_cubic)
        M_cubic = mmn.M
        kpb_k_cubic = mmn.kpb_k
        kpb_G_cubic = mmn.kpb_G
        nbvecs = length(kpb_k_cubic[1])
        (nbvecs == 6) || error("number of b-vectors in mmn_cubic must be 6")
        (length(kpb_k_cubic) == n_kpoints(model)) || error(
            "number of kpoints in mmn_cubic must match the number of kpoints in model"
        )
        bvectors = Wannier.get_bvectors(
            model.recip_lattice, model.kpoints, kpb_k_cubic, kpb_G_cubic
        )
        # Just set weights to 0.0, they are not complete b-vectors
        weights = zeros(Float64, nbvecs)
        kstencil_cubic = Wannier.KspaceStencil{Float64}(
            model.recip_lattice,
            model.kgrid_size,
            model.kpoints,
            bvectors,
            weights,
            kpb_k_cubic,
            kpb_G_cubic,
        )
        model_cubic = Wannier.Model(model, kstencil_cubic, M_cubic)
        models_Us = mrwf(model, indices, model_cubic; kwargs...)
    end

    win = read_win("$prefix.win")
    win = OrderedDict(pairs(win))
    for k in [
            "num_bands",
            "dis_froz_proj",
            "dis_proj_min",
            "dis_proj_max",
            "dis_win_min",
            "dis_win_max",
            "dis_froz_min",
            "dis_froz_max",
            "projections",
            "auto_projections",
        ]
        pop!(win, k, nothing)
    end
    # Just use auto_projections as a placeholder
    win["auto_projections"] = true
    win["num_iter"] = 2000

    for (od, (m, U)) in zip(outdirs, models_Us)
        # This writes to the folder of prefix
        # od = isabspath(od) ? od : joinpath(dirname(prefix), od)
        mkpath(od)
        prefix_i = joinpath(od, basename(prefix))

        # Write files
        # The parallel transport gauge is written in the amn file
        write_w90(prefix_i, m)

        header =
            Wannier.default_header() *
            "    Gauge matrix from the original mmn/eig to the subgroup mmn/eig"
        write_amn("$(prefix_i)_split.amn", U; header)

        # Prepare win file with correct num_wann
        win["num_wann"] = n_wannier(m)
        write_win("$prefix_i.win", win)
    end

    # UNK files for plotting WFs
    if rot_unk
        dir = dirname(prefix)
        isempty(dir) && (dir = ".")
        # This writes to the folder of prefix
        # outdirs = [joinpath(dir, od) for od in outdirs]
        # Probably better just write to the `outdirs` specified by user
        split_unk(dir, [mU[2] for mU in models_Us], outdirs)
    end
    return nothing
end

"""
    $(SIGNATURES)

# Arguments
- `config_file`: config file for subgroups, e.g.
    ```toml
    [groups]
    indices = [[1, 2], [3, 4, 5, 6], [7, 8,]]
    outdirs = ["val_1", "val_2", "cond_3"]
    ```

# Keyword Arguments
See [`mrwf`](@ref) for more details.
"""
function mrwf(
        prefix::AbstractString,
        config_file::AbstractString,
        mmn_cubic::Union{AbstractString, Nothing} = nothing;
        kwargs...,
    )
    @info "reading config file: $config_file"
    groups = parsefile(config)["groups"]
    indices = groups["indices"]
    outdirs = groups["outdirs"]
    return mrwf(prefix, indices, outdirs, mmn_cubic; kwargs...)
end

"""
    $(SIGNATURES)

# Arguments
- `nval`: number of valence WFs.
- `outdir_val`: dirname for output valence `amn`/`mmn`/`eig`. Default is `val`
- `outdir_cond`: dirname for output conduction `amn`/`mmn`/`eig`. Default is `cond`

# Keyword Arguments
See [`mrwf`](@ref) for more details.
"""
function mrwf(
        prefix::AbstractString,
        nval::Integer,
        outdir_val::AbstractString = "val",
        outdir_cond::AbstractString = "cond",
        mmn_cubic::Union{AbstractString, Nothing} = nothing;
        kwargs...,
    )
    win = read_win(joinpath(prefix, ".win"))
    nwan = win["num_wann"]
    (0 < nval < nwan) || @error "nval must > 0 and < n_wannier"
    @info "number of valence WFs = $nval"

    indices = [[1:nval, (nval + 1):nwan]]
    outdirs = [[outdir_val, outdir_cond]]

    return mrwf(prefix, indices, outdirs, mmn_cubic; kwargs...)
end

"""
    $(SIGNATURES)

# Arguments
- `Usplits`: gauge matrices from `prefix_split.amn` for val/cond.
- `Umaxlocs`: gauge matrices from `prefix.chk` or `prefix_u.mat` files for val/cond.
"""
function merge_gauge(Usplits::AbstractVector, Umaxlocs::AbstractVector)
    Us = map(zip(Usplits, Umaxlocs)) do (Usplit, Umaxloc)
        #=
        # Usplit
        From `_split.amn` files, these are the gauge matrices mapping from the
        val+cond mmn/eig to the val or cond mmn/eig, i.e., including the
        dis+maxloc gauge of val+cond, and the gauge due to diagonalizing the
        val+cond Hamiltonian.

        # Umaxloc
        From the `chk` or `u.mat` file, these are the gauge matrices mapping
        from the val/cond mmn/eig to the final val/cond MLWFs, i.e., including
        the parallel transport gauge (which is written in the val/prefix.amn file)
        and the maximal localization gauge from W90.
        =#
        Wannier.merge_gauge(Usplit, Umaxloc)
    end
    Utot = map(zip(Us...)) do Uks
        # Concatenate along the n_wannier indices, at each kpoint
        reduce(hcat, Uks)
    end
    return Utot
end

"""
    $(SIGNATURES)

# Arguments
- `filename`: filename for the output `amn` file. Default is `mrwf.amn`.
- `split_amns`: a vector of filenames for the `_split.amn` files for each subgroup.
    The order of the filenames determines the order of the indices of MRWFs in
    the output total gauge matrix.
- `maxlocs`: a vector of filenames for the `chk` or `u.mat` files for each subgroup.
    You must ensure the order matches the order of `split_amns`.

# Keyword Arguments
- `umat`: if `true`, read the `prefix_u.mat` file instead of the `prefix.chk` file.
    Default is `false`.
"""
function merge_gauge(
        filename::AbstractString,
        split_amns::AbstractVector{<:AbstractString},
        maxlocs::AbstractVector{<:AbstractString};
        umat::Bool = false,
    )
    Usplits = map(split_amns) do f
        read_amn(f).A
    end
    Umaxlocs = map(maxlocs) do f
        if umat
            return WannierIO.read_u_mat(f)[1]
        else
            return WannierIO.gauge_matrices(read_chk(f))
        end
    end
    Utot = merge_gauge(Usplits, Umaxlocs)
    header =
        WannierIO.default_header() *
        "   Block-diagonal gauge matrix from the original mmn/eig to the final MRWFs"
    write_amn(filename, Utot; header)
    return Utot
end

"""
    $(SIGNATURES)

# Arguments
- `filename`: filename for the output `amn` file. Default is `prefix_mrwf.amn`.
- `prefix`: prefix for the `prefix_split.amn` and `prefix.chk` or `prefix_u.mat` files.
- `outdirs`: a vector of output directories for each group of isolated band manifold.
    Note the order determines the order of the indices of MRWFs in the output total
    gauge matrix.

# Keyword Arguments
- `umat`: if `true`, read the `prefix_u.mat` file instead of the `prefix.chk` file.
    Default is `false`.
"""
function merge_gauge(
        filename::AbstractString,
        prefix::AbstractString,
        outdirs::AbstractVector{<:AbstractString};
        umat::Bool = false,
    )
    split_amns = map(outdirs) do od
        joinpath(od, "$(prefix)_split.amn")
    end
    maxlocs = map(outdirs) do od
        if umat
            return joinpath(od, "$(prefix)_u.mat")
        else
            return joinpath(od, "$(prefix).chk")
        end
    end
    return merge_gauge(filename, split_amns, maxlocs; umat)
end

function merge_gauge(
        prefix::AbstractString, outdirs::AbstractVector{<:AbstractString}; kwargs...
    )
    filename = "$(prefix)_mrwf.amn"
    return merge_gauge(filename, prefix, outdirs; kwargs...)
end

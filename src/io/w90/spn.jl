export read_chk_spn, read_w90_tb_chk_spn

"""
    $(SIGNATURES)

Read `prefix.chk` and `prefix.spn` to construct Wannier-gauge k-space spin operator.

# Arguments
- `prefix`: the prefix of `prefix.chk` and `prefix.spn`

# Return
- `kpoints`: fractional kpoint coordinates
- `S`: Wannier-gauge spin operator in k-space
"""
function read_chk_spn(prefix::AbstractString; chk = "$prefix.chk", spn = "$prefix.spn")
    spn_dat = read_spn(spn)
    spin_vecs = map(zip(spn_dat.Sx, spn_dat.Sy, spn_dat.Sz)) do (x, y, z)
        MVec3.(x, y, z)
    end
    # spn file is in Bloch gauge, need to read chk to convert to Wannier gauge
    chk = read_chk(chk)
    gauges = WannierIO.gauge_matrices(chk)
    S = map(zip(gauges, spin_vecs)) do (Uₖ, Sₖ)
        Uₖ' * Sₖ * Uₖ
    end
    return chk.kpoints, S
end

"""
    $(SIGNATURES)

Read `prefix.chk` and `prefix.spn` to construct tight-binding models.

# Arguments
- `Rspace`: R-space domain

# Return
- `spin`: spin operator in R-space

!!! note
    This will call [`simplify`](@ref) to absorb the R-vector degeneracies and
    T-vectors into the operator, leading to faster interpolations.
"""
function read_chk_spn(prefix::AbstractString, Rspace::AbstractRspace; kwargs...)
    kpoints, S_k = read_chk_spn(prefix; kwargs...)
    # fourier transform to Rspace
    S_R = fourier(kpoints, S_k, Rspace)
    bare_Rspace, bare_spin = simplify(Rspace, S_R)
    spin = TBSpin(bare_Rspace, bare_spin)
    return spin
end

"""
    $(SIGNATURES)

Read `prefix_tb.dat`, `prefix_wsvec.dat`, `prefix.chk` and `prefix.spn` to
construct tight-binding models.

# Arguments
- `prefix`: the prefix of `prefix_tb.dat` and `prefix_wsvec.dat`

# Keyword Arguments
- `chk`: the path to `prefix.chk`, default is `"\$prefix.chk"`
- `spn`: the path to `prefix.spn`, default is `"\$prefix.spn"`

# Return
- a [`TBHamiltonian`](@ref)
- a [`TBPosition`](@ref)
- a [`TBSpin`](@ref)

!!! note
    This will call [`simplify`](@ref) to absorb the R-vector degeneracies and
    T-vectors into the operator, leading to faster interpolations.
"""
function read_w90_tb_chk_spn(prefix::AbstractString; kwargs...)
    dat = _raw_read_w90_tb(prefix)
    Rspace = dat.Rspace
    H = dat.hamiltonian
    pos = dat.position

    bare_Rspace, bare_H, bare_pos = simplify(Rspace, H, pos)
    hamiltonian = TBHamiltonian(bare_Rspace, bare_H)
    position = TBPosition(bare_Rspace, bare_pos)

    spin = read_chk_spn(prefix, Rspace; kwargs...)

    return (; hamiltonian, position, spin)
end

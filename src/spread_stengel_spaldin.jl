export StengelSpaldinPenalty

@doc raw"""
    struct StengelSpaldinPenalty

The Stengel-Spaldin (SS) spread functional, i.e. wannier90's `use_ss_functional`.

The SS functional takes the average over ``\bm{k}`` *before* the logarithm,
where the Marzari-Vanderbilt (MV) functional of [`Spread`](@ref) takes it after:

```math
S_{n\bm{b}} = \frac{1}{N_k} \sum_{\bm{k}} M_{nn}^{\bm{k},\bm{b}}, \qquad
\langle \bm{r} \rangle_n = -\sum_{\bm{b}} w_{\bm{b}} \bm{b} \, \Im \log S_{n\bm{b}}
```

Only the diagonal part of the spread changes; ``\Omega_{\mathrm{I}}`` and
``\Omega_{\mathrm{OD}}`` are identical to MV. ``\Omega_{\mathrm{D}}`` becomes the
variance over ``\bm{k}`` of the diagonal overlaps, and is therefore manifestly
non-negative:

```math
\Omega_{\mathrm{D}} = \sum_{n\bm{b}} w_{\bm{b}} \left(
    \frac{1}{N_k} \sum_{\bm{k}} | M_{nn}^{\bm{k},\bm{b}} |^2
    - | S_{n\bm{b}} |^2 \right)
```

The total collapses to ``\Omega = \sum_{n\bm{b}} w_{\bm{b}}
(1 - | S_{n\bm{b}} |^2)``, which [`omega!`](@ref) reports per WF in `ω`.

Both functionals converge to the same result on an infinitely fine
``\bm{k}``-grid, and both are translationally invariant, but only SS is size
consistent. See Stengel and Spaldin, Phys. Rev. B **73**, 075121 (2006).

Use it in place of [`SpreadPenalty`](@ref):

```julia
p = StengelSpaldinPenalty(model.kstencil)
model.gauges .= max_localize(p, model)
```

!!! note

    Because ``S_{n\bm{b}}`` sums over ``\bm{k}`` at fixed ``\bm{b}``, slot `ib`
    must denote the same ``\bm{b}``-vector at every kpoint. Rather than
    reordering the overlaps, this struct stores the permutation `nnord` and
    indexes through it, which is what wannier90's `kmesh_info%nnord` does.

# Fields
$(FIELDS)
"""
struct StengelSpaldinPenalty{T <: Real} <: AbstractPenalty
    """the kpoint stencil this penalty was built from. `nnord` indexes its
    b-vector slots, so the penalty is only valid for this stencil."""
    kstencil::KspaceStencil{T}

    """`nnord[ik][ib]` is the slot at kpoint `ik` holding
    `kstencil.bvectors[ib]`, i.e. wannier90's `kmesh_info%nnord`"""
    nnord::Vector{Vector{Int}}
end

"""
    $(SIGNATURES)

Build the b-vector orderings the SS functional needs from a `KspaceStencil`.

Errors if the stencil is not closed under ``\\bm{b} \\to -\\bm{b}``, which the
gradient requires.
"""
function StengelSpaldinPenalty(kstencil::KspaceStencil{T}) where {T}
    bvectors = kstencil.bvectors
    inv_recip_lattice = inv(reciprocal_lattice(kstencil))
    bvectors_frac = map(b -> inv_recip_lattice * b, bvectors)

    nnord = map(1:n_kpoints(kstencil)) do ik
        bvecs = get_bvectors(kstencil, ik; fractional = true)
        slots = map(enumerate(bvectors_frac)) do (ib, b)
            slot = findfirst(x -> isapprox(x, b; atol = 1.0e-6), bvecs)
            isnothing(slot) && error("bvector $ib not found at kpoint $ik")
            return slot
        end
        # `omega_grad!` indexes through `nnord` under `@inbounds`
        isperm(slots) || error("bvector mapping at kpoint $ik is not a permutation")
        return slots
    end

    # The gradient folds the -b term onto +b using S[n, -b] = conj(S[n, b]), so it
    # needs every -b present. The functional alone would not.
    for (ib, b) in enumerate(bvectors_frac)
        any(x -> isapprox(x, -b; atol = 1.0e-6), bvectors_frac) || error(
            "no -b partner for bvector $ib; the Stengel-Spaldin gradient needs " *
            "a stencil closed under b -> -b"
        )
    end

    return StengelSpaldinPenalty{T}(kstencil, nnord)
end

StengelSpaldinPenalty(model::Model) = StengelSpaldinPenalty(model.kstencil)

n_bvectors(p::StengelSpaldinPenalty) = length(p.kstencil.bvectors)

"""
    $(SIGNATURES)

Check that the penalty was built from the `KspaceStencil` it is being used with.

`nnord` indexes b-vector slots of the stencil the penalty was built from, while
the caller's `UtMU` is built from the stencil handed to `omega!`. If the two
differ the b-vector orderings are mixed and the result is silently wrong.
`kpb_k` is what has to agree: `reorder` returns a stencil sharing the same
`bvectors` and `bweights` but with permuted `kpb_k`, so identity on those two
alone would not catch it.
"""
function check_stencil(p::StengelSpaldinPenalty, kstencil::KspaceStencil)
    (p.kstencil === kstencil || p.kstencil.kpb_k == kstencil.kpb_k) || error(
        "the StengelSpaldinPenalty was built from a different KspaceStencil; " *
        "rebuild it with StengelSpaldinPenalty(kstencil)"
    )
    return nothing
end

"""
    $(SIGNATURES)

Compute ``S_{n\\bm{b}} = \\langle (U^\\dagger M U)_{nn} \\rangle_{\\bm{k}}``,
the k-average taken before the logarithm.

Returns a `n_wann * n_bvectors` matrix.
"""
function compute_ss_overlaps(p::StengelSpaldinPenalty, UtMU::AbstractVector, nwann::Integer)
    nkpts = length(UtMU)
    nbvecs = n_bvectors(p)
    S = zeros(eltype(UtMU[1][1]), nwann, nbvecs)

    for ik in 1:nkpts
        Nk = UtMU[ik]
        nnord = p.nnord[ik]
        for ib in 1:nbvecs
            Nkb = Nk[nnord[ib]]
            for n in 1:nwann
                S[n, ib] += Nkb[n, n]
            end
        end
    end

    S ./= nkpts
    return S
end

function omega!(
        p::StengelSpaldinPenalty, cache::Cache{FT}, bvectors::KspaceStencil{FT}, M
    ) where {FT <: Real}
    check_stencil(p, bvectors)
    UtMU = cache.UtMU
    nwann = n_wann(cache)
    nkpts = n_kpts(cache)
    nbvecs = n_bvectors(p)

    # S[n, ib] = ⟨(U† M U)_nn⟩_k, the k-average taken *before* the logarithm.
    # This is what distinguishes SS from MV, and why slot `ib` must denote the
    # same b-vector at every kpoint.
    S = zeros(Complex{FT}, nwann, nbvecs)
    # ⟨|M_nn|²⟩_k, the other half of the k-variance
    absM² = zeros(FT, nwann, nbvecs)

    ΩI::FT = 0.0
    ΩOD::FT = 0.0

    for ik in 1:nkpts
        Nk = UtMU[ik]
        nnord = p.nnord[ik]
        for ib in 1:nbvecs
            Nkb = Nk[nnord[ib]]
            wᵇ = p.kstencil.bweights[ib]

            # ΩI and ΩOD are unchanged from MV
            ts = zero(FT)
            ts2 = zero(FT)
            for i in axes(Nkb, 2)
                for j in axes(Nkb, 1)
                    a2 = abs2(Nkb[j, i])
                    ts += a2
                    if i == j
                        S[i, ib] += Nkb[i, i]
                        absM²[i, ib] += a2
                    else
                        ts2 += a2
                    end
                end
            end

            ΩI += wᵇ * (nwann - ts)
            ΩOD += wᵇ * ts2
        end
    end

    S ./= nkpts
    absM² ./= nkpts
    ΩI /= nkpts
    ΩOD /= nkpts

    # ΩD is the k-variance of the diagonal overlaps, hence ≥ 0
    ΩD::FT = 0.0
    for ib in 1:nbvecs
        wᵇ = p.kstencil.bweights[ib]
        for n in 1:nwann
            ΩD += wᵇ * (absM²[n, ib] - abs2(S[n, ib]))
        end
    end

    r = cache.r
    fill!(r, zero(eltype(r)))
    ω = zeros(FT, nwann)
    for ib in 1:nbvecs
        wb_b = p.kstencil.bweights[ib] * p.kstencil.bvectors[ib]
        wᵇ = p.kstencil.bweights[ib]
        for n in 1:nwann
            r[n] -= imaglog(S[n, ib]) * wb_b
            ω[n] += wᵇ * (1 - abs2(S[n, ib]))
        end
    end

    Ω̃ = ΩOD + ΩD
    Ω = ΩI + Ω̃
    return Spread(Ω, ΩI, ΩOD, ΩD, Ω̃, ω, copy(r), :StengelSpaldin)
end

function omega(p::StengelSpaldinPenalty, bvectors::KspaceStencil, M, U)
    cache = Cache(bvectors, M, U)
    compute_MU_UtMU!(cache, bvectors, M, U)
    return omega!(p, cache, bvectors, M)
end

function omega(p::StengelSpaldinPenalty, bvectors::KspaceStencil, M, X, Y)
    return omega(p, bvectors, M, X_Y_to_U(X, Y))
end

"""
    $(SIGNATURES)

Stengel-Spaldin spread of a [`Model`](@ref), optionally for a given gauge.
"""
function omega(
        p::StengelSpaldinPenalty, model::Model, gauges::AbstractVector = model.gauges
    )
    return omega(p, model.kstencil, model.overlaps, gauges)
end

@doc raw"""
    $(SIGNATURES)

Gradient of the Stengel-Spaldin spread w.r.t. the gauge.

With ``G = 2 \partial \Omega / \partial \bar{U}``, which is the convention of
[`omega_grad!`](@ref) for the MV functional, and
``\Omega = \sum_{n\bm{b}} w_{\bm{b}} (1 - | S_{n\bm{b}} |^2)``:

```math
G_{\bm{k}}[m,n] = 2 \frac{\partial \Omega}{\partial \bar{U}_{\bm{k}}[m,n]}
    = -\frac{4}{N_k}
    \sum_{\bm{b}} w_{\bm{b}} \, S^{*}_{n\bm{b}}
    \left( M^{\bm{k},\bm{b}} U_{\bm{k}+\bm{b}} \right)[m,n]
```

Differentiating ``|S_{n\bm{b}}|^2`` w.r.t. ``\bar{U}_{\bm{k}}`` gives two
terms, one through the ``S^{*}`` factor at ``+\bm{b}`` and one through the ``S``
factor at ``-\bm{b}``. They combine because ``M^{\bm{k},\bm{b}\dagger} = M^{\bm{k}+\bm{b},-\bm{b}}`` turns
the ``-\bm{b}`` term into another `MU` at ``\bm{k}``, and because
``S_{n,-\bm{b}} = S^{*}_{n\bm{b}}`` — which is why the stencil must be closed
under ``\bm{b} \to -\bm{b}``.

!!! note

    This is *not* wannier90's `wann_domega` expression. That routine computes
    ``\partial \Omega / \partial W`` for a unitary generator
    (``U \to U e^{W}``), which carries an extra antihermitian structure; here
    the unconstrained ``\partial \Omega / \partial U`` is wanted, matching what
    [`max_localize`](@ref) hands to `Optim`.
"""
function omega_grad!(
        p::StengelSpaldinPenalty, cache::Cache{FT}, bvectors::KspaceStencil{FT}, M
    ) where {FT <: Real}
    check_stencil(p, bvectors)
    G = cache.G
    fill!(G, 0)
    UtMU = cache.UtMU
    MU = cache.MU

    nbands, nwann, nkpts = size(G)
    nbvecs = n_bvectors(p)

    S = compute_ss_overlaps(p, UtMU, nwann)

    @inbounds for ik in 1:nkpts
        MUk = MU[ik]
        nnord = p.nnord[ik]
        Gk = view(G, :, :, ik)
        for ib in 1:nbvecs
            wᵇ = p.kstencil.bweights[ib]
            MUkb = MUk[nnord[ib]]
            for n in 1:nwann
                c = -4 * wᵇ * conj(S[n, ib])
                for m in 1:nbands
                    Gk[m, n] += c * MUkb[m, n]
                end
            end
        end
    end

    G ./= nkpts
    return G
end

function omega_grad(p::StengelSpaldinPenalty, bvectors::KspaceStencil, M, U)
    cache = Cache(bvectors, M, U)
    compute_MU_UtMU!(cache, bvectors, M, U)
    return omega_grad!(p, cache, bvectors, M)
end

function omega_grad(p::StengelSpaldinPenalty, bvectors::KspaceStencil, M, X, Y, frozen)
    U = X_Y_to_U(X, Y)
    G = omega_grad(p, bvectors, M, U)
    return GU_to_GX_GY(G, X, Y, frozen)
end

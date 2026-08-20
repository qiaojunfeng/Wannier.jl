# Benchmark of symmetry-constrained (SAWF) localization: full-mesh baseline
# versus the Level-1 (expand + pull back) and Level-2 (IBZ transport) paths.
#
# Uses the Ge4Ru4 dataset (IBZ .iamn/.immn/.ieig + .isym); point RERUN_DIR at
# a directory with these files. The full-mesh overlaps are unfolded from the
# IBZ ones (validated elsewhere against pw2wannier90 + wannier90).
#
#   julia --project dev/benchmark_symmetric.jl

using Wannier, WannierIO, LinearAlgebra, Printf

RERUN_DIR = get(ENV, "RERUN_DIR", "/scratch/git/WannierDatasets/datasets/Ge4Ru4/rerun")
prefix = "Ge4Ru4"

# ---- load & build ----------------------------------------------------------
nnkp = read_nnkp(joinpath(RERUN_DIR, "$prefix.nnkp"))
ks0 = Wannier.KspaceStencil(
    nnkp["recip_lattice"], nnkp["kpoints"], nnkp["kpb_k"], nnkp["kpb_G"]
)
isym = read_isym(joinpath(RERUN_DIR, "$prefix.isym"))
Wannier.rescale!(isym.littlegroup_reps)
centers = [p.center for p in nnkp["projections"]]

tsc = @elapsed sc = Wannier.symmetry_constraint(ks0, isym, centers)
ks = Wannier.globalize_stencil(ks0)

mmn_i = read_mmn(joinpath(RERUN_DIR, "$prefix.immn"))
tunf = @elapsed Mf = Wannier.unfold_overlaps_cached(mmn_i.M, sc)
Ai = read_amn(joinpath(RERUN_DIR, "$prefix.iamn")).A
Ei = read_eig(joinpath(RERUN_DIR, "$prefix.ieig"))
win = read_win(joinpath(RERUN_DIR, "$prefix.win"))
Ef = Wannier.unfold_eigvals(Ei, [collect(t) for t in sc.fbz2ibz])
frozen = Wannier.get_frozen_bands(Ef, win["dis_froz_max"])
Af = Wannier.expand_gauges(Wannier.project_covariant(Ai, sc), sc)
atom_positions = [p.second for p in win["atoms_frac"]]
atom_labels = map(x -> string(x.first), win["atoms_frac"])
model = Wannier.Model(
    win["unit_cell_cart"], atom_positions, atom_labels, ks, Mf, Af, Ef, frozen
)

@printf(
    "Ge4Ru4: %d bands, %d WFs, %d FBZ / %d IBZ kpoints, %d b-vectors, %d symops\n",
    sc.nbands, sc.nwann, sc.nk_fbz, sc.nk_ibz, sc.nbvecs, length(isym.symops)
)
@printf("setup: constraint tables %.2fs, M unfolding %.3fs\n\n", tsc, tunf)

# ---- per-evaluation cost ----------------------------------------------------
ws1 = Wannier.SymmetricWorkspace(model, sc)
ws2 = Wannier.SymmetricWorkspace2(Ef, frozen, sc)
X, Y = Wannier.U_to_X_Y(Ai, ws1.frozen_ibz)
xy = Wannier.X_Y_to_XY(X, Y)
G1, G2 = zero(xy), zero(xy)

prob = Wannier.Problem(Wannier.Variance(), model)
fgfull = Wannier._make_fg!(prob)
x0full = Wannier.initial_x(prob.layout, model)
Gfull = zero(x0full)

mintime(f, n = 20) = minimum((f(); @elapsed f()) for _ in 1:n)

t_full_fg = mintime(() -> fgfull(1.0, Gfull, x0full))
t_full_f = mintime(() -> fgfull(1.0, nothing, x0full))
t_l1_fg = mintime(() -> Wannier.symmetric_fg1!(1.0, G1, xy, model, sc, ws1))
t_l1_f = mintime(() -> Wannier.symmetric_fg1!(1.0, nothing, xy, model, sc, ws1))
t_l2_fg = mintime(() -> Wannier.symmetric_fg2!(1.0, G2, xy, mmn_i.M, sc, ws2))
t_l2_f = mintime(() -> Wannier.symmetric_fg2!(1.0, nothing, xy, mmn_i.M, sc, ws2))

frozen_ibz = frozen[:, sc.ibz2fbz]
tsb = @elapsed sb = Wannier.schur_basis(sc, frozen_ibz)
U0i = Wannier.project_covariant(Ai, sc)
xs = Wannier.schur_initial_x(U0i, sb)
gs = zero(xs)
ws2s = Wannier.SymmetricWorkspace2(Ef, frozen, sc)
fg_schur! = function (F, G, x)
    Wannier.schur_decode!(ws2s.U_ibz, x, sb)
    Ω = Wannier._fg2_core!(F, G === nothing ? nothing : ws2s.G_ibz, mmn_i.M, sc, ws2s)
    G === nothing || Wannier.schur_encode_gradient!(G, ws2s.G_ibz, x, sb)
    return Ω
end
t_sch_fg = mintime(() -> fg_schur!(1.0, gs, xs))
t_sch_f = mintime(() -> fg_schur!(1.0, nothing, xs))

println("per-evaluation wall time (min of 20):")
@printf("  %-28s %8.1f ms   value-only %8.1f ms\n", "full-mesh (unconstrained)", 1e3 * t_full_fg, 1e3 * t_full_f)
@printf("  %-28s %8.1f ms   value-only %8.1f ms\n", "Level 1 (constrained)", 1e3 * t_l1_fg, 1e3 * t_l1_f)
@printf("  %-28s %8.1f ms   value-only %8.1f ms\n", "Level 2 (constrained)", 1e3 * t_l2_fg, 1e3 * t_l2_f)
@printf("  %-28s %8.1f ms   value-only %8.1f ms\n", "Level 2 + Schur blocks", 1e3 * t_sch_fg, 1e3 * t_sch_f)
@printf("  Level-2 speedup: %.2fx vs full mesh, %.2fx vs Level 1\n", t_full_fg / t_l2_fg, t_l1_fg / t_l2_fg)
@printf("  Schur speedup:   %.2fx vs full mesh\n", t_full_fg / t_sch_fg)
@printf(
    "  parameters: XY %d complex, Schur %d complex (%.1fx fewer; basis setup %.1fs)\n",
    (sc.nwann^2 + sc.nbands * sc.nwann) * sc.nk_ibz, sb.nx,
    (sc.nwann^2 + sc.nbands * sc.nwann) * sc.nk_ibz / sb.nx, tsb,
)
# corepresentation (anti-unitary) block statistics
allblk = reduce(vcat, sb.blocks)
nfall = sum(
    count(b -> b.akind == 0, blks)
        for (iki, blks) in enumerate(sb.blocks) if sb.aop[iki] !== nothing;
    init = 0,
)
@printf(
    "  anti-unitary classes: %d derived from pairing partners, %d soft-average fallback\n",
    count(b -> b.akind == 2, allblk), nfall,
)
# Hermiticity-pair halving coverage of the Level-2 pass 0 (only the 2-cycle
# pairs of the partner map are derived; see `_fg2_core!`)
pkey(iki, ibi) = (iki - 1) * sc.nbvecs + ibi
pof(iki, ibi) = (sc.ikb[ibi, iki], sc.ibi_of[sc.opp_b[ibi], sc.ikpb_fbz[ibi, iki]])
nderived = count(
    p -> pkey(pof(p...)...) < pkey(p...) && pof(pof(p...)...) == p,
    ((iki, ibi) for iki in 1:sc.nk_ibz, ibi in 1:sc.nbvecs),
)
@printf(
    "  Level-2 pass 0: %d of %d IBZ pairs derived via Hermiticity pairs\n\n",
    nderived, sc.nbvecs * sc.nk_ibz,
)

# equivalence at the shared starting point
Ω1 = Wannier.symmetric_fg1!(1.0, G1, xy, model, sc, ws1)
Ω2 = Wannier.symmetric_fg2!(1.0, G2, xy, mmn_i.M, sc, ws2)
@printf("L1/L2 equivalence: |ΔΩ| = %.2e, rel grad diff = %.2e\n\n", abs(Ω1 - Ω2), norm(G1 - G2) / norm(G1))

# ---- full optimizations ------------------------------------------------------
NITER = parse(Int, get(ENV, "NITER", "100"))
ts = @elapsed (Us, Usi) = Wannier.localize_symmetric(model, mmn_i.M, sc; level = 2, schur = true, max_iter = NITER)
t2 = @elapsed (U2, U2i) = Wannier.localize_symmetric(model, mmn_i.M, sc; level = 2, max_iter = NITER)
t1 = @elapsed (U1, _) = Wannier.localize_symmetric(model, mmn_i.M, sc; level = 1, max_iter = NITER)
t0 = @elapsed U0 = Wannier.localize(model; max_iter = NITER)

s2 = Wannier.spread(model.kstencil, model.overlaps, U2)
s1 = Wannier.spread(model.kstencil, model.overlaps, U1)
s0 = Wannier.spread(model.kstencil, model.overlaps, U0)

println("optimization ($NITER LBFGS iterations):")
@printf("  %-28s Ω = %.6f Å²   %6.1f s\n", "full-mesh MLWF", s0.Ω, t0)
@printf("  %-28s Ω = %.6f Å²   %6.1f s\n", "SAWF Level 1", s1.Ω, t1)
@printf("  %-28s Ω = %.6f Å²   %6.1f s\n", "SAWF Level 2", s2.Ω, t2)
ss = Wannier.spread(model.kstencil, model.overlaps, Us)
@printf("  %-28s Ω = %.6f Å²   %6.1f s\n", "SAWF Level 2 + Schur", ss.Ω, ts)
@printf("  covariance residual: L2 %.2e, Schur %.2e\n\n",
    Wannier.covariance_residual(U2i, sc), Wannier.covariance_residual(Usi, sc))

# ---- spread symmetry ---------------------------------------------------------
# WF orbits from the joint nonzero pattern of the orbital representations
orbit_of = collect(1:sc.nwann)
for rep in isym.orbital_reps, n in 1:sc.nwann, m in 1:sc.nwann
    if abs(rep.D[m, n]) > 1e-6
        r = min(orbit_of[m], orbit_of[n])
        orbit_of[orbit_of .== orbit_of[m]] .= r
        orbit_of[orbit_of .== orbit_of[n]] .= r
    end
end
println("per-orbit spread uniformity (max - min within symmetry orbit):")
for r in unique(orbit_of)
    ws = findall(orbit_of .== r)
    @printf(
        "  orbit %-12s  SAWF %.2e   Schur %.2e   MLWF %.2e\n",
        string(ws[1]) * "-" * string(ws[end]),
        maximum(s2.ω[ws]) - minimum(s2.ω[ws]),
        maximum(ss.ω[ws]) - minimum(ss.ω[ws]),
        maximum(s0.ω[ws]) - minimum(s0.ω[ws]),
    )
end

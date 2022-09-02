# Algorithms

Here are some detailed explanations of the algorithms used in this package,
however, it might still be too concise, please refer to the references for full discussion.

## Wannierization

- maximal localization for isolated bands, e.g. insulators
  - different from[^MV97], optimize on unitary matrix manifolds (adaptation of [^DLL19] to isolated bands)
- disentanglement for entangled bands, e.g. metal
  - different from[^SMV01], optimize on Stiefel manifolds[^DLL19]
- parallel transport gauge[^GLS19]
  - you can further improve the spread by optimization w.r.t. a single rotation matrix[^QMP21]
- split valence and conduction WFs from a valence + conduction calculation[^QMP21]
  - as a by-product, automated initial projection for valence or conduction WFs
  - for the initial projection of valence + conduction calculation, you can start with either
    conventional spdf projection, SCDM[^DL18], or an automated projection and disentanglement
    from pseudopotential orbitals[^QPM21]
  - different from SCDM, the valence+conduction manifold is chosen by the valence+conduction calculation,
    instead of SCDM μ and σ. Moreover, works in reciprocal space thus more memory-efficient
- constrain WF center for max localization or disentanglement[^QMP21]
  - similar to[^WLPMM14], add an Lagrange multiplier term to spread functional, but optimize
    on matrix manifolds, and applying to both max localization and disentanglement
    (whereas in [^WLPMM14] the center is constrained only during max localization)

## Interpolation

Two algorithms:

- Wigner-Seitz (WS) interpolation
- Minimal-distance replica selection (MDRS) method

for band structure along a kpath or on a grid.

## References

[^MV97]: Marzari, N. & Vanderbilt, D.,
    Maximally localized generalized Wannier functions for composite energy bands,
    [Physical Review B, 1997, 56, 12847](https://doi.org/10.1103/physrevb.56.12847)
[^MMYSV12]: Marzari, N.; Mostofi, A. A.; Yates, J. R.; Souza, I. & Vanderbilt, D.,
    Maximally localized Wannier functions: Theory and applications,
    [Reviews of Modern Physics, 2012, 84, 1419](https://doi.org/10.1103/revmodphys.84.1419)
[^SMV01]: Souza, I.; Marzari, N. & Vanderbilt, D.,
    Maximally localized Wannier functions for entangled energy bands,
    [Physical Review B, 2001, 65, 035109](https://doi.org/10.1103/physrevb.65.035109)
[^DLL19]: Damle, A.; Levitt, A. & Lin, L.,
    Variational Formulation for Wannier Functions with Entangled Band Structure,
    [Multiscale Modeling & Simulation, 2019, 17, 167](https://doi.org/10.1137/18m1167164)
[^GLS19]: Gontier, D.; Levitt, A. & Siraj-dine, S.,
    Numerical construction of Wannier functions through homotopy,
    [Journal of Mathematical Physics, 2019, 60, 031901](https://doi.org/10.1063/1.5085753)
[^QPM21]: Qiao, J.; Pizzi, G. & Marzari, N.,
    Projectability disentanglement for accurate high-throughput Wannierization,
    xxx
[^QMP21]: Qiao, J.; Marzari, N. & Pizzi, G.,
    Automated separate Wannierization for valence and conduction manifolds,
    xxx
[^DL18]: Damle, A. & Lin, L.,
    Disentanglement via Entanglement: A Unified Method for Wannier Localization
    [Multiscale Modeling & Simulation, 2018, 16, 1392](https://doi.org/10.1137/17m1129696)
[^WLPMM14]: Wang, R.; Lazar, E. A.; Park, H.; Millis, A. J. & Marianetti, C. A.,
    Selectively localized Wannier functions,
    [Physical Review B, 2014, 90, 165125](https://doi.org/10.1103/physrevb.90.165125)

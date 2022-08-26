# Normalization convention of WFs

There are several different conventions for Fourier transform, this page describes
the normalization convention used in `Wannier.jl`. In short, we adopt the same convention
as the Section II 3 of [^RMP].

For a unit cell having volume ``V``, the Bloch functions ``| \psi \rangle`` are periodic
in a super cell which is ``N`` times the unit cell (thus having volume ``NV``),
where ``N`` is determined by the discretization of the reciprocal space, i.e. the total
number of kpoints.

## Continuous case

If the super cell is large enough (``N \to \infty``), then the kpoint ``\bm{k}`` is a
continuous variable in the reciprocal cell having volume ``\frac{(2\pi)^3}{V}``. The inner
product between bra and ket ``\langle f | g \rangle`` are integrated over the the ``NV``
super cell that is ``(-\infty, \infty)`` in the continuous case. If we normalize the Bloch
functions over one unit cell, that is,

```math
\langle \psi_{n \bm{k}} | \psi_{n \bm{k}} \rangle_V
= \int_V | u_{n \bm{k}}(\bm{r}) |^2 \bm{r} = 1,
```

where ``\psi_{n \bm{k}} = \exp(i \bm{k} \bm{r}) u_{n \bm{k}}(\bm{r})``, then the inner
product between bra and ket is

```math
\langle \psi_{n \bm{k}} | \psi_{m \bm{k^\prime}} \rangle
= \frac{(2\pi)^3}{V} \delta_{n m} \delta(\bm{k} - \bm{k^\prime}).
```

This can be proven, e.g. in 1D with periodicity ``\bm{R}``, by

```math
\begin{aligned}
& \langle \psi_{n \bm{k}} | \psi_{m \bm{k^\prime}} \rangle
= \int_{-\infty}^{+\infty} \psi_{n \bm{k}}^* \psi_{m \bm{k^\prime}} d\bm{r} \\
= & \int_{-\infty}^{+\infty} \exp(i (\bm{k^\prime} - \bm{k}) \bm{r})
    u_{n \bm{k}}^*(\bm{r}) u_{m \bm{k^\prime}}(\bm{r}) d\bm{r} \\
= & \sum_{l = -\infty}^{+\infty} \int_{l \bm{R}}^{(l+1) \bm{R}} d(\bm{r} + l \bm{R})
    \exp(i (\bm{k^\prime} - \bm{k}) (\bm{r} + l \bm{R}))
    u_{n \bm{k}}^*(\bm{r} + l \bm{R}) u_{m \bm{k^\prime}}(\bm{r} + l \bm{R}) \\
= & \sum_{l = -\infty}^{+\infty} \exp(i (\bm{k^\prime} - \bm{k}) l \bm{R})
    \int_{0}^{\bm{R}} d\bm{r} \exp(i (\bm{k^\prime} - \bm{k}) \bm{r})
    u_{n \bm{k}}^*(\bm{r}) u_{m \bm{k^\prime}}(\bm{r}) d\bm{r} \\
= & \frac{2\pi}{\bm{R}} \operatorname{\text{Ш}}_{\frac{2\pi}{\bm{R}}}(\bm{k^\prime} - \bm{k})
    \int_{0}^{\bm{R}} \psi_{n \bm{k}}^* \psi_{m \bm{k^\prime}} d\bm{r} \\
= & \frac{2\pi}{\bm{R}} \delta(\bm{k} - \bm{k^\prime}) \delta_{n m}.
\end{aligned}
```

In the third line of the above equation, we replace the real space coordinate ``\bm{r}``
by ``\bm{r} + l \bm{R}``, such that ``\bm{r} \in [0, \bm{R})``. In the fourth line,
we utilize the real space periodicity of ``u_{n \bm{k}}``.
In the second last line, the [Dirac comb](https://en.wikipedia.org/wiki/Dirac_comb)
``\operatorname{\text{Ш}}_{T}(t)`` having periodicity ``T`` can be expanded as a
[Fourier series](https://en.wikipedia.org/wiki/Fourier_series)

```math
\operatorname {\text{Ш}} _{T}(t) =
{\frac {1}{T}} \sum _{n=-\infty }^{\infty }e^{i2\pi n{\frac {t}{T}}}.
```

In the last line, the Bloch function ``| \psi \rangle`` is periodic in ``\bm{k}``
thus the Dirac comb is reduced to a Dirac delta function
``\delta(\bm{k} - \bm{k^\prime})``, and of course the normalization of ``| \psi \rangle``
over one unit cell is used in the last line.

Thus, the Bloch functions are orthogonal w.r.t. ``\bm{k}`` only in the **super cell**,
they are **not** orthogonal when only integrating over **one unit cell**.

## Discretized case

For finite sampling of reciprocal space, e.g. a 1D cell of periodicity ``\bm{R}``,
sampled by ``N`` kpoints, then
``\bm{k} = \frac{2\pi}{\bm{R}}\frac{n}{N}, n = 0, 1, ..., N-1``.
The inner product between bra and ket is

```math
\begin{aligned}
& \langle \psi_{n \bm{k}} | \psi_{m \bm{k^\prime}} \rangle
= \int_{0}^{N \bm{R}} \psi_{n \bm{k}}^* \psi_{m \bm{k^\prime}} d\bm{r} \\
= & \sum_{l = 0}^{N-1} \exp(i (\bm{k^\prime} - \bm{k}) l \bm{R})
    \int_{0}^{\bm{R}} \psi_{n \bm{k}}^* \psi_{m \bm{k^\prime}} d\bm{r} \\
= & N \delta_{\bm{k},\bm{k^\prime}} \delta_{n m}.
\end{aligned}
```

The second last line uses the fact that the
[discrete Fourier transform](https://en.wikipedia.org/wiki/Discrete_Fourier_transform)
of Dirac delta function is 1,

```math
{\mathcal {F}}(\delta_{n})_{k} =
    \sum_{n=0}^{N-1} \delta_{n} \cdot e^{-{\frac{i2\pi }{N}} kn}
    = 1, k = 0, 1, ..., N-1.
```

so that the delta function can be represented by

```math
\delta_{n} = \frac{1}{N} \sum_{k=0}^{N-1} e^{{\frac{i2\pi }{N}} kn}
```

## Real space representation

In the discrete case, the WFs are the discrete Fourier transform of the Bloch function,

```math
| \bm{R} n \rangle = \frac{1}{N} \sum_{\bm{k}}
    e^{-i \bm{k} \bm{R}} | \psi_{n \bm{k}} \rangle.
```

Note the factor ``\frac{1}{N}`` is multiplied in the forward transform,
while the inverse transform is

```math
| \psi_{n \bm{k}} \rangle = \sum_{\bm{k}} e^{i \bm{k} \bm{R}} | \bm{R} n \rangle.
```

Then the orthonormalization relation of WFs are (note the inner product is
integrated over the super cell ``NV``)

```math
\begin{aligned}
& \langle \bm{R} n | \bm{R^\prime} m \rangle \\
= & \frac{1}{N^2} \sum_{\bm{k},\bm{k^\prime}} e^{i \bm{k} \bm{R} - i \bm{k^\prime} \bm{R^\prime}}
    \langle \psi_{n \bm{k}} | \psi_{m \bm{k^\prime}} \rangle \\
= & \frac{1}{N^2} \sum_{\bm{k},\bm{k^\prime}} e^{i \bm{k} \bm{R} - i \bm{k^\prime} \bm{R^\prime}}
    N \delta_{\bm{k},\bm{k^\prime}} \delta_{n m} \\
= & \frac{1}{N} \delta_{n m} \sum_{\bm{k}} e^{i \bm{k} (\bm{R} - \bm{R^\prime})} \\
= & \frac{1}{N} \delta_{n m} N \delta_{\bm{R}, \bm{R^\prime}}
= \delta_{n m} \delta_{\bm{R}, \bm{R^\prime}}.
\end{aligned}
```

### Implementation details

- The `UNK` files (generated by `QE` `pw2wannier90.x`) are not normalized to 1
  in the unit cell, we need to multiply a factor of ``\frac{1}{\sqrt{N}}``.
  This is done in the `read_realspace_wf` function.

- Then calculating the WF center in real space can be easily achieved
  in `center(rgrid::RGrid, W::AbstractArray)`.

  !!! note

      This is not as accurate as computing WF center in reciprocal space,
      unless you have a dense real space sampling (usually determined by the energy
      cutoff in plane-wave code), and generate the real space WF in the ``NV``-sized
      super cell (the WF lives in the super cell!).

- However, `Wannier90` does not consider this factor when writing `xsf` or `cube` files,
  so we remove the factor in `write_realspace_wf` function, such that

  - the output file is consistent with `Wannier90` output.
  - another advantage is that WF is normalized in one unit cell,
    so the numbers in the `xsf` are not too small.

  !!! note

      Somehow the `Wannier90` output `xsf` or `cube` is on a slightly shifted grid,
      in `write_realspace_wf` we just output the WF on the unshifted grid so that
      the origin of the WF 3D mesh is the same as the lattice origin
      (or its periodic replicas translated by lattice vectors).

[^RMP]: Marzari, N.; Mostofi, A. A.; Yates, J. R.; Souza, I. & Vanderbilt, D.,
    Maximally localized Wannier functions: Theory and applications,
    [Reviews of Modern Physics, 2012, 84, 1419](https://dx.doi.org/10.1103/revmodphys.84.1419)

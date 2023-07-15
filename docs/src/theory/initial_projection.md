# Initial projection for plane-wave DFT codes

Here I write down the formulae for initial projection in plane-wave DFT codes,
e.g., as implemented in `DFTK.jl` or `Quantum ESPRESSO`.

In essence, the initial projection ``A_{m n \mathbf{k}}`` is the inner product
between the Bloch wavefunction ``\psi_{m \mathbf{k}}(\mathbf{r})`` and some
initial guess ``\phi_{n \mathbf{k}}(\mathbf{r})``,

```math
A_{m n \mathbf{k}} = \langle \psi_{m \mathbf{k}} | \phi_{n \mathbf{k}} \rangle.
```

## Bloch wavefunction

The Bloch wavefunction is expanded in plane waves,

```math
\psi_{m \mathbf{k}}(\mathbf{r}) = \frac{1}{\sqrt{V}} \sum_{\mathbf{G}}
  c_{m \mathbf{k}}(\mathbf{G}) e^{i (\mathbf{k} + \mathbf{G}) \mathbf{r}},
```

where ``V`` is the volume of the unit cell, the sum is over the
plane-wave expansions (integer multiples of reciprocal lattice vectors)
``\mathbf{G}``, which are limited by the energy cutoff ``E_{\mathrm{cut}}``,

```math
\frac{\hbar^2}{2 m_e} | \mathbf{k} + \mathbf{G} |^2 \leq E_{\mathrm{cut}},
```

where ``m_e`` is the electron mass.
The coefficients ``c_{m \mathbf{k}}(\mathbf{G})`` are obtained from DFT calculations.

## Localized atomic orbitals

Usually we choose some atomic orbitals as the initial projections, given by
a product of radial function ``R(r)`` times spherical harmonics
``Y_{l m}(\theta, \phi)``,

```math
f(\mathbf{r}) = R(r) Y_{l m}(\theta, \phi),
```

where ``\mathbf{r}`` is the vector and ``r = |\mathbf{r}|`` is its norm;
``\theta`` and ``\phi`` are the polar and azimuthal angles or ``\mathbf{r}``,
respectively; ``l`` and ``m`` are the angular momentum, and
the ``Y_{l m}(\theta, \phi)`` are the real spherical harmonics
(as a comparison, we use the notation ``Y_{l}^{m}(\theta, \phi)`` for
the complex spherical harmonics, see [^SphHarm]).

## Bloch sum

The atomic orbitals do not have the periodicity of the lattice, so we need
to sum over the lattice to obtain the Bloch sum of the localized atomic orbitals,

```math
\phi_{\mathbf{k}}(\mathbf{r}) = \frac{1}{N \sqrt{V}} \sum_{\mathbf{R}}
  e^{i \mathbf{k} \mathbf{R}} f(\mathbf{r} - \mathbf{R}),
```

where ``N`` is the number of unit cells in the supercell (determined by
the k-point mesh), ``V`` is the volume of unit cell,
and ``\mathbf{R}`` labels the unit cell inside the supercell.

## Fourier transform to plane waves

To obtain the inner product between the Bloch wavefunction and the Bloch
sum of the localized atomic orbitals, we need to Fourier transform the
Bloch sum to plane waves, i.e., compute inner product in ``\mathbf{G}``-vector space.

For a wavefunction having a form of radial functions times spherical harmonics,
i.e., ``f(\mathbf{r}) = R(r) Y_{l m}(\theta, \phi)``,
its Fourier transform ``\hat{f}( \mathbf{q} )`` can be computed by

```math
\begin{aligned}
\hat{f}(\mathbf{q})
&= \int_{{\mathbb R}^{3}} e^{-i \mathbf{q} \cdot \mathbf{r}} f(\mathbf{r}) d\mathbf{r} \\
&= 4 \pi Y_{l m}(\mathbf{q}/q) (-i)^{l}
  \int_{{\mathbb R}^+} r^2 R(r) \ j_{l}(q r) dr,
\end{aligned}
```

where ``q = |\mathbf{q}|``, and
``j_{l}(x)`` is the spherical Bessel function of the first kind [^SphBess].

Here we have used the plane wave expansion in spherical harmonics [^PwExpand]
(valid for both complex and real spherical harmonics),

```math
\begin{aligned}
  e^{i \mathbf{q} \cdot \mathbf{r}}
&= 4 \pi \sum_{l=0}^{\infty} \sum_{m=-l}^{l}
    i^{l} j_{l}(q r) Y_{l}^{m}(\mathbf{q}/q) Y_{l}^{m*}(\mathbf{r}/r) \\
&= 4 \pi \sum_{l=0}^{\infty} \sum_{m=-l}^{l}
    i^{l} j_{l}(q r) Y_{l m}(\mathbf{q}/q) Y_{l m}(\mathbf{r}/r)
\end{aligned}
```

For ``e^{-i \mathbf{q} \cdot \mathbf{r}}``, due to the parity of spherical
harmonics,
``Y_{l}^{m}(- \mathbf{r} ) = (-1)^{l} Y_{l}^{m}( \mathbf{r} )``
(also hold true for real spherical harmonics
``Y_{l m}(- \mathbf{r} ) = (-1)^{l} Y_{l m}( \mathbf{r} )``),
the ``i^l`` in the expression of ``e^{i \mathbf{q} \cdot \mathbf{r}}``
becomes ``(-i)^l`` in that of ``e^{-i \mathbf{q} \cdot \mathbf{r}}``.

!!! note

    To really convince you that there should be a minus sign, I created a
    notebook which includes some numerical tests,
    [click here to download](./pluto_fourier.html).

Now, using the orthogonality of spherical harmonics,

```math
\int_{\theta = 0}^{\pi} \int_{\phi=0}^{2\pi}
Y_{l}^{m} Y_{l^{\prime}}^{m^{\prime} *} d\Omega
= \delta_{l l^{\prime}} \delta_{m m^{\prime}},
```

we reach the expression of ``\hat{f}(\mathbf{q})``.
Note in the above equation ``d\Omega = \sin \theta d\theta d\phi`` is the
solid angle element, and ``dV = r^2 dr d\Omega`` is the volume element.

Note that alternatively, we can use Hankel transform to reach the same result,
see [Appendix: Hankel transform](#Appendix:-Hankel-transform).

Thus, ``f(\mathbf{r})`` can be expanded as

```math
f(\mathbf{r}) = \frac{1}{(2\pi)^3} \int d\mathbf{q} \hat{f}(\mathbf{q})
\exp(i \mathbf{q} \cdot \mathbf{r}).
```

Then the Fourier transform of the Bloch sum is

```math
\begin{aligned}
\phi_{\mathbf{k}}(\mathbf{q})
&= \mathcal{F}_{\mathbf{r} \rightarrow \mathbf{q}}
  \{ \phi_{\mathbf{k}} (\mathbf{r} - \mathbf{\tau}) \} \\
&= \mathcal{F} \{ \frac{1}{N \sqrt{V}}
  \sum_{\mathbf{R}} \exp(i \mathbf{k} \mathbf{R})
  f(\mathbf{r} - \mathbf{\tau} - \mathbf{R}) \} \\
&= \frac{1}{N \sqrt{V}} \sum_{\mathbf{R}}
  \exp(i \mathbf{k} \mathbf{R})
  \mathcal{F} \{ f(\mathbf{r} - \mathbf{\tau} - \mathbf{R}) \} \\
&= \frac{1}{N \sqrt{V}} \sum_{\mathbf{R}}
  \exp(i \mathbf{k} \mathbf{R})
  \hat{f}(\mathbf{q}) \exp(-i \mathbf{q} (\mathbf{\tau} + \mathbf{R}) ),
\end{aligned}
```

where ``\mathbf{\tau} = (\tau_x, \tau_y, \tau_z)`` is the center of the orbital.

## Inner product

Now we can compute the initial projection ``A_{m n \mathbf{k}}``.
Note that the Bloch state ``\psi_{n \mathbf{k}}(\mathbf{r})`` are actually
a sum over ``e^{i (\mathbf{k} + \mathbf{G}) \mathbf{r}}``,
so we need to set the ``\mathbf{q}`` in ``\phi_{\mathbf{k}}(\mathbf{q})``
to ``\mathbf{k} + \mathbf{G}`` so that the inner product is on the same basis
functions.

```math
\phi_{\mathbf{k}}(\mathbf{k + G})
= \frac{1}{\sqrt{V}} \exp(-i (\mathbf{k + G}) \mathbf{\tau} )
  \hat{f}(\mathbf{k + G}).
```

Since the localized orbitals ``f(\mathbf{r})`` might not be normalized,
often we apply Löwdin orthonormalization on the ``\mathbf{q}``-space wavefunction
``\phi_{\mathbf{k}}(\mathbf{q})`` at each ``\mathbf{k}``-point.

Finally, the initial projection is

```math
\begin{aligned}
A_{m n \mathbf{k}}
&= \psi_{m \mathbf{k}}^* \phi_{n \mathbf{k}} \\
&= \sum_{\mathbf{G}} c_{m \mathbf{k}}(\mathbf{G})^* \phi_{n \mathbf{k}}(\mathbf{k + G}).
\end{aligned}
```

## Appendix: Hankel transform

The Fourier transform of ``f(\mathbf{r})`` can also be obtained using
Hankel transform [^Hankel],

```math
\begin{aligned}
F_{\nu }(k) &= \int _{0}^{\infty } f(r) J_{\nu }(kr) \,r \, \mathrm{d} r, \\
f(r) &= \int _{0}^{\infty } F_{\nu }(k) J_{\nu }(kr) \,k \, \mathrm{d} k,
\end{aligned}
```

where ``\nu \in \mathbb{C}``, ``J_{\nu }(x)`` is the Bessel function;
for integer ``n``, the spherical Bessel function ``j_{n}(x)``
is related to the Bessel function by

```math
j_{n}(x) = {\sqrt {\frac {\pi }{2x}}} J_{n + {\frac{1}{2}}}(x).
```

By using the relation between Fourier transform and Hankel transform [^Hankel],
we can write the Fourier transform as

```math
\hat{f}(\mathbf{q}) = (2\pi)^{\frac{3}{2}} (-i)^{l} R_{l}(q) Y_{l m}(\mathbf{q}/q),
```

where

```math
\begin{aligned}
R_{l}(q)
&= \frac{1}{\sqrt{q}} \int_{{\mathbb R}^+} \sqrt{r} R(r) J_{l+\frac{1}{2}}(q r) r dr \\
&= \sqrt{\frac{2}{\pi}} \int_{{\mathbb R}^+} r^2 R(r) \ j_{l}(q r) dr,
\end{aligned}
```

restoring the same equation.

## References

[^SphHarm]: <https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form>
[^SphBess]: <https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions>
[^PwExpand]: <https://en.wikipedia.org/wiki/Plane-wave_expansion>
[^Hankel]: <https://en.wikipedia.org/wiki/Hankel_transform>

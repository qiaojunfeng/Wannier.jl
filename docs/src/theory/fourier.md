# Fourier transforms

The Fourier-space frequencies for Wannier interpolation, also called
``\\mathbf{R}``-vectors.

The Wannier interpolation is a Fourier transform of a (theoretically continuous
periodic function (Bloch wavefunction) into a Fourier series.
However, since the wavefunctions are discretized in ``\\mathbf{k}``-space on
a grid, the Fourier series of a coutinuous functionis turns into
a discrete-time Fourier transform (DTFT).
For the usual DTFT, the signal is transformed from time domain to frequency
domain. Here, the time domain is the ``\\mathbf{k}``-space, while the
frequency domain is the ``\\mathbf{R}``-space.

To improve interpolation accuracy, the frequencies (often they are called the
``\\mathbf{R}``-vectors, since they live in 3D space) are not chosen inside
a parallelepiped, but inside a Wigner-Seitz cell to make sure most of the
large-magnitude interactions (e.g., Hamiltonian ``H(\\mathbf{R})``) are
included. Therefore, we have some special algorithms to generate the
``\\mathbf{R}``-space domain.

Th function [`generate_Rspace`] can generate two kinds of R-space domains, by providing an argument
of type [`AbstractRspace`](@ref):
- [`WignerSeitzRspace`](@ref): Wigner-Seitz interpolation,
    returns a [`WignerSeitzRspace`](@ref)
- [`MDRSRspace`](@ref): minimal-distance replica selection interpolation,
    returns a [`MDRSRspace`](@ref)

However, the details of how to place the R-vectors and setting their degeneracies
are irrelevant to the Fourier / inverse Fourier transforms, therefore, we
provide a function [`simplify`](@ref) that can convert each kind of RspaceDomain
to a [`BareRspace`](@ref), enabling faster Fourier transforms.

Depending on the type of RspaceDomain, there are three kinds of Fourier transform:
- `BareRspace`: simple Fourier sum
    ```math
    O_{mn}(\\mathbf{R}) = \\frac{1}{N_{\\mathbf{k}}}
    \\sum_{\\mathbf{k}} \\exp(-i {\\mathbf{k}} \\mathbf{R}) O_{mn}(\\mathbf{k}),
    ```
- `WignerSeitzRspace`: Wigner-Seitz interpolation, the forward Fourier transform
    is the same as `BareRspace`,
- `MDRSRspace`: minimal-distance replica selection interpolation
    ```math
    O_{mn}(\\widetilde{ \\mathbf{R} }) =
    \\sum_{ \\mathbf{R} } \\frac{1}{\\mathcal{N}_{\\mathbf{R}} \\mathcal{N}_{mn \\mathbf{R}}}
    O_{mn}(\\mathbf{R})
    \\sum_{i=1}^{\\mathcal{N}_{mn \\mathbf{R}}}
    \\delta_{ \\widetilde{\\mathbf{R}}, \\mathbf{R} + \\mathbf{T}_{mn \\mathbf{R}}^{(i)} },
    ```

where
- ``N_{\\mathbf{k}}``: the total number of kpoints
- ``\\mathcal{N}_{\\mathbf{R}}``: the degeneracy of R-vectors
- ``\\mathcal{N}_{mn \\mathbf{R}}`` is the degeneracy of ``\\mathbf{T}_{mn \\mathbf{R}}`` vectors


## invfourier
### WS

```math
O_{mn}(\\mathbf{k}) = \\sum_{\\mathbf{R}} \\frac{1}{\\mathcal{N}_{\\mathbf{R}}}
\\exp(i \\mathbf{k} \\mathbf{R}) O_{mn}(\\mathbf{R}),
```
where ``\\mathcal{N}_{\\mathbf{R}}`` is the degeneracy of R vectors (not the total number of R vectors).


### MDRS


```math
X_{mn}(\\mathbf{k}) = \\sum_{\\mathbf{R}}
\\frac{1}{ \\mathcal{N}_{mn \\mathbf{R}} } X_{mn}(\\mathbf{R})
\\sum_{j=1}^{\\mathcal{N}_{ mn \\mathbf{R} }}
\\exp\left( i \\mathbf{k} \cdot \left( \\mathbf{R} + \\mathbf{T}_{ mn \\mathbf{R} }^{(j)} \right) \right)
```
where ``\\mathcal{N}_{\\mathbf{R}}`` is the degeneracy of R vectors (not the total number of R vectors),
and ``\\mathcal{N}_{ mn \\mathbf{R} }`` is the degeneracy of ``\\mathbf{T}_{ mn \\mathbf{R} }`` vectors.


#### MDRSv2

```math
X_{mn}(\\mathbf{k}) = \\sum_{ \\widetilde{\\mathbf{R}} }
\\exp\left( i \\mathbf{k} \\widetilde{ \\mathbf{R} } \right) \\widetilde{X}_{mn}(\\widetilde{\\mathbf{R}})
```

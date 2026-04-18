# # Constraining Wannier function centers

#=
```@meta
CurrentModule = Wannier
```
=#

#=
The matrix manifold based view point of Wannierization provides a transparent
way to apply constraint on WF centers, during both disentanglement and maximal
localization. Instead of localizing the MV spread functional, we optimize

```math
\Omega_r = \Omega + \lambda(\langle \bm{r} \rangle - \bm{r}_0)^{2},
```

where ``\lambda`` is a Lagrange multiplier, and ``r_0`` is the desired center.

In this tutorial, we will start from the `s, p` initial projections of
silicon valence + conduction bands, add WF center penalty, to force the
WFs to be centered at the bond centers, i.e., bonding and anti-bonding orbitals.

## Outline

1. construct a [`Model`](@ref), by reading the `win`, `amn`, `mmn`, and `eig` files
2. localize without WF center penalty — `localize(model)`
3. localize with WF center penalty — `localize(CenteredVariance(r₀, λ), model)`
=#

# ## Preparation
# Load the package
using Wannier
using Wannier.Datasets
using StaticArrays: SVector

#=
## Model construction

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct a [`Model`](@ref) that abstracts the calculation
=#
model = load_dataset("Si2")

#=
## Localization without center penalty

First let's localize the valence + conduction manifold, without any WF
center penalty.
=#
U = localize(model);
# Now we arrive at `s, p` WFs centered at atom centers,
spread(model, U)

#=
## Localization with center penalty

As has been done in the [1. Maximal localization of isolated manifold](@ref)
tutorial, we can use the [`find_nearests`](@ref) function to find the bond
centers. Here we just hard-code the eight targets (Cartesian, Å):
=#
r₀ = [
    SVector( 0.67882, -0.67882, -0.67882),
    SVector(-0.67882, -0.67882,  0.67882),
    SVector(-0.67882,  0.67882, -0.67882),
    SVector( 0.67882,  0.67882,  0.67882),
    SVector( 0.67882, -0.67882, -0.67882),
    SVector(-0.67882, -0.67882,  0.67882),
    SVector(-0.67882,  0.67882, -0.67882),
    SVector( 0.67882,  0.67882,  0.67882),
]

#=
We use `1.0` as the Lagrange multiplier,
=#
λ = 1.0

#=
The fused penalty-aware spread + gradient is exposed through the
[`CenteredVariance`](@ref) objective. Running [`localize`](@ref) on it
dispatches to the same LBFGS driver used for the plain variance case, only
with the WF-center penalty folded into every `fg!` sweep.
=#
U1 = localize(CenteredVariance(r₀, λ), model);

#=
Inspect the final base spread and centers
=#
spread(model, U1)

#=
The WFs are now positioned where we want, giving 4 degenerate bonding WFs
and 4 degenerate anti-bonding WFs. The base spread is typically smaller than
the plain `s, p` disentanglement case because these are
bonding/anti-bonding combinations. 🎉
=#

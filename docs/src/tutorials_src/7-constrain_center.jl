# # 7. Constraining the WF center

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
2. disentangle, without WF center penalty
3. disentangle, with WF center penalty

!!! tip

    This is a HTML version of the tutorial, you can download corresponding
    - Jupyter notebook: [`7-constrain_center.ipynb`](./7-constrain_center.ipynb)
    - Julia script: [`7-constrain_center.jl`](./7-constrain_center.jl)
=#

# ## Preparation
# Load the package
using Wannier

# Path of current tutorial
CUR_DIR = "7-constrain_center"

#=
## Model generation

We will use the [`read_w90`](@ref) function to read the
`win`, `amn`, `mmn`, and `eig` files, and construct a [`Model`](@ref) that abstracts the calculation
=#
model = read_w90("$CUR_DIR/si2")

#=
## Disentanglement

First let's disentangle the valence + conduction manifold, without WF center penalty
=#
A = disentangle(model);
# Now we arrive at `s, p` WFs centered at atom centers,
omega(model, A)

#=
## Disentanglement with constraint

As has been done in the [1. Maximal localization of isolated manifold](@ref) tutorial,
we can use the [`find_nearests`](@ref) function to find the bond centers, which are
=#
r₀ = [
    0.67882 -0.67882 -0.67882 0.67882 0.67882 -0.67882 -0.67882 0.67882
    -0.67882 -0.67882 0.67882 0.67882 -0.67882 -0.67882 0.67882 0.67882
    -0.67882 0.67882 -0.67882 0.67882 -0.67882 0.67882 -0.67882 0.67882
]
#=
note each column is a target center.

We use `1.0` as the Lagrange multiplier,
=#
λ = 1.0

#=
First, we calculate the initial spread with WF center penalty,
now we need to use [`omega_center`](@ref) instead of [`omega`](@ref),
=#
Wannier.omega_center(model, r₀, λ)
#=
where the last three columns, `ω`, `ωc`, `ωt` are the
MV spread, the penalty, and the total spread for each WF.
As expected, there are large penalties on WFs due to there centers.

Then disentangle with [`disentangle_center`](@ref) function,
=#
A1 = disentangle_center(model, r₀, λ);
# the final spreads are
Wannier.omega_center(model, A1, r₀, λ)
#=
Look, the the WF centers are positioned where we want, and now
these are 4 degenerate bonding WFs, and 4 degenerate anti-bonding WFs!

Note now the total spread is smaller than the `s, p` disentanglement case,
since they are the bonding/anti-bonding combinations. 🎉
=#

"""Abstract type for terms in spread functional."""
abstract type AbstractSpreadTerm end

struct MarzariVanderbiltSpread <: AbstractSpreadTerm end

struct CenterPenalty{T} <: AbstractSpreadTerm
    centers::Vector{Vec3{T}}
    weights::Vector{T}
end

function (cache::Cache)(::Val{:Cost}, term::MarzariVanderbiltSpread)
    return sum(Ω -> Ω.Ω, cache.UtMU)
end

function (cache::Cache)(::Val{:Gradient}, term::MarzariVanderbiltSpread)
    return sum(Ω -> Ω.Ω, cache.UtMU)
end

@inline function build_cost_function(
    cache::Cache, terms::AbstractVector{<:AbstractSpreadTerm}
)
    return sum(term -> cache(Val(:Cost), term), terms)
end

@inline function build_gradient(cache::Cache, terms::AbstractVector{<:AbstractSpreadTerm})
    return sum(term -> cache(Val(:Gradient), term), terms)
end

@inline function build_optimization_problem(model::Model)
    cache = Cache(model)
    cost = build_cost_function(cache, model.terms)
    gradient = build_gradient(cache, model.terms)
    if isentanged(model)
        manifold = ProductManifold()
        initial_point = model.U
    else
    end
    return manifold, cost, gradient, initial_point
end

# https://manoptjl.org/stable/tutorials/CountAndCache/#Summary

using Manopt
using Random
using LinearAlgebra
using Manifolds
using LRUCache

m = 25
Random.seed!(42)
A = randn(m + 1, m + 1)
A = Symmetric(A)
p_star = eigvecs(A)[:, end] # minimizer (or similarly -p)
f_star = -eigvals(A)[end] # cost (note that we get - the largest Eigenvalue)

N = Sphere(m);

g(M, p) = -p' * A * p
∇g(p) = -2 * A * p
grad_g(M, p) = project(M, p, ∇g(p))
grad_g!(M, X, p) = project!(M, X, p, ∇g(p))

struct StorageG{T,M}
    A::M
    Ap::T
    p::T
end
function (g::StorageG)(::Val{:Cost}, M::AbstractManifold, p)
    if !(p == g.p) #We are at a new point -> Update
        g.Ap .= g.A * p
        g.p .= p
    end
    return -g.p' * g.Ap
end
function (g::StorageG)(::Val{:Gradient}, M::AbstractManifold, X, p)
    if !(p == g.p) #We are at a new point -> Update
        g.Ap .= g.A * p
        g.p .= p
    end
    X .= -2 .* g.Ap
    project!(M, X, p, X)
    return X
end

p0 = [(1 / 5 .* ones(5))..., zeros(m - 4)...];
#Define the new functor
storage_g = StorageG(A, zero(p0), zero(p0))
# and cost and gradient that use this functor as
g3(M, p) = storage_g(Val(:Cost), M, p)
grad_g3!(M, X, p) = storage_g(Val(:Gradient), M, X, p)
@time s3 = gradient_descent(
    N,
    g3,
    grad_g3!,
    p0;
    stopping_criterion=StopWhenGradientNormLess(1e-5),
    evaluation=InplaceEvaluation(),
    count=[:Cost, :Gradient],
    cache=(:LRU, [:Cost, :Gradient], 2),
    return_objective=true,#, return_state=true
)

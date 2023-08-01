
function localize(model::Model)
    manifold, cost, gradient, initial_point = build_optimization_problem(model)
    return gradient_descent(cost, gradient, manifold, initial_point)
end

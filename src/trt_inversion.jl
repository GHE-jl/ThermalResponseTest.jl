using Optimization, OptimizationOptimJL, ForwardDiff, SpecialFunctions, DataFrames

"""
This script uses the direct models defined in GroundHeatExchanger.jl.
"""

"""
    fit_ils(t, T, q, r_b, T₀; lb=[0.5, 0.01], ub=[5.0, 0.5])
"""
function fit_ils(t, T, q, rb, T₀; lb=[0.5, 0.01], ub=[5.0, 0.5])
    # Average heat rate per meter (assuming you know borehole depth L)
    q_avg = sum(q) / length(q)
    
    # Define Loss Function (Mean Squared Error)
    function loss(p, _)
        pred = ils(t, rb, p, Cs)
        return sum((pred .- T).^2)
    end

    # Set up Optimization.jl
    opt_func = OptimizationFunction(loss, Optimization.AutoForwardDiff())
    p_initial = [2.0, 0.1] # Initial guesses for [λ, Rb]
    
    prob = OptimizationProblem(opt_func, p_initial, lb=lb, ub=ub)
    
    # Solve using L-BFGS
    sol = solve(prob, LBFGS())
    
    return (λ = sol.u[1], Rb = sol.u[2], original_sol = sol)
end
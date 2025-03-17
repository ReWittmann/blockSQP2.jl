using blockSQP
using Test
using Optimization
using ForwardDiff
using Test

@testset "Optimization.jl " begin

    rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
    x0 = zeros(2)
    _p = [1.0, 1.0]

    @testset "With constraints" begin
        cons(res, x, p) = (res .= [x[1]^2 + x[2]^2, x[1] * x[2]])
        optprob_wcons = OptimizationFunction(rosenbrock, Optimization.AutoForwardDiff(), cons = cons)
        prob_wcons = OptimizationProblem(optprob_wcons, x0, _p, lcons = [-Inf, -1.0], ucons = [0.8, 2.0])
        sol = solve(prob_wcons, BlockSQPOpt())
        @test isapprox(sol.u, [0.751964453; 0.484303067])
    end

    @testset "Without constraints" begin
        optprob_wocons = OptimizationFunction(rosenbrock, Optimization.AutoForwardDiff())
        prob_wocons = OptimizationProblem(optprob_wocons, x0, _p)
        sol = solve(prob_wocons, BlockSQPOpt())
        @test isapprox(sol.u, ones(2))
    end

end

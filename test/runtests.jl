using blockSQP
using Test
using Optimization
using ForwardDiff
using OptimizationMOI, Ipopt
using Test

@testset "Optimization.jl " begin

    @testset "Rosenbrock" begin
        rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
        x0 = zeros(2)
        _p = [1.0, 1.0]

        cons(res, x, p) = (res .= [x[1]^2 + x[2]^2, x[1] * x[2]])
        optprob_wcons = OptimizationFunction(rosenbrock, Optimization.AutoForwardDiff(), cons = cons)
        prob_wcons = OptimizationProblem(optprob_wcons, x0, _p, lcons = [-Inf, -1.0], ucons = [0.8, 2.0])
        sol_bsqp_wcons = solve(prob_wcons, BlockSQPOpt())
        sol_ipopt_wcons = solve(prob_wcons, Ipopt.Optimizer())

        optprob_wocons = OptimizationFunction(rosenbrock, Optimization.AutoForwardDiff())
        prob_wocons = OptimizationProblem(optprob_wocons, x0, _p)
        sol_bsqp_wocons = solve(prob_wocons, BlockSQPOpt())
        sol_ipopt_wocons = solve(prob_wocons, Ipopt.Optimizer())
        @testset "Primal solution" begin
            @test isapprox(sol_bsqp_wcons.u, sol_ipopt_wcons.u)
            @test isapprox(sol_bsqp_wocons.u, sol_ipopt_wocons.u)
        end

        @testset "Lagrange multiplier" begin
            @test isapprox(sol_bsqp_wcons.original.multiplier,
                            sol_ipopt_wcons.original.inner.mult_g, atol=1e-8)
        end
    end
    @testset "LP on unit circle" begin
        lin_ex(x,p) = sum(x)
        function cons_unit(res,x,p)
            res .= sum(x.^2)
        end
        x0 = ones(2)

        optprob_lin = OptimizationFunction(lin_ex, Optimization.AutoForwardDiff(), cons=cons_unit)
        prob_lin = OptimizationProblem(optprob_lin, x0, [], lcons=[1.0], ucons=[1.0])

        sol_bsqp = solve(prob_lin, BlockSQPOpt())
        sol_ipopt = solve(prob_lin, Ipopt.Optimizer())

        @testset "Primal solution" begin
            @test isapprox(sol_bsqp.u, sol_ipopt.u)
        end
        @testset "Lagrange multiplier" begin
            @test isapprox(sol_bsqp.original.multiplier,
                        sol_ipopt.original.inner.mult_g, atol=1e-8)
        end

    end
end

using blockSQP
using Optimization
using Test
using SafeTestsets
using ForwardDiff
using OptimizationMOI, Ipopt
using LinearAlgebra
@testset "blockSQP.jl " begin
    @safetestset "Code quality (Aqua.jl)" begin
        using Aqua
        using blockSQP
        Aqua.test_all(blockSQP)
    end
    @safetestset "Rosenbrock" begin
        using blockSQP
        using Optimization
        using Test
        using ForwardDiff
        using OptimizationMOI, Ipopt
        using Symbolics
        rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
        x0 = zeros(2)
        _p = [1.0, 1.0]

        cons(res, x, p) = (res .= [x[1]^2 + x[2]^2, x[1] * x[2]])
        optprob_wcons = OptimizationFunction(rosenbrock, Optimization.AutoForwardDiff(), cons = cons)
        prob_wcons = OptimizationProblem(optprob_wcons, x0, _p, lcons = [-Inf, -1.0], ucons = [0.8, 2.0])
        sol_bsqp_wcons = solve(prob_wcons, blockSQP.Optimizer())
        sol_ipopt_wcons = solve(prob_wcons, Ipopt.Optimizer())

        optprob_wocons = OptimizationFunction(rosenbrock, Optimization.AutoForwardDiff())
        prob_wocons = OptimizationProblem(optprob_wocons, x0, _p)
        sol_bsqp_wocons = solve(prob_wocons, blockSQP.Optimizer())
        sol_ipopt_wocons = solve(prob_wocons, Ipopt.Optimizer())

        @testset "Successful optimization" begin
            @test sol_bsqp_wcons.retcode == ReturnCode.Success
            @test sol_ipopt_wcons.retcode == ReturnCode.Success
            @test sol_bsqp_wocons.retcode == ReturnCode.Success
            @test sol_ipopt_wocons.retcode == ReturnCode.Success
        end
        @testset "Primal solution" begin
            @test isapprox(sol_bsqp_wcons.u, sol_ipopt_wcons.u; atol = 1e-5)
            @test isapprox(sol_bsqp_wocons.u, sol_ipopt_wocons.u; atol = 1e-5)
        end

        @testset "Lagrange multiplier" begin
            @test isapprox(sol_bsqp_wcons.original.multiplier,
                            sol_ipopt_wcons.original.inner.mult_g; atol=1e-5)
        end
        
        @testset "Cache" begin
            cache = Optimization.init(prob_wocons, blockSQP.Optimizer())
            sol = Optimization.solve!(cache)
            @test sol.retcode == ReturnCode.Success
        end

        @testset "Callback" begin
            cb(x,l) = begin
                println(l)
                return false
            end

            cb1(x,l) = begin
                println(l)
                return x.iter > 5 # Stop after five iterations
            end

            sol_bsqp_wcons_cb = solve(prob_wcons, blockSQP.Optimizer(); callback=cb)
            sol_bsqp_wcons_cb1 = solve(prob_wcons, blockSQP.Optimizer(); callback=cb1)

            @test sol_bsqp_wcons_cb.retcode == ReturnCode.Success
            @test sol_bsqp_wcons_cb1.retcode == ReturnCode.Default
        end
    end

    @safetestset "Sparse problems" begin
        using blockSQP
        using Optimization
        using Test
        using ForwardDiff
        using OptimizationMOI, Ipopt
        _f(x,p) = x[1]^2 - 0.5*x[2]^2
        _g(res, x,p) = res .=  [x[1] - x[2]]

        optprob_wcons = OptimizationFunction(_f, Optimization.AutoForwardDiff(), cons = _g)

        prob = OptimizationProblem(optprob_wcons, [10.0, 10.0], Float64[], lcons = [0.0], ucons = [0.0])
        using Symbolics
        blockIdx_calc = blockSQP.compute_hessian_blocks(prob)
        sol_sparse_1 = solve(prob, blockSQP.Optimizer(); blockIdx=blockIdx_calc)
        sol_sparse_2 = solve(prob, blockSQP.Optimizer(); blockIdx=[0,1,2])
        options = blockSQP.Options(sparse=true, hess_approx=:SR1)
        sol_sparse_3 = solve(prob, blockSQP.Optimizer(); options=options)
        @test SciMLBase.successful_retcode(sol_sparse_1) && SciMLBase.successful_retcode(sol_sparse_2) &&
                SciMLBase.successful_retcode(sol_sparse_3)
        @test isapprox(sol_sparse_1.u, sol_sparse_2.u; atol = 1e-5)
        @test isapprox(sol_sparse_2.u, sol_sparse_3.u; atol = 1e-5)
    end
    @safetestset "LP on unit circle" begin
        using blockSQP
        using Optimization
        using Test
        using ForwardDiff
        using OptimizationMOI, Ipopt
        
        lin_ex(x,p) = sum(x)
        function cons_unit(res,x,p)
            res .= sum(x.^2)
        end
        x0 = ones(2)

        optprob_lin = OptimizationFunction(lin_ex, Optimization.AutoForwardDiff(), cons=cons_unit)
        prob_lin = OptimizationProblem(optprob_lin, x0, Float64[], lcons=[1.0], ucons=[1.0])

        sol_bsqp = solve(prob_lin, blockSQP.Optimizer())
        sol_ipopt = solve(prob_lin, Ipopt.Optimizer())

        @testset "Primal solution" begin
            @test isapprox(sol_bsqp.u, sol_ipopt.u; atol = 1e-5)
        end
        @testset "Lagrange multiplier" begin
            @test isapprox(sol_bsqp.original.multiplier,
                        sol_ipopt.original.inner.mult_g; atol = 1e-5)
        end

    end
    @safetestset "compute_hessian_blocks" begin
        using blockSQP
        using Optimization
        using Symbolics
        using Test
        using LinearAlgebra
        H1 = [1 0 ; 0 1]
        blocks_h1 = blockSQP.compute_hessian_blocks(H1)
        @test blocks_h1 == [0,1,2]

        block1, block2 = diagm(ones(3)), ones(2,2)
        H2 = [block1 zeros(3,2); zeros(2,3) block2]
        blocks_h2 = blockSQP.compute_hessian_blocks(H2)
        @test blocks_h2 == [0,1,2,3,5]

        H3 = ones(4,4)
        blocks_h3 = blockSQP.compute_hessian_blocks(H3)
        @test blocks_h3 == [0,4]
    end
    @safetestset "Lotka Volterra Fishing" begin
        include("examples/test_lotka.jl")
    end
end
#include("../src/blockSQP.jl")
#using .blockSQP
using blockSQP

###Example problem taken from blockSQP paper (Janka 2016) (sparse version)###

#min  x1^2 - 0.5x2^2
#       s.t. 0 <= x1 - x2 <= 0
#           -Inf < x1, x2 < Inf

#Initial point: x1 = 10, x2 = 10, lambda = [0., 0., 0.]


f = x::Array{Float64, 1} -> x[1]^2 - 0.5*x[2]^2
g = x::Array{Float64, 1} -> Float64[x[1] - x[2]]
grad_f = x::Array{Float64, 1} -> Float64[2*x[1], -x[2]]
jac_g = x::Array{Float64, 1} -> Float64[1 -1]
nnz = 2
jac_g_nz = x::Array{Float64, 1} -> Float64[1, -1]
jac_g_row = Int32[0, 0]
jac_g_colind = Int32[0, 1, 2]

#Set bounds
lb_var = Float64[-Inf, -Inf]
ub_var = Float64[Inf, Inf]
lb_con = Float64[0.0]
ub_con = Float64[0.0]

#Set initial values
x0 = Float64[10.0, 10.0]
lambda0 = Float64[0., 0., 0.]

prob = blockSQP.blockSQPProblem(f,g, grad_f, jac_g,
                            lb_var, ub_var, lb_con, ub_con,
                            x0, lambda0, blockIdx = Int32[0, 1, 2])
blockSQP.make_sparse!(prob, Int32(nnz), jac_g_nz, jac_g_row, jac_g_colind)


QPopts = qpOASES_options(sparsityLevel = 2)
opts = blockSQPOptions(
                       maxiters = 100,
                       opt_tol = 1.0e-12,
                       feas_tol = 1.0e-12,
                       enable_linesearch = false,
                       hess_approx = :SR1,
                       fallback_approx = :BFGS,
                       print_level = 2,
                       indef_delay = 1
)


stats = blockSQP.SQPstats("./")

meth = blockSQP.Solver(prob, opts, stats)

blockSQP.init!(meth)
ret = blockSQP.run!(meth, Int32(100), Int32(1))
blockSQP.finish!(meth)

x_opt = blockSQP.get_primal_solution(meth)
lam_opt = blockSQP.get_dual_solution(meth)

print("Primal solution\n", x_opt, "\nDual solution\n", lam_opt, "\n")

#include("../src/blockSQP.jl")
#using .blockSQP
using blockSQP

#Set objective, constraints and their first derivatives
f = x::Array{Float64, 1} -> x[1]^2 - 0.5*x[2]^2
g = x::Array{Float64, 1} -> Float64[x[1] - x[2]]
grad_f = x::Array{Float64, 1} -> Float64[2*x[1], -x[2]]
jac_g = x::Array{Float64, 1} -> Float64[1 -1]

#Set bounds
lb_var = Float64[-Inf, -Inf]
ub_var = Float64[Inf, Inf]
lb_con = Float64[0.0]
ub_con = Float64[0.0]


#Set initial values
x0 = Float64[10.0, 10.0]
lambda0 = Float64[0., 0., 0.]

prob = blockSQP.Problem(f,g, grad_f, jac_g,
                        lb_var, ub_var, lb_con, ub_con,
                        x0, lambda0, blockIdx = Int32[0, 1, 2])



QPopts = qpOASES_options(sparsityLevel = 0,
                         printLevel = 0,
                         terminationTolerance=1.0e-10
                         )                            
                            
#Set options
opts = blockSQP.Options(opt_tol = 1e-12,
                       feas_tol = 1e-12,
                       enable_linesearch = false,
                       hess_approx = :BFGS,         #Either :BFGS, "BFGS" or map(c->Cchar(c),collect("BFGS"))
                       lim_mem = true,
                       mem_size = 20,
                       sizing = :None,
                       fallback_approx = :BFGS,
                       fallback_sizing = :None,
                       sparse = false,
                       print_level = 2,
                       qpsol = "qpOASES",
                       qpsol_options = QPopts,
                       indef_delay = 1
                       )

stats = blockSQP.SQPstats("./")

meth = blockSQP.Solver(prob, opts, stats)
blockSQP.init!(meth)

ret = blockSQP.run!(meth, 100, 1)

blockSQP.finish!(meth)

x_opt = blockSQP.get_primal_solution(meth)
lam_opt = blockSQP.get_dual_solution(meth)

print("Primal solution\n", x_opt, "\nDual solution\n", lam_opt, "\n")
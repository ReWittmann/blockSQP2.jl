using blockSQP
using Optimization
#Set objective, constraints and their first derivatives
f(x) =  x[1]^2 - 0.5*x[2]^2
g(x) = [x[1] - x[2]]
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



prob = blockSQP.blockSQPProblem(f,g, grad_f, jac_g,
                            lb_var, ub_var, lb_con, ub_con,
                            x0, lambda0; blockIdx = Int32[0,1,2], nnz=2,
                            jac_g_row = [0, 0], jac_g_colind=[0,1,2], jac_g_nz= (x) -> [1.0, -1.0])



#Set options
BlockSQPOptions().blockHess
BlockSQPOptions().hessUpdate
BlockSQPOptions().printLevel

opts = BlockSQPOptions(opttol=1e-12, sparseQP=2, hessUpdate=1)

stats = blockSQP.SQPstats("./")

meth = blockSQP.Solver(prob, opts, stats)
blockSQP.init(meth)

ret = blockSQP.run(meth, 100, 1);

blockSQP.finish(meth)

x_opt = blockSQP.get_primal_solution(meth)
lam_opt = blockSQP.get_dual_solution(meth)
print("Primal solution\n", x_opt, "\nDual solution\n", lam_opt, "\n")



using Optimization, ForwardDiff

_f(x,p) = x[1]^2 - 0.5*x[2]^2
_g(res, x,p) = res .=  [x[1] - x[2]]

optprob_wcons = OptimizationFunction(_f, Optimization.AutoForwardDiff(), cons = _g)

prob = OptimizationProblem(optprob_wcons, x0, [], lcons = [0.0], ucons = [0.0])
sol_bsqp_wcons = solve(prob, BlockSQPOpt(); sparsity=true)

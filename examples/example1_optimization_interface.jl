
using blockSQP, Optimization, ForwardDiff

_f(x,p) = x[1]^2 - 0.5*x[2]^2
_g(res, x,p) = res .=  [x[1] - x[2]]
x0 = Float64[10.0, 10.0]

optprob_wcons = OptimizationFunction(_f, Optimization.AutoForwardDiff(), cons = _g)
prob = OptimizationProblem(optprob_wcons, x0, [], lcons = [0.0], ucons = [0.0])
sol_bsqp_wcons = solve(prob, BlockSQPOpt(); sparsity=true)

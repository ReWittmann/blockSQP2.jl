using blockSQP
using Optimization, OptimizationMOI
using ModelingToolkit
using Ipopt

rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)

cons(res, x, p) = (res .= [x[1]^2 + x[2]^2, x[1] * x[2]])
optprob_wcons = OptimizationFunction(rosenbrock, Optimization.AutoModelingToolkit(), cons = cons)
prob_wcons = OptimizationProblem(optprob_wcons, x0, [1.0,1.0]; lcons = [-Inf, -1.0], ucons = [0.8, 2.0])
sol = solve(prob_wcons, BlockSQPOpt(); max_conv_QPs = 1)
#sol = solve(prob_wcons, Ipopt.Optimizer())
print(sol)

using NLPstructures

const lotka_params = Dict{Symbol, Float64}(
    :c0 => 0.4,
    :c1 => 0.2,
    :t0 => 0.,
    :t1 => 12.
)

const lotka_vparams = Dict{Symbol, Vector{Float64}}(
    :x_init => [0.5, 0.7]
)

function lotka_rhs(x,u)
    c0, c1 = (lotka_params[x] for x in (:c0, :c1))
    return [x[1] - x[1]*x[2] - u[1]*c0*x[1], -x[2] + x[1]*x[2] - u[1]*c1*x[2]]
end
lotka_nx = 2
lotka_nu = 1

function lotka_quad(x,u)
    return [(x[1] - 1)^2 + (x[2] - 1)^2]
end
lotka_nu = 1

function RK4_step(rhs, quad, x0, u, DT = 1.0)
    K1, K1q = rhs(x0, u), quad(x0, u)
    K2, K2q = rhs(x0 + DT/2 * K1, u), quad(x0 + DT/2 * K1, u)
    K3, K3q = rhs(x0 + DT/2 * K2, u), quad(x0 + DT/2 * K2, u)
    K4, K4q = rhs(x0 + DT * K3, u), quad(x0 + DT * K3, u)
    x1 = x0 + DT/6*(K1 + 2*K2 + 2*K3 + K4)
    Q  = DT/6*(K1q + 2*K2q + 2*K3q + K4q)
    return x1, Q
end

function RK4_M(rhs, quad, x0, u, DT = 1.0, M = 2)
    x, q = RK4_step(rhs, quad, x0, u, DT/M)
    for i = 2:M
        x, _q = RK4_step(rhs, quad, x, u, DT/M)
        q += _q
    end
    return x, q
end

function ODEsol_multi(rhs, quad, xk::AbstractArray{DTP}, uk::AbstractArray{DTP}, T = 12.0, M = 2) where DTP <: Number
    _N_stages = size(xk, 2)
    
    DT = T/_N_stages
    # TMP = map(1:_N_stages) do i
    #     RK4_M(rhs, quad, xk[:,i], uk[:,i], DT, M)
    # end
    
    TMP = Vector{Tuple{Vector{DTP}, Vector{DTP}}}(undef, _N_stages)
    Threads.@threads for i in 1:_N_stages
        TMP[i] = RK4_M(rhs, quad, xk[:,i], uk[:,i], DT, M)
    end
    
    q_out = sum(last.(TMP))
    x_out = reduce(hcat, first.(TMP))
    return x_out, q_out
end

function ODEsol_single(rhs, quad, x0::AbstractArray{DTP}, u::AbstractArray{DTP}, T = 12.0, M = 2, arg_N = 100) where DTP <: Number
    DT = T/arg_N
    traj = Vector{DTP}[x0]
    Q = zeros(quad(x0,u))
    xk = x0
    for i = 1:arg_N
        xk, q = RK4_M(rhs, quad, xk, u[:,i], DT, M)
        push!(traj, xk)
        Q += q
    end
    return reduce(hcat, traj), Q
end


N = 100
nx = 2
nu = 1
nq = 1
_subscript(i::Integer) = (i |> digits |> reverse .|> dgt->Char(0x2080+dgt)) |> join

vLayoutPre = TupleBD[]
cLayoutPre = TupleBD[]

states = BlockDescriptor[]
controls = BlockDescriptor[]

matchings = BlockDescriptor{nlpStateMatchings}(tag = :matchings)
MSsys = BlockDescriptor{nlpMultipleShootingSystemSC}(matchings = matchings, tag = :MSsys)

h0 = BlockDescriptor{nlpHess}(parent = MSsys, tag = :h₀)
u0 = BlockDescriptor{nlpVariables}(parent = h0, tag = :u₀)

push!(vLayoutPre, (h0, [(u0, nu)]))
push!(controls, u0)

c1 = BlockDescriptor{nlpMatching}(tag = :m₁, parent = matchings, input = [u0])
push!(cLayoutPre, (c1, nx))

hk, uk, ck = h0, u0, c1
xk = BlockDescriptor()
for i = 1:(N-1)
    global hk = BlockDescriptor{nlpHess}(parent = MSsys, tag = Symbol(:h, _subscript(i)))
    global xk = BlockDescriptor{nlpVariables}(parent = hk, matching = ck, tag = Symbol(:x, _subscript(i)))
    global uk = BlockDescriptor{nlpVariables}(parent = hk, tag = Symbol(:u, _subscript(i)))
    
    push!(vLayoutPre, (hk, [(xk, nx), (uk, nu)]))
    
    push!(states, xk)
    push!(controls, uk)
    
    global ck = BlockDescriptor{nlpMatching}(parent = matchings, input = [xk, uk], tag = Symbol(:m, _subscript(i+1)))
    push!(cLayoutPre, (ck, nx))
end

hN = BlockDescriptor{nlpHess}(parent = MSsys, tag = Symbol(:h, _subscript(N)))
xN = BlockDescriptor{nlpVariables}(parent = hN, matching = ck, tag = Symbol(:x, _subscript(N)))
uN = BlockDescriptor{nlpVariables}(parent = hN, tag = Symbol(:u, _subscript(N)))
push!(vLayoutPre, (hN, [(xN, nx), (uN, 0)]))
push!(states, xN)

vLayoutPre = [(MSsys, vLayoutPre)]
cLayoutPre = [(matchings, cLayoutPre)]

vLayout = to_Axis(vLayoutPre)
cLayout = to_Axis(cLayoutPre)
vBlockD = get_BlockDescriptors(vLayoutPre)
cBlockD = get_BlockDescriptors(cLayoutPre)

NLPstruc = NLPstructure((vBlockD...,), vLayout, (cBlockD...,), cLayout)

nVar = axlength(vLayout)


lotka_MS(x,u) = ODEsol_multi(lotka_rhs, lotka_quad, x, u, lotka_params[:t1] - lotka_params[:t0], 2)


stateidx = splat(vcat)(states .|> st->collect(axsubrange(vLayout, st)))
controlidx = splat(vcat)(controls .|> ctrl->collect(axsubrange(vLayout, ctrl)))


function objective(arg_x)
    let _sidx = stateidx, _cidx = controlidx
        x0 = lotka_vparams[:x_init]
        x = hcat(x0, reshape(arg_x[_sidx], 2, :))
        u = reshape(arg_x[_cidx], 1, :)
        
        _, Q = lotka_MS(x[:,1:end-1], u)
        return first(Q)
    end
end

function shooting_constraints(arg_x)
    let _sidx = stateidx, _cidx = controlidx
        x0 = lotka_vparams[:x_init]
        x = hcat(x0, reshape(arg_x[_sidx], 2, :))
        u = reshape(arg_x[_cidx], 1, :)
        
        S, _ = lotka_MS(x[:, 1:end-1], u)
        return collect(Iterators.flatten(x[:,2:end] - S))
    end
end

x_init = lotka_vparams[:x_init]
x_start = ComponentArray(zeros(axlength(vLayout)), vLayout)
u_start = zeros(1,N)

lb_var = ComponentArray(repeat([0.], axlength(vLayout)), vLayout)
ub_var = ComponentArray(repeat([Inf], axlength(vLayout)), vLayout)
for state in states
    view(x_start, state)[:] = x_init
end
for control in controls
    view(ub_var, control)[:] .= 1
end

lb_con, ub_con = zeros(axlength(cLayout)), zeros(axlength(cLayout))

using blockSQP
using ForwardDiff

f = objective
g = shooting_constraints
grad_f = x->ForwardDiff.gradient(f, x)


using SparseConnectivityTracer
using SparseMatrixColorings
using DifferentiationInterface
using SparseArrays

sparse_forward_backend = AutoSparse(
    AutoForwardDiff();
    sparsity_detector = TracerSparsityDetector(),
    coloring_algorithm = GreedyColoringAlgorithm()
)


using Optimization
using OptimizationMOI
using Ipopt

Opt_objective(u, ::Any) = objective(u)
Opt_cons(res, u, ::Any) = (res[:] .= shooting_constraints(u))

Opt_f = OptimizationFunction(
    Opt_objective, AutoForwardDiff(); cons = Opt_cons
)
Opt_prob = OptimizationProblem(
    Opt_f, collect(x_start), SciMLBase.NullParameters(); lb = collect(lb_var), ub = collect(ub_var), lcons = lb_con, ucons = ub_con
)

print("Note: We are using dense Jacobians for Ipopt, so the runtime comparison will be off.\n")
Ipoptsol = solve(Opt_prob, Ipopt.Optimizer(),
     tol = 1e-6,
     constr_viol_tol = 1e-6,
     hessian_approximation = "limited-memory",
     max_iter = 300, # 165
)


jac_g0 = jacobian(g, sparse_forward_backend, x_start)
ROW = jac_g0.rowval .-1
COLIND = jac_g0.colptr .-1

jac_gNZ(x) = jacobian(g, sparse_forward_backend, x).nzval

("BlockSQP allows passing sparse Jacobians, so it will have a runtime advantage.\n")
prob = blockSQP.blockSQPProblem(
    f, g, grad_f, blockSQP.fnothing,
    collect(lb_var), collect(ub_var), lb_con, ub_con,
    collect(x_start), zeros(NLPstructures.axlength(vLayout) + NLPstructures.axlength(cLayout));
    blockIdx = hessBlockZeroBasedIndex(NLPstruc), jac_g_row = ROW, jac_g_colind = COLIND, jac_g_nz = jac_gNZ,
    nnz = length(ROW), vblocks = blockSQP.create_vBlocks(NLPstruc)
)
opts = blockSQP.sparse_options()
opts.max_extra_steps = 0
opts.automatic_scaling = true
opts.max_conv_QPs = 4
opts.conv_strategy = 2
stats = blockSQP.SQPstats("./")

meth = blockSQP.Solver(prob, opts, stats)

blockSQP.init!(meth)
blockSQP.run!(meth, 200, 0)
blockSQP.finish!(meth)

blockSQP.get_primal_solution(meth)
x_opt = blockSQP.get_primal_solution(meth)
x_opt = ComponentArray(x_opt, NLPstruc.vLayout)

using CairoMakie

opt_x = hcat(x_init, (states .|> st -> x_opt[st]) |> splat(hcat))
opt_u = (controls .|> ctrl -> x_opt[ctrl]) |> splat(hcat)

fig, _ = Figure(), nothing
ax = CairoMakie.Axis(fig[1,1])
Tgrid = range(0,12,N+1)
lines!(Tgrid, opt_x[1,:], label = "x₁")
lines!(Tgrid, opt_x[2,:], label = "x₂")
stairs!(Tgrid[1:end-1], opt_u[1,:], label = "u", color = :red3)
fig[1, 2] = Legend(fig, ax, framevisible = false)
display(fig)

_ = nothing #Any better way to prevent the Julia REPL from printing the last returned value?





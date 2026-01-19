using NLPstructures

const lotka_params = Dict{Symbol, Float64}(
    :c0 => 0.4,
    :c1 => 0.2,
    :x0_init => 0.5,
    :x1_init => 0.7,
    :t0 => 0.,
    :t1 => 12.
)

function lotka_rhs(x,u)
    c0, c1 = (lotka_params[x] for x in (:c0, :c1))
    return [x[1] - x[1]*x[2] - u*c0*x[1], -x[2] + x[1]*x[2] - u*c1*x[2]]
end

function lotka_quad(x,u)
    return (x[1] - 1)^2 + (x[2] - 1)^2
end

function RK4_step(rhs, quad, x0, u, DT = 1.0)
    K1, K1q = rhs(x0, u), quad(x0, u)
    K2, K2q = rhs(x0 + DT/2 * K1, u), quad(x0 + DT/2 * K1, u)
    K3, K3q = rhs(x0 + DT/2 * K2, u), quad(x0 + DT/2 * K2, u)
    K4, K4q = rhs(x0 + DT * K3, U), quad(x0 + DT * K3, u)
    x1=x0+DT/6*(K1 + 2*K2 + 2*K3 + K4)
    Q = DT/6*(K1q + 2*K2q + 2*K3q + K4q)
    return x1, Q
end

function RK4_M(rhs, quad, x0, u, DT = 1.0, M = 2)
    Q = 0
    x = x0
    for i = 1:M
        x, q = RK4_step(rhs, quad, x, u, DT)
        Q += q
    end
    return x, Q
end

function ODEsol_multi(rhs, quad, xk, uk, T = 12.0, M = 2)
    _N_stages = size(xk, 2)
    xkp = repeat(Vector{Float64}[zeros(size(xk, 1))], _N_stages)
    qk = repeat(Vector{Float64}[zeros(size(qk, 1)), _N_stages])
    
    DT = T/_N_stages
    # Threads.@threads 
    for i = 1:_N_stages
        xkp[i][:], qk[i][:] = RK4_M(rhs, quad, xk[:,i], uk[:,i], DT, M)
    end
    x_out = reduce(hcat, xkp)
    q_out = reduce(+, qk)
    return x_out, q_out
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

matchings = BlockDescriptor{nlpMultipleShootingMatchings}(tag = :matchings)
MSsys = BlockDescriptor{nlpMultipleShootingSystemSC}(matchings = matchings, tag = :MSsys)

h0 = BlockDescriptor{nlpHess}(parent = MSsys, tag = :h0)
u0 = BlockDescriptor{nlpVariables}(parent = h0, tag = :u0)

push!(vLayoutPre, (h0, [(u0, nu)]))

ck = BlockDescriptor{nlpMatching}(tag = :m1, parent = matchings, input = [u0])
hk, uk = h0, u0
xk = BlockDescriptor()
for i = 1:(100-1)
    push!(cLayoutPre, ck)
    
    hk = BlockDescriptor{nlpHess}(parent = MSsys, tag = Symbol(:h, _subscript(i)))
    xk = BlockDescriptor{nlpVariables}(parent = hk, matching = ck, tag = Symbol(:x, _subscript(i)))
    uk = BlockDescriptor{nlpVariables}(parent = hk, tag = Symbol(:u, _subscript(i)))
    
    push!(vLayoutPre, (hk, [(xk, nx), (uk, nu)]))
    push!(cLayoutPre, (ck, nx))
    
    push!(states, xk)
    push!(controls, uk)
    
    ck = BlockDescriptor{nlpMatching}(parent = matchings, input = [xk, uk], tag = Symbol(:m, _subscript(i+1)))
end

hN = BlockDescriptor{nlpHess}(parent = MSsys)
xN = BlockDescriptor{nlpVariables}(parent = hN, matching = ck)
push!(vLayoutPre, (hN, [(xN, nx)]))

vLayoutPre = [(MSsys, vLayoutPre)]
cLayoutPre = [(matchings, cLayoutPre)]

vLayout = to_Axis(vLayoutPre)
cLayout = to_Axis(cLayoutPre)
vBlockD = get_BlockDescriptors(vLayoutPre)
cBlockD = get_BlockDescriptors(cLayoutPre)

NLPstruc = NLPstructure((vBlockD...,), vLayout, (cBlockD...,), cLayout)

nVar = axlength(vLayout)


lotka_MS(x,u) = ODEsol_multi(lotka_rhs, lotka_quad, x, u, lotka_params[:t1] - lotka_params[:t0], 2)

function objective(x)
    let ax = vLayout, nVar = nVar, st = states, ct = controls
        cArr = ComponentArray(x, ax)
        x0 = [lotka_params[:x0], lotka_params[:x1]]
        x = vcat([x0], (states .|> x->cArr[x]))
        u = reduce(vcat, (controls .|> x->cArr[x]))
        
        _, Q = lotka_MS(x[:,1:end-1], u)
        return Q
    end
end

function shooting_constraints(x)
    let ax = vLayout, nVar = nVar, st = states, ct = controls
        cArr = ComponentArray(x, ax)
        x0 = [lotka_params[:x0], lotka_params[:x1]]
        x = vcat([x0], (states .|> x->cArr[x]))
        u = reduce(vcat, (controls .|> x->cArr[x]))
        
        S, _ = lotka_MS(x[:, 1:end-1], u[:, 1:end-1])
        return collect(Iterators.flatten(x[:,2:end] - S))
    end
end

x_init = [lotka_params[:x0], lotka_params[:x1]]
x0 = ComponentArray(zeros(axlength(vLayout)), vLayout)
lb_var = ComponentArray(repeat([Inf], axlength(vLayout)), vLayout)
ub_var = ComponentArray(repeat([-Inf], axlength(vLayout)), vLayout)
for state in states
    view(x0, state)[:] = x_init
end
for control in controls
    view(lb_var, control) .= 0
    view(ub_var, control) .= 1
end

lb_con, ub_con = zeros(axlength(cLayout)), zeros(axlength(cLayout))

using blockSQP
using ForwardDiff

f = objective
g = shooting_constraints
grad_f = x->ForwardDiff.gradient(f, x)
jac_g = x->ForwardDiff.gradient(g, x)

prob = BlockSQPProblem(
    f, g, grad_f, jac_g,
    lb_var, ub_var, lb_con, ub_con,
    x0, zeros(axlength(vLayout) + axlenth(cLayout))
)
opts = sparse_options()
stats = SQPstats("./")

meth = SQPmethod(prob, opts, stats)

blockSQP.init!(meth)
blockSQP.run!(meth, 200)
blockSQP.finish!(meth)




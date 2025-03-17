module OptimizationExtension

using Optimization, blockSQP, Optimization.SciMLBase

SciMLBase.allowsbounds(::BlockSQPOpt) = true
SciMLBase.allowsconstraints(::BlockSQPOpt) = true
SciMLBase.allowscallback(::BlockSQPOpt) = false
SciMLBase.supports_opt_cache_interface(opt::BlockSQPOpt) = false

@info "Loading extension for blockSQP..."

function __map_optimizer_args!(prob::OptimizationProblem, opt::blockSQP.BlockSQPOptions;
    callback = nothing,
    maxiters::Union{Number, Nothing} = nothing,
    maxtime::Union{Number, Nothing} = nothing,
    abstol::Union{Number, Nothing} = nothing,
    reltol::Union{Number, Nothing} = nothing,
    kwargs...)
    for j in kwargs
        setproperty!(opt.options, j.first, j.second)
    end

    if !isnothing(maxiters)
        opt.maxiters = maxiters
    end

    if !isnothing(maxtime)
        @warn "common maxtime is currently not used by $(opt)"
    end

    if !isnothing(reltol)
        @warn "common reltol is currently not used by $(opt)"
    end

    if !isnothing(abstol)
        opt.opttol = abstol
    end

    return nothing
end

function SciMLBase.__solve(prob::OptimizationProblem,
    opt::BlockSQPOpt;
    callback = nothing,
    maxiters::Union{Number, Nothing} = nothing,
    maxtime::Union{Number, Nothing} = nothing,
    abstol::Union{Number, Nothing} = nothing,
    reltol::Union{Number, Nothing} = nothing,
    progress = false,
    options::Union{Nothing,blockSQP.BlockSQPOptions} = nothing,
    kwargs...)

    local x

    maxiters = Optimization._check_and_convert_maxiters(maxiters)
    maxtime = Optimization._check_and_convert_maxtime(maxtime)
    num_cons = prob.ucons === nothing ? 0 : length(prob.ucons)
    num_x = length(prob.u0)
    T = eltype(prob.u0)
    reinit_cache = OptimizationBase.ReInitCache(prob.u0, prob.p) # everything that can be changed via `reinit`
    f = Optimization.instantiate_function(
        prob.f, reinit_cache, prob.f.adtype, num_cons;
        g = true, h = true, cons_j = true, cons_h = true)


    _loss = function (θ)
        x = f.f(θ, prob.p)
        return first(x)
    end

    _g = function(θ)
        g_ = similar(θ)
        f.grad(g_, θ)
        return g_
    end

    _cons = num_cons > 0  ? function(θ)
        c_ = zeros(eltype(θ), num_cons)
        f.cons(c_, θ)
        return c_
    end : x -> zeros(eltype(x), 1) # dummy constraint that it doesnt crash, TODO: FIX THIS

    _jac_cons = num_cons > 0 ? function(θ)
        J = zeros(eltype(θ), num_cons, num_x)
        f.cons_j(J, θ)
        return J
    end : x -> zeros(eltype(x), 1, num_x)

    num_cons = max(1, num_cons)
    _lb = isnothing(prob.lb) ? -Inf * ones(T, num_x) : prob.lb
    _ub = isnothing(prob.ub) ? Inf * ones(T, num_x) : prob.ub
    _lb_cons = isnothing(prob.lcons) ? -Inf * ones(T, num_cons) : prob.lcons
    _ub_cons = isnothing(prob.ucons) ? Inf * ones(T, num_cons) : prob.ucons

    _lambda_0 = zeros(T, num_x+num_cons)
    opts = isnothing(options) ? blockSQP.BlockSQPOptions() : options

    __map_optimizer_args!(prob, opts, callback = callback,
                maxiters = maxiters, maxtime = maxtime,
                abstol = abstol, reltol = reltol,
                ; kwargs...)
    stats = blockSQP.SQPstats("./")

    _lb = blockSQP.__lowerbounds(_lb)
    _ub = blockSQP.__upperbounds(_ub)
    _u0 = blockSQP.__initial_values(prob.u0)

    sqp_prob = blockSQP.blockSQPProblem(_loss, _cons, _g, _jac_cons,
                            _lb, _ub, _lb_cons, _ub_cons,
                            _u0, _lambda_0)

    t0 = time()
        meth = blockSQP.Solver(sqp_prob, opts, stats)

        blockSQP.init(meth)

        ret = blockSQP.run(meth, opts.maxiters, 1)

        blockSQP.finish(meth)
    t1 = time()

    x_opt = blockSQP.get_primal_solution(meth)
    f_opt = _loss(x_opt)
    retcode = ret == 0 ? SciMLBase.ReturnCode.Success : SciMLBase.ReturnCode.Default

    SciMLBase.build_solution(SciMLBase.DefaultOptimizationCache(prob.f, prob.p), opt,
    x_opt, f_opt;
         (; original = ret, retcode = retcode,
            solve_time = t1 - t0)...)
end

end # module OptimizationExtension

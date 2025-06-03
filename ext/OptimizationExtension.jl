module OptimizationExtension

using blockSQP, blockSQP.SparseArrays
using Optimization, Optimization.SciMLBase
using Pkg

SciMLBase.allowsbounds(::BlockSQPOpt) = true
SciMLBase.allowsconstraints(::BlockSQPOpt) = true
SciMLBase.allowscallback(::BlockSQPOpt) = true
SciMLBase.requiresgradient(::BlockSQPOpt) = true
SciMLBase.requiresconsjac(::BlockSQPOpt) = true
SciMLBase.supports_opt_cache_interface(opt::BlockSQPOpt) = true

@info "Loading Optimization.jl extension for blockSQP..."

function SciMLBase.__init(prob::SciMLBase.OptimizationProblem, opt::BlockSQPOpt,
        ;
        callback = Optimization.DEFAULT_CALLBACK,
        progress = false, maxiters=nothing,
        options::BlockSQPOptions=BlockSQPOptions(),
        sparsity::Union{Vector{<:Integer}, Bool}=false, kwargs...)

        use_sparse_functions = sparsity != false || options.sparseQP == 2
        num_x = prob.u0 |> length
        num_cons = prob.ucons === nothing ? 0 : length(prob.ucons)
        @info sparsity
        blocks_hess = begin
            if use_sparse_functions
                if isa(sparsity, Bool)
                    function cons_ip(cons,x)
                        prob.f.cons(cons, x, prob.p)
                        return cons
                    end
                    blockSQP.compute_hessian_blocks(prob.f.f, cons_ip, num_x, num_cons; parameters=prob.p)
                else
                    @assert (sparsity[1] == 0) && (sparsity[end] == num_x) "sparsity[1] must be 0, sparsity[num_vars+1] must be num_vars"
                    sparsity
                end
            else
                [0, num_x]
            end
        end
        @info blocks_hess
    return OptimizationCache(prob, opt; callback = callback, progress = progress,
        use_sparse_functions = use_sparse_functions, sparsity = blocks_hess, options=options,
        maxiters = maxiters, kwargs...)
end

function __map_optimizer_args!(cache::OptimizationCache,
    opt::blockSQP.BlockSQPOptions;
    callback = nothing,
    maxiters::Union{Number, Nothing} = nothing,
    maxtime::Union{Number, Nothing} = nothing,
    abstol::Union{Number, Nothing} = nothing,
    reltol::Union{Number, Nothing} = nothing,
    sparsity::AbstractVector{Int},
    kwargs...)

    for j in kwargs
        if j.first in fieldnames(blockSQP.BlockSQPOptions)
            setproperty!(opt, j.first, j.second)
        end
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

    if length(sparsity) > 2
        opt.hessUpdate = 1
        opt.sparseQP = 2
    end

    return nothing
end


function SciMLBase.__solve(
    cache::OptimizationCache{F,RC,LB,UB,LC,UC,S,O,D,P,C}
 ) where {F,RC,LB,UB,LC,UC,S,O <:BlockSQPOpt,D,P,C}

    local x

    maxiters = Optimization._check_and_convert_maxiters(cache.solver_args.maxiters)
    maxtime = Optimization._check_and_convert_maxtime(cache.solver_args.maxtime)
    num_cons = cache.ucons === nothing ? 0 : length(cache.ucons)
    num_x = length(cache.u0)
    T = eltype(cache.u0)

    use_sparse_functions = cache.solver_args.use_sparse_functions
    sparsity = cache.solver_args.sparsity

    callback = cache.callback

    global iter = 0
    global x_prev = similar(cache.u0)
    _loss = function (θ)
        x = cache.f(θ, cache.p)
        if !isequal(x_prev, θ)
            state = Optimization.OptimizationState(iter = iter,
                u = θ,
                objective = x[1])
            cache.callback(state, x[1])
            iter +=1
            x_prev .= θ
        end
        return only(x)
    end

    _g = function(θ)
        g_ = similar(θ)
        cache.f.grad(g_, θ)
        return g_
    end

    _cons = begin
        if num_cons > 0
            function(θ)
                c_ = zeros(eltype(θ), num_cons)
                cache.f.cons(c_, θ)
                return c_
            end
        else
            (x) -> zeros(eltype(x), 1) # dummy constraint that it doesnt crash, TODO: FIX THIS
        end
    end

    _jac_cons = begin
        if num_cons > 0
            function(θ)
                J = zeros(eltype(θ), num_cons, num_x)
                cache.f.cons_j(J, θ)
                return J
            end
        else
            (x) -> zeros(eltype(x), 1, num_x)
        end
    end

    jac_g_row, jac_g_col, nnz = begin
        if use_sparse_functions
            J_sparse = sparse(_jac_cons(cache.u0))
            J_sparse.rowval .- 1, J_sparse.colptr .- 1, length(J_sparse.nzval)
        else
            Int32[], Int32[], -1
        end
    end

    sparse_jac(x) = let row = jac_g_row, col = jac_g_col, nnz=nnz
        _J = sparse(_jac_cons(x))
        newrow = _J.rowval .- 1
        newcol = _J.colptr .- 1
        newnnz = length(_J.nzval)
        @assert all(row .== newrow)
        @assert all(col .== newcol)
        @assert nnz == newnnz
        _J.nzval
    end

    num_cons = max(1, num_cons)
    _lb = isnothing(cache.lb) ? -Inf * ones(T, num_x) : cache.lb
    _ub = isnothing(cache.ub) ? Inf * ones(T, num_x) : cache.ub
    _lb_cons = isnothing(cache.lcons) ? -Inf * ones(T, num_cons) : cache.lcons
    _ub_cons = isnothing(cache.ucons) ? Inf * ones(T, num_cons) : cache.ucons

    _lambda_0 = zeros(T, num_x+num_cons)
    opts = cache.solver_args.options

    __map_optimizer_args!(cache, opts, callback = callback,
                maxiters = maxiters, maxtime = maxtime,
                abstol = cache.solver_args.abstol, reltol = cache.solver_args.reltol,
                sparsity=sparsity; cache.solver_args...)

    stats = blockSQP.SQPstats("./")

    _lb = blockSQP.__lowerbounds(_lb)
    _ub = blockSQP.__upperbounds(_ub)
    _u0 = blockSQP.__initial_values(cache.u0)

    sqp_prob = blockSQP.blockSQPProblem(_loss, _cons, _g, _jac_cons,
                            _lb, _ub, _lb_cons, _ub_cons,
                            _u0, _lambda_0; blockIdx = sparsity,
                            jac_g_row = jac_g_row, jac_g_colind = jac_g_col,
                            nnz = nnz,
                            jac_g_nz = use_sparse_functions ? sparse_jac : blockSQP.fnothing
                            )

    t0 = time()
        meth = blockSQP.Solver(sqp_prob, opts, stats)

        blockSQP.init(meth)

        ret = blockSQP.run(meth, opts.maxiters, 1)

        blockSQP.finish(meth)
    t1 = time()

    x_opt = blockSQP.get_primal_solution(meth)
    lambda = blockSQP.get_dual_solution(meth)
    f_opt = _loss(x_opt)
    retcode = ret == 0 ? SciMLBase.ReturnCode.Success : SciMLBase.ReturnCode.Default

    SciMLBase.build_solution(cache, cache.opt,
    x_opt, f_opt;
         (; original = (ret = ret, multiplier = lambda) , retcode = retcode,
            solve_time = t1 - t0)...)
end

end # module OptimizationExtension

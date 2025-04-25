module OptimizationExtension

using blockSQP
using Optimization, Optimization.SciMLBase
using Pkg
using Symbolics, SparseDiffTools

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
    sparsity::Union{Bool, AbstractVector{Int}},
    kwargs...)

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

    if sparsity != false
        opt.hessUpdate = 1
        opt.sparseQP = 2
    end

    return nothing
end

function find_blocks_from_csc(A::SparseDiffTools.SparseArrays.SparseMatrixCSC)
    n = size(A,1) # assume quadratic matrix here, since it is a hessian
    blockIdx = [0]
    max_row = zeros(n)
    for i=1:n
        col_start, col_end = A.colptr[i], A.colptr[i+1]-1
        max_row[i] = maximum(A.rowval[col_start:col_end])
    end

    for i=1:n
        if max_row[i] <= i
            append!(blockIdx, i)
        end
    end
    return blockIdx
end

function SciMLBase.__solve(prob::OptimizationProblem,
    opt::BlockSQPOpt;
    callback = nothing,
    maxiters::Union{Number, Nothing} = nothing,
    maxtime::Union{Number, Nothing} = nothing,
    abstol::Union{Number, Nothing} = nothing,
    reltol::Union{Number, Nothing} = nothing,
    progress = false,
    sparsity::Union{Bool, AbstractVector{Int}}=false,
    options::Union{Nothing,blockSQP.BlockSQPOptions} = nothing,
    kwargs...)

    local x

    maxiters = Optimization._check_and_convert_maxiters(maxiters)
    maxtime = Optimization._check_and_convert_maxtime(maxtime)
    num_cons = prob.ucons === nothing ? 0 : length(prob.ucons)
    num_x = length(prob.u0)
    T = eltype(prob.u0)
    # Hacky solution for backward compatibility with DynamicOED. Will be removed shortly.
    f = begin
        if pkgversion(Optimization) < v"4"
            Optimization.instantiate_function(prob.f, prob.u0, prob.f.adtype, prob.p, num_cons)
        else
            reinit_cache = Optimization.ReInitCache(prob.u0, prob.p) # everything that can be changed via `reinit`
            Optimization.instantiate_function(
                prob.f, reinit_cache, prob.f.adtype, num_cons;
                g = true, h = false, cons_j = true, cons_h = false)
        end
    end

    sparsity_defined = sparsity != false
    blocks_hess = begin
        if isa(sparsity, Bool)
            lag(x, mu) = begin
                fx = f.f(x, prob.p)
                g = zeros(eltype(x), num_cons)
                f.cons(g, x)
                return fx + sum(mu .* g)
            end
            input = Vector{Float64}(undef, num_x);
            sparse_hess = Symbolics.hessian_sparsity(x -> lag(x, ones(num_cons)), input)
            find_blocks_from_csc(sparse_hess)
        else
            @assert length(sparsity) == num_x + 1
            @assert (sparsity[0] == 0) && (sparsity[end] == num_x) "sparsity[1] must be 0, sparsity[num_vars+1] must be num_vars"
            sparsity
        end
    end

    sparse_jac(x) = begin
        sparsity_defined || return nothing
        sparse_ad = AutoSparse(prob.f.adtype)
        sd = SymbolicsSparsityDetection()
        sparse_jacobian(sparse_ad, sd, f.cons, zeros(num_cons), x)
    end


    jac_g_row, jac_g_col, nnz = begin
        if !sparsity_defined
            Int32[], Int32[], -1
        else
            J_sparse = sparse_jac(prob.u0)
            J_sparse.rowval .- 1, J_sparse.colptr .- 1, length(J_sparse.nzval)
        end
    end

    _loss = function (θ)
        x = f.f(θ, prob.p)
        return first(x)
    end

    _g = function(θ)
        g_ = similar(θ)
        f.grad(g_, θ)
        return g_
    end

    _cons = begin
        if num_cons > 0
            function(θ)
                c_ = zeros(eltype(θ), num_cons)
                f.cons(c_, θ)
                return c_
            end
        else
            (x) -> zeros(eltype(x), 1) # dummy constraint that it doesnt crash, TODO: FIX THIS
        end
    end

    _jac_cons = begin
        if !sparsity_defined
            if num_cons > 0
                function(θ)
                    J = zeros(eltype(θ), num_cons, num_x)
                    f.cons_j(J, θ)
                    return J
                end
            else
                (x) -> zeros(eltype(x), 1, num_x)
            end
        else
            (x) -> sparse_jac(x).nzval
        end
    end

    num_cons = max(1, num_cons)
    _lb = isnothing(prob.lb) ? -Inf * ones(T, num_x) : prob.lb
    _ub = isnothing(prob.ub) ? Inf * ones(T, num_x) : prob.ub
    _lb_cons = isnothing(prob.lcons) ? -Inf * ones(T, num_cons) : prob.lcons
    _ub_cons = isnothing(prob.ucons) ? Inf * ones(T, num_cons) : prob.ucons

    _lambda_0 = zeros(T, num_x+num_cons)
    opts = isnothing(options) ? blockSQP.BlockSQPOptions() : options

    __map_optimizer_args!(prob, opts, callback = callback,
                maxiters = maxiters, maxtime = maxtime,
                abstol = abstol, reltol = reltol, sparsity=sparsity,
                ; kwargs...)

    stats = blockSQP.SQPstats("./")

    _lb = blockSQP.__lowerbounds(_lb)
    _ub = blockSQP.__upperbounds(_ub)
    _u0 = blockSQP.__initial_values(prob.u0)

    sqp_prob = blockSQP.blockSQPProblem(_loss, _cons, _g, _jac_cons,
                            _lb, _ub, _lb_cons, _ub_cons,
                            _u0, _lambda_0; blockIdx = blocks_hess,
                            jac_g_row = jac_g_row, jac_g_colind = jac_g_col,
                            nnz = nnz, jac_g_nz = sparsity_defined ? _jac_cons : blockSQP.fnothing
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

    SciMLBase.build_solution(SciMLBase.DefaultOptimizationCache(prob.f, prob.p), opt,
    x_opt, f_opt;
         (; original = (ret = ret, multiplier = lambda) , retcode = retcode,
            solve_time = t1 - t0)...)
end

end # module OptimizationExtension

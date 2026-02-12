module OptimizationExtension

using blockSQP2, blockSQP2.SparseArrays, blockSQP2.NLPlayouts
using OptimizationBase, OptimizationBase.SciMLBase

SciMLBase.allowsbounds(::BlockSQP2Optimizer) = true
SciMLBase.allowsconstraints(::BlockSQP2Optimizer) = true
SciMLBase.allowscallback(::BlockSQP2Optimizer) = true
SciMLBase.requiresgradient(::BlockSQP2Optimizer) = true
SciMLBase.requiresconsjac(::BlockSQP2Optimizer) = true
SciMLBase.supports_opt_cache_interface(::BlockSQP2Optimizer) = true
SciMLBase.has_init(::BlockSQP2Optimizer) = true

@info "Loading Optimization(Base).jl extension for blockSQP2..."

function SciMLBase.__init(
            prob::SciMLBase.OptimizationProblem, opt::BlockSQP2Optimizer,
            ;
            callback = OptimizationBase.DEFAULT_CALLBACK,
            progress = false, maxiters = nothing,
            options::blockSQP2.Options=blockSQP2.Options(),
            blockIdx::Union{Vector{<:Integer}, Nothing} = nothing,
            vblocks::Union{Vector{blockSQP2.vblock}, Nothing} = blockSQP2.vblock[],
            condenser::Union{blockSQP2.Condenser, Nothing} = nothing,
            layout::Union{NLPlayouts.NLPlayout, Nothing} = nothing,
            kwargs...
            )
    if isnothing(vblocks)
        vblocks = blockSQP2.vblock[]
    end
    return OptimizationCache(prob, opt; 
        callback = callback, progress = progress, maxiters = maxiters,
        options=options, blockIdx = blockIdx, vblocks = vblocks, layout = layout,
        condenser = condenser, kwargs...)
end

function __map_optimizer_args!(cache::OptimizationCache,
    opt::blockSQP2.Options;
    callback = nothing,
    maxiters::Union{Number, Nothing} = nothing,
    maxtime::Union{Number, Nothing} = nothing,
    abstol::Union{Number, Nothing} = nothing,
    reltol::Union{Number, Nothing} = nothing,
    kwargs...)

    for j in kwargs
        if j.first in fieldnames(blockSQP2.Options)
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
        opt.opt_tol = abstol
    end

    return nothing
end


function SciMLBase.__solve(
    cache::OptimizationCache{O,IIP,F,RC,LB,UB,LC,UC,S,P,C,M}
 ) where {O<:BlockSQP2Optimizer,IIP,F,RC,LB,UB,LC,UC,S,P,C,M}
    local x
    
    num_cons = cache.ucons === nothing ? 0 : length(cache.ucons)
    num_x = length(cache.u0)
    T = eltype(cache.u0)
    
    opts = cache.solver_args.options
    callback = cache.callback
    maxiters = OptimizationBase._check_and_convert_maxiters(cache.solver_args.maxiters)
    maxtime = OptimizationBase._check_and_convert_maxtime(cache.solver_args.maxtime)
    
    __map_optimizer_args!(cache, opts; callback = callback,
                maxiters = maxiters, maxtime = maxtime,
                abstol = cache.solver_args.abstol, reltol = cache.solver_args.reltol,
                cache.solver_args...)
    
    _blockIdx = Cint[0, num_x]
    _vblocks = blockSQP2.vblock[]
    _condenser = nothing
    
    if !isnothing(cache.solver_args.layout) 
        _layout = cache.solver_args.layout
        _blockIdx = hessBlockIndexZeroBased(_layout)
        _vblocks = create_vblocks(_layer)
        
        #Deactivate this for now, requiring explicit passing of a condenser.
        # _condenser = blockSQP2.Condenser(_layout)
    end
    if !isnothing(cache.solver_args.condenser)
        _condenser = cache.solver_args.condenser
    end
    if length(cache.solver_args.vblocks) > 0
        _vblocks = cache.solver_args.vblocks
    end
    if !isnothing(cache.solver_args.blockIdx)
        _blockIdx = cache.solver_args.blockIdx
    end
    
    _loss = function (θ)
        x = cache.f(θ, cache.p)
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
    
    jac_g_row, jac_g_col, nnz, jac_row, jac_col = begin
        # if use_sparse_functions
        if opts.sparse
            #Hacky: Calculate constraint Jacobian for perturbed points to find structural nonzero elements
            u0_pert_1 = [x + 1e-6*rand() for x in cache.u0]
            u0_pert_2 = [x + 1e-5*rand() for x in cache.u0]
            u0_pert_3 = [x + 1e-4*rand() for x in cache.u0]
            J_sparse = sparse(_jac_cons(cache.u0) + _jac_cons(u0_pert_1) + _jac_cons(u0_pert_2) + _jac_cons(u0_pert_3))
            jac_row, jac_col, jac_val = findnz(J_sparse)
            J_sparse.rowval .- 1, J_sparse.colptr .- 1, length(jac_val), jac_row, jac_col
        else
            Cint[], Cint[], -1, nothing, nothing
        end
    end
    
    sparse_jac(x) = let jac_row = jac_row, jac_col = jac_col
        _J = _jac_cons(x)
        return [_J[i,j] for (i,j) in zip(jac_row,jac_col)]
    end
    
    num_cons = max(1, num_cons)
    _lb = isnothing(cache.lb) ? -Inf * ones(T, num_x) : cache.lb
    _ub = isnothing(cache.ub) ? Inf * ones(T, num_x) : cache.ub
    _lb_cons = isnothing(cache.lcons) ? -Inf * ones(T, num_cons) : cache.lcons
    _ub_cons = isnothing(cache.ucons) ? Inf * ones(T, num_cons) : cache.ucons
    
    _lambda_0 = zeros(T, num_x+num_cons)
    
    
    stats = blockSQP2.Stats("./")
    
    _lb = blockSQP2.__lowerbounds(_lb)
    _ub = blockSQP2.__upperbounds(_ub)
    _u0 = blockSQP2.__initial_values(cache.u0)
    
    sqp_prob = blockSQP2.Problem(_loss, _cons, _g, _jac_cons,
                            _lb, _ub, _lb_cons, _ub_cons,
                            _u0, _lambda_0; 
                            blockIdx = _blockIdx,
                            vblocks = _vblocks,
                            condenser = _condenser,
                            jac_g_row = jac_g_row, jac_g_colind = jac_g_col,
                            nnz = nnz,
                            jac_g_nz = opts.sparse ? sparse_jac : blockSQP2.fnothing
                            )
    
    meth = blockSQP2.Solver(sqp_prob, opts, stats)
    
    blockSQP2.init!(meth)
    t0 = time()
    if cache.callback == OptimizationBase.DEFAULT_CALLBACK
        ret = blockSQP2.run!(meth, opts.maxiters, 1)
    else
        for i=1:opts.maxiters
            ret = blockSQP2.run!(meth, 1, 1)
            _iterate = blockSQP2.get_primal_solution(meth)
            _obj = _loss(_iterate)[1]
            state = OptimizationBase.OptimizationState(iter = i,
                u = _iterate,
                objective = _obj[1])
            cret = cache.callback(state, _obj)
            (cret || blockSQP2.is_success(ret)) && break
        end
    end
    blockSQP2.finish!(meth)
    t1 = time()
    
    x_opt = blockSQP2.get_primal_solution(meth)
    lambda = blockSQP2.get_dual_solution(meth)
    f_opt = _loss(x_opt)
    retcode = blockSQP2.is_success(ret) ? SciMLBase.ReturnCode.Success : SciMLBase.ReturnCode.Default
    
    SciMLBase.build_solution(cache, cache.opt,
    x_opt, f_opt;
         (; original = (ret = ret, multiplier = lambda, solve_time = t1 - t0, solve_it = Int64(blockSQP2.get_itCount(meth))) , retcode = retcode,
            )...)
end

end # module OptimizationExtension

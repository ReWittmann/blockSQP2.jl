abstract type QPsolver_options end

mutable struct blockSQPOptions
    maxiters::Cint
    eps::Cdouble
    inf::Cdouble
    print_level::Cint
    result_print_color::Cint
    debug_level::Cint
    opt_tol::Cdouble
    feas_tol::Cdouble
    sparse::Bool
    enable_rest::Bool
    lim_mem::Bool
    mem_size::Cint
    block_hess::Cint
    exact_hess::Cint
    hess_approx::Cint
    fallback_approx::Cint
    initial_hess_scale::Cdouble
    sizing::Cint
    fallback_sizing::Cint
    COL_eps::Cdouble
    COL_tau_1::Cdouble
    COL_tau_2::Cdouble
    OL_eps::Cdouble
    BFGS_damping_factor::Cdouble
    conv_strategy::Cint
    max_conv_QPs::Cint
    enable_linesearch::Bool
    max_linesearch_steps::Cint
    max_consec_reduced_steps::Cint
    max_consec_skipped_updates::Cint
    skip_first_linesearch::Cint
    max_SOC::Cint
    qpsol::Union{String, Symbol, Vector{Cchar}}
    qpsol_options::Union{QPsolver_options, Nothing}
    max_QP_it::Cint
    max_QP_secs::Cdouble
    max_extra_steps::Cint
    max_filter_overrides::Cint
    par_QPs::Bool
    enable_QP_cancellation::Bool
    automatic_scaling::Bool
    enable_premature_termination::Bool
    function blockSQPOptions(;
        maxiters::Integer = 100,
        eps::AbstractFloat = 1.0e-16,
        inf::AbstractFloat = Inf,
        print_level::Integer = 2,
        result_print_color::Integer = 2,
        debug_level::Integer = 0,
        opt_tol::AbstractFloat = 1.0e-6,
        feas_tol::AbstractFloat = 1.0e-6,
        sparse::Bool = true,
        enable_rest::Bool = true,
        lim_mem::Bool = true,
        mem_size::Integer = 20,
        block_hess::Integer = 1,
        exact_hess::Integer = 0,
        hess_approx::Integer = 1,
        fallback_approx::Integer = 2,
        initial_hess_scale::AbstractFloat = 1.0,
        sizing::Integer = 2,
        fallback_sizing::Integer = 4,
        COL_eps::AbstractFloat = 0.1,
        COL_tau_1::AbstractFloat = 0.5,
        COL_tau_2::AbstractFloat = 1.0e4,
        OL_eps::AbstractFloat = 1.0e-4,
        BFGS_damping_factor::AbstractFloat = 1/3,
        conv_strategy::Integer = 1,
        max_conv_QPs::Integer = 4,
        enable_linesearch::Bool = true,
        max_linesearch_steps::Integer = 10,
        max_consec_reduced_steps::Integer = 8,
        max_consec_skipped_updates::Integer = 100,
        skip_first_linesearch::Bool = false,
        max_SOC::Integer = 3,
        qpsol::String = "qpOASES",
        qpsol_options::Union{QPsolver_options, Nothing} = nothing,
        max_QP_it::Integer = 5000,
        max_QP_secs::AbstractFloat = 3600.0,
        max_extra_steps::Integer = 0,
        max_filter_overrides::Integer = 2,
        par_QPs::Bool = false,
        enable_QP_cancellation::Bool = true,
        automatic_scaling::Bool = false,
        enable_premature_termination::Bool = false
    )
        return new(
            maxiters,
            eps,
            inf,
            print_level,
            result_print_color,
            debug_level,
            opt_tol,
            feas_tol,
            sparse,
            enable_rest,
            lim_mem,
            mem_size,
            block_hess,
            exact_hess,
            hess_approx,
            fallback_approx,
            initial_hess_scale,
            sizing,
            fallback_sizing,
            COL_eps,
            COL_tau_1,
            COL_tau_2,
            OL_eps,
            BFGS_damping_factor,
            conv_strategy,
            max_conv_QPs,
            enable_linesearch,
            max_linesearch_steps,
            max_consec_reduced_steps,
            max_consec_skipped_updates,
            skip_first_linesearch,
            max_SOC,
            qpsol,
            qpsol_options,
            max_QP_it,
            max_QP_secs,
            max_extra_steps,
            max_filter_overrides,
            par_QPs,
            enable_QP_cancellation,
            automatic_scaling,
            enable_premature_termination
        )
    end
end


mutable struct qpOASES_options <: QPsolver_options
    sparsityLevel::Cint
    printLevel::Cint
    terminationTolerance::Cdouble
    function qpOASES_options(;
        sparsityLevel::Integer = 2,
        printLevel::Integer = 0,
        terminationTolerance::AbstractFloat = 5.0e6*2.221e-16
        )
        return new(sparsityLevel, printLevel, terminationTolerance)
    end
end


function create_cxx_options(opts::blockSQPOptions)
    BSQP = libblockSQP[]
    SQPoptions_obj::Ptr{Cvoid} = ccall(@dlsym(BSQP, "create_SQPoptions"), Ptr{Cvoid}, ())
    QPsolver_options_obj::Ptr{Cvoid} = Ptr{Cvoid}()
    # Constants
    ccall(@dlsym(BSQP, "SQPoptions_set_eps"), Cvoid, (Ptr{Cvoid}, Cdouble), SQPoptions_obj, Cdouble(opts.eps))
    ccall(@dlsym(BSQP, "SQPoptions_set_inf"), Cvoid, (Ptr{Cvoid}, Cdouble), SQPoptions_obj, Cdouble(opts.inf))
    
    # Output
    ccall(@dlsym(BSQP, "SQPoptions_set_print_level"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.print_level))
    ccall(@dlsym(BSQP, "SQPoptions_set_result_print_color"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.result_print_color))
    ccall(@dlsym(BSQP, "SQPoptions_set_debug_level"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.debug_level))
    
    # Termination criteria
    ccall(@dlsym(BSQP, "SQPoptions_set_opt_tol"), Cvoid, (Ptr{Cvoid}, Cdouble), SQPoptions_obj, Cdouble(opts.opt_tol))
    ccall(@dlsym(BSQP, "SQPoptions_set_feas_tol"), Cvoid, (Ptr{Cvoid}, Cdouble), SQPoptions_obj, Cdouble(opts.feas_tol))
    ccall(@dlsym(BSQP, "SQPoptions_set_enable_premature_termination"), Cvoid, (Ptr{Cvoid}, Cchar), SQPoptions_obj, Cchar(opts.automatic_scaling))
    ccall(@dlsym(BSQP, "SQPoptions_set_max_extra_steps"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.max_extra_steps))
    
    # Line search heuristics
    ccall(@dlsym(BSQP, "SQPoptions_set_max_filter_overrides"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.max_filter_overrides))

    # Derivative evaluation
    ccall(@dlsym(BSQP, "SQPoptions_set_sparse"), Cvoid, (Ptr{Cvoid}, Cchar), SQPoptions_obj, Cchar(opts.sparse))
    
    # Restoration phase
    ccall(@dlsym(BSQP, "SQPoptions_set_enable_rest"), Cvoid, (Ptr{Cvoid}, Cchar), SQPoptions_obj, Cchar(opts.enable_rest))
    
    # Full/limited memory quasi newton
    ccall(@dlsym(BSQP, "SQPoptions_set_lim_mem"), Cvoid, (Ptr{Cvoid}, Cchar), SQPoptions_obj, Cchar(opts.lim_mem))
    ccall(@dlsym(BSQP, "SQPoptions_set_mem_size"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.mem_size))
    
    # Hessian approximation
    ccall(@dlsym(BSQP, "SQPoptions_set_block_hess"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.block_hess))
    ccall(@dlsym(BSQP, "SQPoptions_set_exact_hess"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.exact_hess))
    ccall(@dlsym(BSQP, "SQPoptions_set_hess_approx"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.hess_approx))
    ccall(@dlsym(BSQP, "SQPoptions_set_fallback_approx"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.fallback_approx))
    
    # Hessian sizing
    ccall(@dlsym(BSQP, "SQPoptions_set_initial_hess_scale"), Cvoid, (Ptr{Cvoid}, Cdouble), SQPoptions_obj, Cdouble(opts.initial_hess_scale))
    ccall(@dlsym(BSQP, "SQPoptions_set_sizing"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.sizing))
    ccall(@dlsym(BSQP, "SQPoptions_set_fallback_sizing"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.fallback_sizing))
    ccall(@dlsym(BSQP, "SQPoptions_set_COL_eps"), Cvoid, (Ptr{Cvoid}, Cdouble), SQPoptions_obj, Cdouble(opts.COL_eps))
    ccall(@dlsym(BSQP, "SQPoptions_set_COL_tau_1"), Cvoid, (Ptr{Cvoid}, Cdouble), SQPoptions_obj, Cdouble(opts.COL_tau_1))
    ccall(@dlsym(BSQP, "SQPoptions_set_COL_tau_2"), Cvoid, (Ptr{Cvoid}, Cdouble), SQPoptions_obj, Cdouble(opts.COL_tau_2))
    ccall(@dlsym(BSQP, "SQPoptions_set_OL_eps"), Cvoid, (Ptr{Cvoid}, Cdouble), SQPoptions_obj, Cdouble(opts.OL_eps))
    
    # Quasi-Newton
    ccall(@dlsym(BSQP, "SQPoptions_set_BFGS_damping_factor"), Cvoid, (Ptr{Cvoid}, Cdouble), SQPoptions_obj, Cdouble(opts.BFGS_damping_factor))
    
    # Convexification strategy
    ccall(@dlsym(BSQP, "SQPoptions_set_conv_strategy"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.conv_strategy))
    ccall(@dlsym(BSQP, "SQPoptions_set_max_conv_QPs"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.max_conv_QPs))
    ccall(@dlsym(BSQP, "SQPoptions_set_par_QPs"), Cvoid, (Ptr{Cvoid}, Cchar), SQPoptions_obj, Cchar(opts.par_QPs))
    ccall(@dlsym(BSQP, "SQPoptions_set_enable_QP_cancellation"), Cvoid, (Ptr{Cvoid}, Cchar), SQPoptions_obj, Cchar(opts.enable_QP_cancellation))
    
    # Scaling
    ccall(@dlsym(BSQP, "SQPoptions_set_automatic_scaling"), Cvoid, (Ptr{Cvoid}, Cchar), SQPoptions_obj, Cchar(opts.automatic_scaling))
    
    # Filter line search
    ccall(@dlsym(BSQP, "SQPoptions_set_enable_linesearch"), Cvoid, (Ptr{Cvoid}, Cchar), SQPoptions_obj, Cchar(opts.enable_linesearch))
    ccall(@dlsym(BSQP, "SQPoptions_set_max_linesearch_steps"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.max_linesearch_steps))
    ccall(@dlsym(BSQP, "SQPoptions_set_max_consec_reduced_steps"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.max_consec_reduced_steps))
    ccall(@dlsym(BSQP, "SQPoptions_set_max_consec_skipped_updates"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.max_consec_skipped_updates))
    ccall(@dlsym(BSQP, "SQPoptions_set_skip_first_linesearch"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.skip_first_linesearch))
    ccall(@dlsym(BSQP, "SQPoptions_set_max_SOC"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.max_SOC))
    
    #qpsol and qpsol_options below
    ccall(@dlsym(BSQP, "SQPoptions_set_max_QP_it"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(opts.max_QP_it))
    ccall(@dlsym(BSQP, "SQPoptions_set_max_QP_secs"), Cvoid, (Ptr{Cvoid}, Cdouble), SQPoptions_obj, Cdouble(opts.max_QP_secs))
    if (opts.qpsol == "qpOASES" || opts.qpsol == :qpOASES || opts.qpsol == Cchar['q', 'p', 'O', 'A', 'S', 'E', 'S', '\0'])
        ccall(@dlsym(BSQP, "SQPoptions_set_qpsol"), Cvoid, (Ptr{Cvoid}, Cint), SQPoptions_obj, Cint(0))
    end
    if typeof(opts.qpsol_options) == qpOASES_options
        QPsolver_options_obj = ccall(@dlsym(BSQP, "create_qpOASES_options"), Ptr{Cvoid}, ())
        ccall(@dlsym(BSQP, "qpOASES_options_set_sparsityLevel"), Cvoid, (Ptr{Cvoid}, Cint), QPsolver_options_obj, Cint(opts.qpsol_options.sparsityLevel))
        ccall(@dlsym(BSQP, "qpOASES_options_set_printLevel"), Cvoid, (Ptr{Cvoid}, Cint), QPsolver_options_obj, Cint(opts.qpsol_options.printLevel))
        ccall(@dlsym(BSQP, "qpOASES_options_set_terminationTolerance"), Cvoid, (Ptr{Cvoid}, Cdouble), QPsolver_options_obj, Cdouble(opts.qpsol_options.terminationTolerance))
        ccall(@dlsym(BSQP, "SQPoptions_set_qpsol"), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), SQPoptions_obj, QPsolver_options_obj)
    # elseif typeof(opts.qpsol_options) == ...
    end
    
    return SQPoptions_obj, QPsolver_options_obj
end

function sparse_options()
    opts = BlockSQPOptions()
    opts.sparse = true
    opts.enable_linesearch = true
    opts.hess_approx = 1
    opts.sizing = 2
    opts.fallback_approx = 2
    opts.fallback_sizing = 4
    opts.opt_tol = 1e-6
    opts.feas_tol = 1e-6
    return opts
end

#=
mutable struct blockSQPOptions
    opttol::Float64
    nlinfeastol::Float64
    maxiters::Int32
    globalization::Int32
    hessUpdate::Int32
    fallbackUpdate::Int32
    hessScaling::Int32
    fallbackScaling::Int32
    hessLimMem::Int32
    hessMemsize::Int32
    maxConsecSkippedUpdates::Int32
    blockHess::Int32
    whichSecondDerv::Int32
    sparseQP::Int32
    printLevel::Int32
    printColor::Int32
    debugLevel::Int32
    which_QPsolver::String
    maxSOCiter::Int32
    maxConvQP::Int32
    convStrategy::Int32
    skipFirstGlobalization::Int32
    hessDampFac::Float64
    hessDamp::Int32
    colEps::Float64
    colTau1::Float64
    colTau2::Float64
    eps::Float64
    inf::Float64
    restoreFeas::Int32
    maxLineSearch::Int32
    maxConsecReducedSteps::Int32
    maxItQP::Int32
    maxTimeQP::Float64
    iniHessDiag::Float64
    function BlockSQPOptions(;
        blockHess::Integer = Int32(1),

        colEps::Float64=0.0,
        colTau1::Float64=0.0,
        colTau2::Float64=0.0,
        convStrategy::Integer=Int32(0),

        debugLevel::Integer = Int32(0),

        fallbackScaling::Integer = Int32(0),
        fallbackUpdate::Integer = Int32(2),

        globalization::Integer = Int32(0),
        hessDamp::Integer=Int32(1),

        hessDampFac::Float64=0.0,
        hessLimMem::Integer = Int32(1),
        hessMemsize::Integer = Int32(20),
        hessScaling::Integer = Int32(0),
        hessUpdate::Integer = Int32(1),

                            maxiters::Integer=100,
                            maxLineSearch::Integer=Int32(20),
                            maxConsecReducedSteps::Integer=Int32(100),
                            maxItQP::Integer=Int32(5000),
                            maxTimeQP::Float64=10000.0,
                            maxSOCiter::Integer=Int32(3),
                            maxConvQP::Integer=Int32(1),
                            maxConsecSkippedUpdates::Int32 = Int32(200),

                            nlinfeastol::Float64 = 1.0e-12,
                            opttol::Float64 = 1.0e-5,


                            printLevel::Integer = Int32(2),
                            printColor::Integer=Int32(0),
                            restoreFeas::Integer=Int32(1),

                            skipFirstGlobalization::Integer=Int32(1),
                            sparseQP::Integer = Int32(0),

                            which_QPsolver::String = "qpOASES",
                            whichSecondDerv::Integer = Int32(0),

                            eps::Float64=1e-16,
                            inf::Float64=1e20,


                            iniHessDiag::Float64=1.0
                            )
        return new(
                opttol,
                nlinfeastol,
                maxiters,
                globalization,
                hessUpdate,
                fallbackUpdate,
                hessScaling,
                fallbackScaling,
                hessLimMem,
                hessMemsize,
                maxConsecSkippedUpdates,
                blockHess,
                whichSecondDerv,
                sparseQP,
                printLevel,
                printColor,
                debugLevel,
                which_QPsolver,
                maxSOCiter,
                maxConvQP,
                convStrategy,
                skipFirstGlobalization,
                hessDampFac,
                hessDamp,
                colEps,
                colTau1,
                colTau2,
                eps,
                inf,
                restoreFeas,
                maxLineSearch,
                maxConsecReducedSteps,
                maxItQP,
                maxTimeQP,
                iniHessDiag)
    end
end
=#
abstract type QPsolver_options end

#=
mutable struct blockSQPOptions
    opt_tol::Float64
    feas_tol::Float64
    maxiters::Int32
    enable_linesearch::CxxBool
    hess_approx::Int32
    fallback_approx::Int32
    sizing::Int32
    fallback_sizing::Int32
    lim_mem::CxxBool
    mem_size::Int32
    max_consec_skipped_updates::Int32
    block_hess::Int32
    exact_hess::Int32
    sparse::Bool
    print_level::Int32
    result_print_color::Int32
    debug_level::Int32
    qpsol::String
    qpsol_options::Union{QPsolver_options, Nothing}
    max_SOC::Int32
    max_conv_QPs::Int32
    conv_strategy::Int32
    skip_first_linesearch::CxxBool
    BFGS_damping_factor::Float64
    COL_eps::Float64
    COL_tau_1::Float64
    COL_tau_2::Float64
    OL_eps::Float64
    eps::Float64
    inf::Float64
    enable_rest::CxxBool
    max_linesearch_steps::Int32
    max_consec_reduced_steps::Int32
    max_QP_it::Int32
    max_QP_secs::Float64
    initial_hess_scale::Float64
    function blockSQPOptions(;
                            opt_tol::Float64 = 1.0e-6,
                            feas_tol::Float64 = 1.0e-6,
                            maxiters::Integer=100,
                            enable_linesearch::Bool = true,
                            hess_approx::Integer = Int32(1),
                            fallback_approx::Integer = Int32(2),
                            sizing::Integer = Int32(1),
                            fallback_sizing::Integer = Int32(2),
                            lim_mem::Bool = true,
                            mem_size::Integer = Int32(20),
                            max_consec_skipped_updates::Int32 = Int32(200),
                            block_hess::Integer = Int32(1),
                            exact_hess::Integer = Int32(0),
                            sparse::Bool = true,
                            print_level::Integer = Int32(2),
                            result_print_color::Integer=Int32(2),
                            debug_level::Integer = Int32(0),
                            qpsol::String = "qpOASES",
                            qpsol_options::Union{QPsolver_options, Nothing} = nothing,
                            max_SOC::Integer=Int32(3),
                            max_conv_QPs::Integer=Int32(1),
                            conv_strategy::Integer=Int32(0),
                            skip_first_linesearch::Bool=false,
                            BFGS_damping_factor::Float64=1/3,
                            COL_eps::Float64=0.1,
                            COL_tau_1::Float64=0.5,
                            COL_tau_2::Float64=1.0e4,
                            OL_eps::Float64=1.0e-4,
                            eps::Float64=1e-16,
                            inf::Float64=Inf,
                            enable_rest::Bool=true,
                            max_linesearch_steps::Integer=Int32(20),
                            max_consec_reduced_steps::Integer=Int32(100),
                            max_QP_it::Integer=Int32(5000),
                            max_QP_secs::Float64=120.0,
                            initial_hess_scale::Float64=1.0
                            )
        return new(
                opt_tol,
                feas_tol,
                maxiters,
                enable_linesearch,
                hess_approx,
                fallback_approx,
                sizing,
                fallback_sizing,
                lim_mem,
                mem_size,
                max_consec_skipped_updates,
                block_hess,
                exact_hess,
                sparse,
                print_level,
                result_print_color,
                debug_level,
                qpsol,
                qpsol_options,
                max_SOC,
                max_conv_QPs,
                conv_strategy,
                skip_first_linesearch,
                BFGS_damping_factor,
                COL_eps,
                COL_tau_1,
                COL_tau_2,
                OL_eps,
                eps,
                inf,
                enable_rest,
                max_linesearch_steps,
                max_consec_reduced_steps,
                max_QP_it,
                max_QP_secs,
                initial_hess_scale)
    end
end

function sparse_options()
    opts = BlockSQPOptions()
    opts.sparseQP = 2
    opts.globalization = 1
    opts.hessUpdate = 1
    opts.hessScaling = 2
    opts.fallbackUpdate = 2
    opts.fallbackScaling = 4
    opts.opttol = 1e-6
    opts.nlinfeastol = 1e-6
    return opts
end
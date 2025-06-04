mutable struct BlockSQPOptions
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

function set_cxx_options(opts::BlockSQPOptions)
    cxx_opts = SQPoptions()
    opt_keys = string.(fieldnames(typeof(opts)))

    if "printLevel" in opt_keys
        set_printLevel(cxx_opts, Int32(opts.printLevel))
    end
    if "printColor" in opt_keys
        set_printColor(cxx_opts, Int32(opts.printColor))
    end
    if "debugLevel" in opt_keys
        set_debugLevel(cxx_opts, Int32(opts.debugLevel))
    end
    if "eps" in opt_keys
        set_eps(cxx_opts, Float64(opts.eps))
    end
    if "inf" in opt_keys
        set_inf(cxx_opts, Float64(opts.inf))
    end
    if "opttol" in opt_keys
        set_opttol(cxx_opts, Float64(opts.opttol))
    end
    if "nlinfeastol" in opt_keys
        set_nlinfeastol(cxx_opts, Float64(opts.nlinfeastol))
    end
    if "sparseQP" in opt_keys
        set_sparseQP(cxx_opts, Int32(opts.sparseQP))
    end
    if "globalization" in opt_keys
        set_globalization(cxx_opts, Int32(opts.globalization))
    end
    if "restoreFeas" in opt_keys
        set_restoreFeas(cxx_opts, Int32(opts.restoreFeas))
    end
    if "maxLineSearch" in opt_keys
        set_maxLineSearch(cxx_opts, Int32(opts.maxLineSearch))
    end
    if "maxConsecReducedSteps" in opt_keys
        set_maxConsecReducedSteps(cxx_opts, Int32(opts.maxConsecReducedSteps))
    end
    if "maxConsecSkippedUpdates" in opt_keys
        set_maxConsecSkippedUpdates(cxx_opts, Int32(opts.maxConsecSkippedUpdates))
    end
    if "maxItQP" in opt_keys
        set_maxItQP(cxx_opts, Int32(opts.maxItQP))
    end
    if "blockHess" in opt_keys
        set_blockHess(cxx_opts, Int32(opts.blockHess))
    end
    if "hessScaling" in opt_keys
        set_hessScaling(cxx_opts, Int32(opts.hessScaling))
    end
    if "fallbackScaling" in opt_keys
        set_fallbackScaling(cxx_opts, Int32(opts.fallbackScaling))
    end
    if "maxTimeQP" in opt_keys
        set_maxTimeQP(cxx_opts, Float64(opts.maxTimeQP))
    end
    if "iniHessDiag" in opt_keys
        set_iniHessDiag(cxx_opts, Float64(opts.iniHessDiag))
    end
    if "colEps" in opt_keys
        set_colEps(cxx_opts, Float64(opts.colEps))
    end
    if "colTau1" in opt_keys
        set_colTau1(cxx_opts, Float64(opts.colTau1))
    end
    if "colTau2" in opt_keys
        set_colTau2(cxx_opts, Float64(opts.colTau2))
    end
    if "hessDamp" in opt_keys
        set_hessDamp(cxx_opts, Int32(opts.hessDamp))
    end
    if "hessDampFac" in opt_keys
        set_hessDampFac(cxx_opts, Float64(opts.hessDampFac))
    end
    if "hessUpdate" in opt_keys
        set_hessUpdate(cxx_opts, Int32(opts.hessUpdate))
    end
    if "fallbackUpdate" in opt_keys
        set_fallbackUpdate(cxx_opts, Int32(opts.fallbackUpdate))
    end
    if "hessLimMem" in opt_keys
        set_hessLimMem(cxx_opts, Int32(opts.hessLimMem))
    end
    if "hessMemsize" in opt_keys
        set_hessMemsize(cxx_opts, Int32(opts.hessMemsize))
    end
    if "whichSecondDerv" in opt_keys
        set_whichSecondDerv(cxx_opts, Int32(opts.whichSecondDerv))
    end
    if "skipFirstGlobalization" in opt_keys
        set_skipFirstGlobalization(cxx_opts, Int32(opts.skipFirstGlobalization))
    end
    if "convStrategy" in opt_keys
        set_convStrategy(cxx_opts, Int32(opts.convStrategy))
    end
    if "maxConvQP" in opt_keys
        set_maxConvQP(cxx_opts, Int32(opts.maxConvQP))
    end
    if "maxSOCiter" in opt_keys
        set_maxSOCiter(cxx_opts, Int32(opts.maxSOCiter))
    end

    return cxx_opts

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
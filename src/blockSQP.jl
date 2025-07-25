module blockSQP
    using SparseArrays, Symbolics

	import Base.setproperty!, Base.getproperty
    
    #See Julia documentation of Indirect Calls
    macro dlsym(lib, func)
        z = Ref{Ptr{Cvoid}}(C_NULL)
        quote
            let zlocal = $z[]
                if zlocal == C_NULL
                    zlocal = Base.Libc.Libdl.dlsym($(esc(lib))::Ptr{Cvoid}, $(esc(func)))::Ptr{Cvoid}
                    $z[] = zlocal
                end
                zlocal
            end
        end
    end
    
    const hasjll::Bool =
        try
            import blockSQP_jll
            true
        catch
            false
        end
        
    const libblockSQP = Ref{Ptr{Nothing}}(Ptr{Nothing}())
    function __init__()
        libblockSQP[] = try
            Base.Libc.Libdl.dlopen(joinpath(Base.@__DIR__, "..", "bin", "libblockSQP_jl"))
        catch
            @info "Could not load blockSQP dynamic library from bin folder, loading blockSQP_jll instead\n"
            if !hasjll
                error("Nether local blockSQP dynamic library nor blockSQP_jll are available")
            end
            Base.Libc.Libdl.dlopen(blockSQP_jll.libblockSQP)
        end
    end
    
    
	function fnothing(args...)
	end

	export setindex!
    
    function __lowerbounds(x::AbstractVector)
        return x
    end
    function __upperbounds(x::AbstractVector)
        return x
    end
    function __initial_values(x::AbstractVector)
        return x
    end

    """
    Struct to hold Optimization.jl solver.
    """
    struct BlockSQPOpt end
    export BlockSQPOpt
    
    # Structs to hold structure data used for scaling and condensing
    struct vblock
        size::Integer
        dependent::Bool 
    end

    struct cblock
        size::Integer
    end

    
    include("eval.jl")
    
    include("condenser.jl")
    
    include("problem.jl")

    include("options.jl")
    export blockSQPOptions, qpOASES_options

    include("utils.jl")

    include("solver.jl")
    

end # module blockSQP

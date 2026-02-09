module blockSQP
    using EnumX, SparseArrays
    using Reexport
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
    
    const hasjll::Bool = false        
    const libblockSQP = Ref{Ptr{Nothing}}(Ptr{Nothing}())
    function __init__()
        libblockSQP[] = try
            Base.Libc.Libdl.dlopen(joinpath(Base.@__DIR__, "..", "bin", "libblockSQP_jl"))
        catch blockSQP_load_error
            @info "Could not load blockSQP dynamic library from bin folder." blockSQP_load_error "\nLoading blockSQP_jll instead\n"
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
    struct blockSQPOptimizer end
    Optimizer() = blockSQPOptimizer()
    export blockSQPOptimizer
    
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
    blockSQPoptions = Options
    export blockSQPoptions, qpOASESoptions
    
    include("utils.jl")
    # Some utilities for computing the block structure using Symbolics.
    function compute_hessian_blocks(args...;kwargs...) 
        error("Symbolics.jl being loaded is required for Hessian block computation")
    end

    include("solver.jl")
    
    include("NLPstructures/NLPstructures.jl")
    
    include("structures.jl")
    export create_vblocks, create_condenser_args
    
    #ComponentArrays is used by NLPstructures submodule, so we always use ComponentArrays for now.
    include("ComponentArraysExtension.jl")
end # module blockSQP

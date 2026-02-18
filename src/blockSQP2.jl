module blockSQP2
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
    
    # const hasjll::Bool = try
    import LinearAlgebra
    import OpenBLAS32_jll
    import blockSQP2_jll


    const libblockSQP2 = Ref{Ptr{Nothing}}(Ptr{Nothing}())
    function __init__()
        libblockSQP2[] = try
            throw(ErrorException("."))
            Base.Libc.Libdl.dlopen(joinpath(Base.@__DIR__, "..", "bin", "libblockSQP2_jl"))
        catch blockSQP2_load_error
            @info "Could not load blockSQP dynamic library from bin folder." blockSQP2_load_error "\nTrying blockSQP2_jll instead\n"
            # if !hasjll
            #     error("Nether local blockSQP dynamic library nor blockSQP2_jll are available")
            # end
            LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path; verbose = true)
            Base.Libc.Libdl.dlopen(blockSQP2_jll.libblockSQP2_jl)
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
    struct Optimizer end
    BlockSQP2Optimizer = Optimizer
    export BlockSQP2Optimizer
    
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
    BlockSQP2Problem = Problem
    export BlockSQP2Problem

    include("options.jl")
    BlockSQP2Options = Options
    export BlockSQP2Options
    
    include("utils.jl")
    # Some utilities for computing the block structure using Symbolics.
    function compute_hessian_blocks(args...;kwargs...) 
        error("Symbolics.jl being loaded is required for Hessian block computation")
    end

    include("solver.jl")
    BlockSQP2Stats = Stats
    BlockSQP2Solver = Solver
    export BlockSQP2Stats, BlockSQP2Solver, 
           init!, run!, finish!, get_itCount, 
           get_primal_solution, get_dual_solution, get_dual_solution_full

    include("NLPlayouts/NLPlayouts.jl")
    
    include("layouts.jl")
    export create_vblocks, create_condenser_args
    
    #ComponentArrays is used by NLPlayouts submodule, so we always use ComponentArrays for now.
    include("ComponentArraysExtension.jl")
end # module blockSQP

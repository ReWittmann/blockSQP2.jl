module blockSQP
	import Base.setproperty!, Base.getproperty
    using CxxWrap
    using SparseArrays, Symbolics
    using blockSQP_mumps_jll
    if Sys.iswindows()
        using MKL
    end

	@readmodule(()->libblockSQP_MUMPS_wrapper)
	@wraptypes
	@wrapfunctions

	function __init__()
		@initcxx
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

    include("eval.jl")

    include("problem.jl")

    include("options.jl")
    export BlockSQPOptions

    include("utils.jl")

    include("solver.jl")

end # module blockSQP

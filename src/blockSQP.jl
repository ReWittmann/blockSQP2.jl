module blockSQP
	import Base.setproperty!, Base.getproperty

	using CxxWrap
    using blockSQP_jll
	@readmodule(()->libblockSQP_wrapper)
	@wraptypes
	@wrapfunctions

	function __init__()
		@initcxx
	end

	function fnothing(args...)
	end

	export setindex!


    include("eval.jl")

    include("problem.jl")

    include("options.jl")
    export BlockSQPOptions

    include("utils.jl")

    include("solver.jl")



end # module blockSQP

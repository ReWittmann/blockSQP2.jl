struct Solver
    BSQP_solver::SQPmethod
    Cxx_Problem::Problemform
    BSQP_options::SQPoptions
    BSQP_stats::SQPstats
    Jul_Problem::blockSQPProblem
    Options::BlockSQPOptions
    Solver(J_prob::blockSQPProblem, opts::BlockSQPOptions, cxx_stats::SQPstats) = begin
        #Create problem class on the C++ side
        cxx_prob = Problemform(J_prob.nVar, J_prob.nCon)
        set_scope(cxx_prob, pointer_from_objref(J_prob))
        set_dense_init(cxx_prob, @safe_cfunction(initialize_dense, Nothing, (Ptr{Nothing}, CxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{Float64})))
        set_dense_eval(cxx_prob, @safe_cfunction(evaluate_dense, Nothing, (Ptr{Nothing}, ConstCxxPtr{Float64}, ConstCxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{CxxPtr{Float64}}, Int32, CxxPtr{Int32})))
        set_simple_eval(cxx_prob, @safe_cfunction(evaluate_simple, Nothing, (Ptr{Nothing}, ConstCxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{Int32})))
        set_sparse_init(cxx_prob, @safe_cfunction(initialize_sparse, Nothing, (Ptr{Nothing}, CxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{Int32}, CxxPtr{Int32})))
        set_sparse_eval(cxx_prob, @safe_cfunction(evaluate_sparse, Nothing, (Ptr{Nothing}, ConstCxxPtr{Float64}, ConstCxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{Float64}, CxxPtr{Int32}, CxxPtr{Int32}, CxxPtr{CxxPtr{Float64}}, Int32, CxxPtr{Int32})))
        set_continuity_restoration(cxx_prob, @safe_cfunction(reduceConstrVio, Nothing, (Ptr{Nothing}, CxxPtr{Float64}, CxxPtr{Int32})))
        set_blockIdx(cxx_prob, J_prob.blockIdx)
        set_nnz(cxx_prob, J_prob.nnz)
        set_bounds(cxx_prob, J_prob.lb_var, J_prob.ub_var, J_prob.lb_con, J_prob.ub_con, J_prob.objLo, J_prob.objUp)

        #Create options class on the C++ side
        cxx_opts = set_cxx_options(opts)

        #Create method class on the C++ side
        cxx_method = SQPmethod(CxxPtr(cxx_prob), CxxPtr(cxx_opts), CxxPtr(cxx_stats))

        new(cxx_method, cxx_prob, cxx_opts, cxx_stats, J_prob, opts)
    end
end


function init(meth::Solver)
    cpp_init(meth.BSQP_solver)
end

function run(meth::Solver, maxIt::Integer, warmStart::Integer)
    return cpp_run(meth.BSQP_solver, Int32(maxIt), Int32(warmStart))
end

function finish(meth::Solver)
    cpp_finish(meth.BSQP_solver)
end

function get_primal_solution(meth::Solver)
    xi_ptr = get_primal(meth.BSQP_solver)
    xi_arr = unsafe_wrap(Array{Float64, 1}, xi_ptr.cpp_object, meth.Jul_Problem.nVar, own = false)
    return copy(xi_arr)
end

function get_dual_solution(meth::Solver)
    lam_ptr = get_dual(meth.BSQP_solver)
    lam_arr = unsafe_wrap(Array{Float64, 1}, lam_ptr.cpp_object, meth.Jul_Problem.nVar + meth.Jul_Problem.nCon, own = false)
    return copy(lam_arr[meth.Jul_Problem.nVar+1:end])
end

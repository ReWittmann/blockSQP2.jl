# Holder for C++ side SQPstats object
mutable struct __SQPstats
    obj::Ptr{Cvoid}
    function __SQPstats(arg_obj::Ptr{Cvoid})
        newobj = new(arg_obj)
        function __SQPstats_finalizer!(arg_stats::__SQPstats)
            BSQP = libblockSQP[]
            ccall(@dlsym(BSQP, "delete_SQPstats"), Cvoid, (Ptr{Cvoid},), arg_stats.obj)
        end
        finalizer(__SQPstats_finalizer!, newobj)
    end
end

function SQPstats(outpath::String)
    BSQP = libblockSQP[]
    Ctrans = transcode(UInt8, outpath)
    if !all(Ctrans .<= 0x7f)
        error("SQPstats outpath may not contain non-ASCII characters")
    end
    return __SQPstats(ccall(@dlsym(BSQP, "create_SQPstats"), Ptr{Cvoid}, (Ptr{Cchar},), pointer(reinterpret(Cchar, Ctrans))))
end


@enumx SQPresults::Cint begin
    it_finished = Cint(0)
    partial_success = Cint(1)
    success = Cint(2)
    super_success = Cint(3)
    local_infeasibility = Cint(-1)
    restoration_failure = Cint(-2)
    linesearch_failure = Cint(-3)
    qp_failure = Cint(-4)
    eval_failure = Cint(-5)
    misc_error = Cint(-10)
end

function is_success(ret::SQPresults.T)
    return ret == SQPresults.partial_success || ret == SQPresults.success || ret == SQPresults.super_success
end

mutable struct Solver
    #C++ side objects
    SQPmethod_obj::Ptr{Cvoid}
    Problemspec_obj::Ptr{Cvoid}
    SQPoptions_obj::Ptr{Cvoid}
    QPsolver_options_obj::Ptr{Cvoid}
    
    #Julia side objects
    Jul_Problem::blockSQPProblem
    Jul_Opts::blockSQPOptions
    Jul_Stats::__SQPstats
    
    Solver(J_prob::blockSQPProblem, J_opts::blockSQPOptions, J_stats::__SQPstats) = begin        
        BSQP = libblockSQP[]
        new_Problemspec_obj = ccall(@dlsym(BSQP, "create_Problemspec"), Ptr{Cvoid}, (Cint, Cint), Cint(J_prob.nVar), Cint(J_prob.nCon))
        
        #Shared closure (blockSQPProblem instance) of all callbacks
        ccall(@dlsym(BSQP, "Problemspec_set_closure"), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), new_Problemspec_obj, pointer_from_objref(J_prob))
        
        ccall(@dlsym(BSQP, "Problemspec_set_dense_init"), Cvoid,
            (Ptr{Cvoid}, Ptr{Cvoid}),
            new_Problemspec_obj,
            @cfunction(initialize_dense, Nothing, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble})))

        ccall(@dlsym(BSQP, "Problemspec_set_dense_eval"), Cvoid,
            (Ptr{Cvoid}, Ptr{Cvoid}),
            new_Problemspec_obj,
            @cfunction(evaluate_dense, Nothing,
                        (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble},
                        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                        Ptr{Cdouble}, Ptr{Ptr{Cdouble}}, Cint, Ptr{Cint})))

        ccall(@dlsym(BSQP, "Problemspec_set_simple_eval"), Cvoid,
            (Ptr{Cvoid}, Ptr{Cvoid}),
            new_Problemspec_obj,
            @cfunction(evaluate_simple, Nothing,
                        (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint})))

        ccall(@dlsym(BSQP, "Problemspec_set_sparse_init"), Cvoid,
            (Ptr{Cvoid}, Ptr{Cvoid}),
            new_Problemspec_obj,
            @cfunction(initialize_sparse, Nothing,
                        (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint})))
        
        ccall(@dlsym(BSQP, "Problemspec_set_sparse_eval"), Cvoid,
            (Ptr{Cvoid}, Ptr{Cvoid}),
            new_Problemspec_obj,
            @cfunction(evaluate_sparse, Nothing,
                        (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble},
                        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                        Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint},
                        Ptr{Ptr{Cdouble}}, Cint, Ptr{Cint})))
        
        
        ccall(@dlsym(BSQP, "Problemspec_set_continuity_restoration"), Cvoid,
            (Ptr{Cvoid}, Ptr{Cvoid}),
            new_Problemspec_obj,
            @cfunction(reduceConstrVio, Nothing, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cint})))

        ccall(@dlsym(BSQP, "Problemspec_set_blockIdx"), Cvoid,
            (Ptr{Cvoid}, Ptr{Cint}, Cint), new_Problemspec_obj, pointer(J_prob.blockIdx), Cint(length(J_prob.blockIdx) - 1))
        
        ccall(@dlsym(BSQP, "Problemspec_set_nnz"), Cvoid,
            (Ptr{Cvoid}, Cint), new_Problemspec_obj, Cint(J_prob.nnz))
        
        ccall(@dlsym(BSQP, "Problemspec_set_bounds"), Cvoid,
            (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble},
            Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble),
            new_Problemspec_obj,
            pointer(J_prob.lb_var), pointer(J_prob.ub_var),
            pointer(J_prob.lb_con), pointer(J_prob.ub_con),
            J_prob.lb_obj, J_prob.ub_obj)
        
        
        if length(J_prob.vblocks) > 0
            #Allocate vblocks
            vblock_array_obj = ccall(@dlsym(BSQP, "create_vblock_array"), Ptr{Cvoid}, (Cint, ), Cint(length(J_prob.vblocks)))
            for i = 1:length(J_prob.vblocks)
                ccall(@dlsym(BSQP, "vblock_array_set"), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cchar), vblock_array_obj, Cint(i - 1), Cint(J_prob.vblocks[i].size), Cchar(J_prob.vblocks[i].dependent))
            end
            #Pass ownership of C++ allocated vblocks
            ccall(@dlsym(BSQP, "Problemspec_pass_vblocks"), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Cint), new_Problemspec_obj, vblock_array_obj, Cint(length(J_prob.vblocks)))
        end
        
        if !isnothing(J_prob.cond)
            ccall(@dlsym(BSQP, "Problemspec_set_cond"), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), new_Problemspec_obj, J_prob.cond.Condenser_obj)
        end
        
        #Create blockSQP and QPsolver options classes on the C++ side
        new_SQPoptions_obj, new_QPsolver_options_obj = create_cxx_options(J_opts)
        
        #Create method class on the C++ side. Return nullpointer if an exception is thrown, in which case an error message will be available
        new_SQPmethod_obj = ccall(@dlsym(BSQP, "create_SQPmethod"), Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), new_Problemspec_obj, new_SQPoptions_obj, J_stats.obj)
        if new_SQPmethod_obj == C_NULL
            error(unsafe_string(ccall(@dlsym(BSQP, "get_error_message"), Ptr{Cchar}, ())))
        end
        
        sol = new(new_SQPmethod_obj, new_Problemspec_obj, new_SQPoptions_obj, new_QPsolver_options_obj, J_prob, J_opts, J_stats)
        function Solver_finalizer!(arg_sol::Solver)
            ccall(@dlsym(BSQP, "delete_SQPmethod"), Cvoid, (Ptr{Cvoid},), arg_sol.SQPmethod_obj)
            ccall(@dlsym(BSQP, "delete_SQPoptions"), Cvoid, (Ptr{Cvoid},), arg_sol.SQPoptions_obj)           
            ccall(@dlsym(BSQP, "delete_QPsolver_options"), Cvoid, (Ptr{Cvoid},), arg_sol.QPsolver_options_obj)
            ccall(@dlsym(BSQP, "delete_Problemspec"), Cvoid, (Ptr{Cvoid},), arg_sol.Problemspec_obj)
        end
        finalizer(Solver_finalizer!, sol)
    end
end


function init!(sol::Solver)
    BSQP = libblockSQP[]
    ccall(@dlsym(BSQP, "SQPmethod_init"), Cvoid, (Ptr{Cvoid},), sol.SQPmethod_obj)
end

function run!(sol::Solver, maxIt::Integer, warmStart::Integer)
    BSQP = libblockSQP[]
    ret::Cint = ccall(@dlsym(BSQP, "SQPmethod_run"), Cint, (Ptr{Cvoid}, Cint, Cint), sol.SQPmethod_obj, Cint(maxIt), Cint(warmStart))
    if ret == -1000 #Code for raised exception
        error(unsafe_string(ccall(@dlsym(BSQP, "get_error_message"), Ptr{Cchar}, ())))
    end
    return SQPresults.T(ret)    
end

function finish!(sol::Solver)
    BSQP = libblockSQP[]
    ccall(@dlsym(BSQP, "SQPmethod_finish"), Cvoid, (Ptr{Cvoid},), sol.SQPmethod_obj)
end

function get_itCount(sol::Solver)
    BSQP = libblockSQP[]
    return ccall(@dlsym(BSQP, "SQPstats_get_itCount"), Cint, (Ptr{Cvoid},), sol.SQPstats_obj)
end

#Allocate space for solution on julia side and call C method to fill it
function get_primal_solution(sol::Solver)
    BSQP = libblockSQP[]
    xi_arr = Array{Cdouble, 1}(undef, sol.Jul_Problem.nVar)
    ccall(@dlsym(BSQP, "SQPmethod_get_xi"), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), sol.SQPmethod_obj, pointer(xi_arr))
    return xi_arr
end

function get_dual_solution(sol::Solver)
    BSQP = libblockSQP[]
    lam_arr = Array{Cdouble, 1}(undef, sol.Jul_Problem.nVar + sol.Jul_Problem.nCon)
    ccall(@dlsym(BSQP, "SQPmethod_get_lambda"), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), sol.SQPmethod_obj, pointer(lam_arr))
    return -lam_arr[sol.Jul_Problem.nVar + 1 : end]
end

function get_dual_solution_full(sol::Solver)
    BSQP = libblockSQP[]
    lam_arr = Array{Cdouble, 1}(undef, sol.Jul_Problem.nVar + sol.Jul_Problem.nCon)
    ccall(@dlsym(BSQP, "SQPmethod_get_lambda"), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), sol.SQPmethod_obj, pointer(lam_arr))
    return -lam_arr
end


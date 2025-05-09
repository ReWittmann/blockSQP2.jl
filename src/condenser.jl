
mutable struct Condenser
    cxx_Condenser::Cxx_Condenser
    cxx_vblocks::vblock_array
    cxx_cblocks::cblock_array
    cxx_hsizes::int_array
    cxx_targets::condensing_targets
    Condenser(VBLOCKS::Array{vblock, 1}, CBLOCKS::Array{cblock, 1}, HSIZES::Array{Int32, 1}, TARGETS::Array{condensing_target, 1}, DEP_BOUNDS::Int32 = Int32(2)) = begin
        vblock_arr = vblock_array(Int32(length(VBLOCKS)))
        for i = 1:length(VBLOCKS)
            array_set(vblock_arr, i, VBLOCKS[i])
        end
        
        cblock_arr = cblock_array(Int32(length(CBLOCKS)))
        for i = 1:length(CBLOCKS)
            array_set(cblock_arr, i, CBLOCKS[i])
        end

        hsize_arr = int_array(Int32(length(HSIZES)))
        for i = 1:length(HSIZES)
            array_set(hsize_arr, i, HSIZES[i])
        end

        target_arr = condensing_targets(Int32(length(TARGETS)))
        for i = 1:length(TARGETS)
            array_set(target_arr, i, TARGETS[i])
        end
        
        cond = construct_Condenser(CxxPtr(vblock_arr), CxxRef(cblock_arr), CxxRef(hsize_arr), CxxRef(target_arr), DEP_BOUNDS)
        new(cond, vblock_arr, cblock_arr, hsize_arr, target_arr)
    end

end

function print_info(arg_C::Condenser)
    print_debug(arg_C.cxx_Condenser)
end

function condensed_num_hessblocks(cond::Condenser)
    return get_condensed_num_hessblocks(cond.cxx_Condenser)
end

mutable struct sparse_Matrix
    m::Int32
    n::Int32
    nz::Array{Float64, 1}
    row::Array{Int32, 1}
    colind::Array{Int32, 1}
end
sparse_Matrix() = sparse_Matrix(Int32(0), Int32(0), Float64[], Int32[], Int32[])


function full_condense!(J_cond::Condenser, grad_obj::Array{Float64, 1}, constr_jac::sparse_Matrix, hess::Array{Array{Float64, 2}, 1}, lb_var::Array{Float64, 1}, ub_var::Array{Float64, 1}, lb_con::Array{Float64, 1}, ub_con::Array{Float64, 1})
    #condensed_h::Array{Float64, 1}, condensed_jacobian::sparse_Matrix, condensed_hess::Array{Array{Float64, 2}, 1}, condensed_lb_var::Array{Float64, 1}, condensed_ub_var::Array{Float64, 1}, condensed_lb_con::Array{Float64, 1}, condensed_ub_con::Array{Float64, 1})
    
    cond = J_cond.cxx_Condenser
    nVar = get_num_vars(cond)
    nCon = get_num_cons(cond)
    nnz = length(constr_jac.nz)
    num_hessblocks = get_num_hessblocks(cond)

    M_grad_obj = BSQP_Matrix(nVar, Int32(1))
    grad_obj_data = unsafe_wrap(Array{Float64, 1}, show_ptr(M_grad_obj).cpp_object, nVar, own = false)
    grad_obj_data[:] = grad_obj

    M_constr_jac = alloc_Cxx_Sparse_Matrix(nCon, nVar, nnz)
    constr_nz_data = unsafe_wrap(Array{Float64, 1}, show_nz(M_constr_jac).cpp_object, nnz, own = false)
    constr_row_data = unsafe_wrap(Array{Int32, 1}, show_row(M_constr_jac).cpp_object, nnz, own = false)
    constr_colind_data = unsafe_wrap(Array{Int32, 1}, show_colind(M_constr_jac).cpp_object, nVar+1, own = false)
    constr_nz_data[:] = constr_jac.nz
    constr_row_data[:] = constr_jac.row
    constr_colind_data[:] = constr_jac.colind


    hsize_ptr = get_hess_block_sizes(cond)
    hess_block_sizes = unsafe_wrap(Array{Int32, 1}, hsize_ptr.cpp_object, num_hessblocks, own = false)

    M_hess = SymMat_array(num_hessblocks)
    for i = 1:num_hessblocks
        hsize = hess_block_sizes[i]
        M_hblock = array_get_ptr(M_hess, i)
        set_size!(M_hblock, hsize)
        
        hblock_data = unsafe_wrap(Array{Float64, 1}, show_ptr(M_hblock).cpp_object, Int64(hsize*(hsize + 1)//2), own = false)
        full_to_lower!(reshape(hess[i], Int64(hsize^2)), hblock_data, hsize)
    end

    M_lb_var = BSQP_Matrix(nVar, Int32(1))
    lb_var_data = unsafe_wrap(Array{Float64, 1}, show_ptr(M_lb_var).cpp_object, nVar, own = false)
    lb_var_data[:] = lb_var
    
    M_ub_var = BSQP_Matrix(nVar, Int32(1))
    ub_var_data = unsafe_wrap(Array{Float64, 1}, show_ptr(M_ub_var).cpp_object, nVar, own = false)
    ub_var_data[:] = ub_var

    M_lb_con = BSQP_Matrix(nCon, Int32(1))
    lb_con_data = unsafe_wrap(Array{Float64, 1}, show_ptr(M_lb_con).cpp_object, nCon, own = false)
    lb_con_data[:] = lb_con

    M_ub_con = BSQP_Matrix(nCon, Int32(1))
    ub_con_data = unsafe_wrap(Array{Float64, 1}, show_ptr(M_ub_con).cpp_object, nCon, own = false)
    ub_con_data[:] = ub_con


    #Return arguments, condensed QP
    condensed_num_hessblocks = get_condensed_num_hessblocks(cond)

    M_condensed_h = BSQP_Matrix()
    M_condensed_jacobian = Cxx_Sparse_Matrix()
    M_condensed_hess = SymMat_array(condensed_num_hessblocks)
    M_condensed_lb_var = BSQP_Matrix()
    M_condensed_ub_var = BSQP_Matrix()
    M_condensed_lb_con = BSQP_Matrix()
    M_condensed_ub_con = BSQP_Matrix()

    Cxx_full_condense!(cond, M_grad_obj, M_constr_jac, M_hess, M_lb_var, M_ub_var, M_lb_con, M_ub_con,
        M_condensed_h, M_condensed_jacobian, M_condensed_hess, M_condensed_lb_var, M_condensed_ub_var, M_condensed_lb_con, M_condensed_ub_con)

    condensed_num_vars = get_condensed_num_vars(cond)
    condensed_num_cons = get_condensed_num_cons(cond)
    condensed_nnz = get_nnz(M_condensed_jacobian)

    condensed_h = unsafe_wrap(Array{Float64, 1}, release!(M_condensed_h).cpp_object, condensed_num_vars, own = true)

    c_nz_ptr = show_nz(M_condensed_jacobian)
    c_row_ptr = show_row(M_condensed_jacobian)
    c_colind_ptr = show_colind(M_condensed_jacobian)
    disown!(M_condensed_jacobian)

    condensed_jacobian = sparse_Matrix()
    condensed_jacobian.m = condensed_num_cons
    condensed_jacobian.n = condensed_num_vars
    condensed_jacobian.nz = unsafe_wrap(Array{Float64, 1}, c_nz_ptr.cpp_object, condensed_nnz, own = true)
    condensed_jacobian.row = unsafe_wrap(Array{Int32, 1}, c_row_ptr.cpp_object, condensed_nnz, own = true)
    condensed_jacobian.colind = unsafe_wrap(Array{Int32, 1}, c_colind_ptr.cpp_object, condensed_num_vars + 1, own = true)

    condensed_hess = Array{Array{Float64, 2}, 1}(undef, condensed_num_hessblocks)

    for i = 1:condensed_num_hessblocks
        chblock = array_get_ptr(M_condensed_hess, i)
        b_size = size_1(chblock)
        chblock_data = unsafe_wrap(Array{Float64, 1}, show_ptr(chblock).cpp_object, Int64(b_size*(b_size+1)//2), own = false)
        condensed_hess[i] = Array{Float64, 2}(undef, b_size, b_size)
        lower_to_full!(chblock_data, reshape(condensed_hess[i], Int64(b_size)^2), b_size)
    end

    condensed_lb_var = unsafe_wrap(Array{Float64, 1}, release!(M_condensed_lb_var).cpp_object, condensed_num_vars, own = true)
    condensed_ub_var = unsafe_wrap(Array{Float64, 1}, release!(M_condensed_ub_var).cpp_object, condensed_num_vars, own = true)
    condensed_lb_con = unsafe_wrap(Array{Float64, 1}, release!(M_condensed_lb_con).cpp_object, condensed_num_cons, own = true)
    condensed_ub_con = unsafe_wrap(Array{Float64, 1}, release!(M_condensed_ub_con).cpp_object, condensed_num_cons, own = true)
    
    return condensed_h, condensed_jacobian, condensed_hess, condensed_lb_var, condensed_ub_var, condensed_lb_con, condensed_ub_con

end

function recover_var_mult(J_cond::Condenser, xi_cond::Array{Float64, 1}, lambda_cond::Array{Float64, 1})
    cond = J_cond.cxx_Condenser

    nVar = get_num_vars(cond)
    nCon = get_num_cons(cond)
    condensed_num_vars = get_condensed_num_vars(cond)
    condensed_num_cons = get_condensed_num_cons(cond)

    M_xi_cond = BSQP_Matrix(condensed_num_vars, Int32(1))
    M_lambda_cond = BSQP_Matrix(condensed_num_vars + condensed_num_cons, Int32(1))
    
    xi_cond_data = unsafe_wrap(Array{Float64, 1}, show_ptr(M_xi_cond).cpp_object, condensed_num_vars, own = false)
    xi_cond_data[:] = xi_cond

    lambda_cond_data = unsafe_wrap(Array{Float64, 1}, show_ptr(M_lambda_cond).cpp_object, condensed_num_vars + condensed_num_cons, own = false)
    lambda_cond_data[:] = lambda_cond

    M_xi_rest = BSQP_Matrix()
    M_lambda_rest = BSQP_Matrix()
    Cxx_recover_var_mult!(cond, M_xi_cond, M_lambda_cond, M_xi_rest, M_lambda_rest)

    xi_rest = unsafe_wrap(Array{Float64, 1}, release!(M_xi_rest).cpp_object, nVar, own = true)
    lambda_rest = unsafe_wrap(Array{Float64, 1}, release!(M_lambda_rest).cpp_object, nVar + nCon, own = true)

    return xi_rest, lambda_rest
end






struct condensing_Solver
    #C++ side objects
    BSQP_solver::Cxx_SCQPmethod
    BSQP_problem::Problemform
    BSQP_options::Cxx_SQPoptions
    QPsol_options::Cxx_QPsolver_options
    BSQP_stats::SQPstats
    
    #Julia side objects
    Jul_Problem::blockSQPProblem
    Options::blockSQPOptions
    Jul_Condenser::Condenser
end

condensing_Solver(J_prob::blockSQPProblem, opts::blockSQPOptions, cxx_stats::SQPstats, J_cond::Condenser) = begin
    print("WARNING: This code is untested!\n")
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
    cxx_opts, cxx_QPsol_opts = set_cxx_options(opts)

    #Create method class on the C++ side
    cxx_method = SCQPmethod(CxxPtr(cxx_prob), CxxPtr(cxx_opts), CxxPtr(cxx_stats), CxxPtr(J_cond.cxx_Condenser))

    condensing_Solver(cxx_method, cxx_prob, cxx_opts, cxx_QPsol_opts, cxx_stats, J_prob, opts, J_cond)
end











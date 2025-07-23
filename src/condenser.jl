struct condensing_target
    n_stages::Integer
    first_free::Integer
    vblock_end::Integer
    first_cond::Integer
    cblock_end::Integer
end

mutable struct Condenser
    #Condenser constructor data
    vblock_array_obj::Ptr{Cvoid}
    cblock_array_obj::Ptr{Cvoid}
    hsize_array_obj::Ptr{Cvoid}
    target_array_obj::Ptr{Cvoid}
    
    #Actual C++ side condenser class instance
    Condenser_obj::Ptr{Cvoid}
    
    #C++ side input
    Matrix_grad_obj::Ptr{Cvoid}
    Sparse_Matrix_constr_jac::Ptr{Cvoid}
    SymMatrix_array_hess::Ptr{Cvoid}
    Matrix_lb_var::Ptr{Cvoid}
    Matrix_ub_var::Ptr{Cvoid}
    Matrix_lb_con::Ptr{Cvoid}
    Matrix_ub_con::Ptr{Cvoid}
    
    #C++ side output
    Matrix_condensed_grad_obj::Ptr{Cvoid}
    Sparse_Matrix_condensed_constr_jac::Ptr{Cvoid}
    SymMatrix_array_condensed_hess::Ptr{Cvoid}
    Matrix_condensed_lb_var::Ptr{Cvoid}
    Matrix_condensed_ub_var::Ptr{Cvoid}
    Matrix_condensed_lb_con::Ptr{Cvoid}
    Matrix_condensed_ub_con::Ptr{Cvoid}
    
    #C++ side [cond]ensed QP solution and [rest]ored uncondensed QP solution
    Matrix_xi_cond::Ptr{Cvoid}
    Matrix_lambda_cond::Ptr{Cvoid}
    Matrix_xi_rest::Ptr{Cvoid}
    Matrix_lambda_rest::Ptr{Cvoid}
    
    Condenser(arg_vblocks::Vector{vblock}, arg_cblocks::Vector{cblock}, arg_hsizes::Vector{INT_T}, arg_targets::Vector{condensing_target}, arg_dep_bounds::INT_T = Int32(2)) where INT_T <: Integer = begin
        new_vblock_array_obj = ccall(@dlsym(libblockSQP, "create_vblock_array"), Ptr{Cvoid}, (Cint,), Cint(length(arg_vblocks)))
        for i = 1:length(arg_vblocks)
            ccall(@dlsym(libblockSQP, "vblock_array_set"), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cchar), new_vblock_array_obj, Cint(i - 1), Cint(arg_vblocks[i].size), Cchar(arg_vblocks[i].dependent))
        end
        new_cblock_array_obj = ccall(@dlsym(libblockSQP, "create_cblock_array"), Ptr{Cvoid}, (Cint,), Cint(length(arg_cblocks)))
        for i = 1:length(arg_cblocks)
            ccall(@dlsym(libblockSQP, "cblock_array_set"), Cvoid, (Ptr{Cvoid}, Cint, Cint), new_cblock_array_obj, Cint(i - 1), Cint(arg_cblocks[i].size))
        end
        new_hsize_array_obj = ccall(@dlsym(libblockSQP, "create_hsize_array"), Ptr{Cvoid}, (Cint,), Cint(length(arg_hsizes)))
        for i = 1:length(arg_hsizes)
            ccall(@dlsym(libblockSQP, "hsize_array_set"), Cvoid, (Ptr{Cvoid}, Cint, Cint), new_hsize_array_obj, Cint(i - 1), Cint(arg_hsizes[i]))
        end
        new_target_array_obj = ccall(@dlsym(libblockSQP, "create_target_array"), Ptr{Cvoid}, (Cint,), Cint(length(arg_targets)))
        for i = 1:length(arg_targets)
            ccall(@dlsym(libblockSQP, "target_array_set"), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint, Cint, Cint, Cint), new_target_array_obj, Cint(i - 1), Cint(arg_targets[i].n_stages), Cint(arg_targets[i].first_free), Cint(arg_targets[i].vblock_end), Cint(arg_targets[i].first_cond), Cint(arg_targets[i].cblock_end))
        end
        
        # Pass ownership of vblock_array, cblock_array, hsize_array, target_array
        Condenser_obj = ccall(@dlsym(libblockSQP, "create_Condenser"), Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Ptr{Cvoid}, Cint, Ptr{Cvoid}, Cint, Ptr{Cvoid}, Cint, Cint), new_vblock_array_obj, Cint(length(arg_vblocks)), new_cblock_array_obj, Cint(length(arg_cblocks)), new_hsize_array_obj, Cint(length(arg_hsizes)), new_target_array_obj, Cint(length(arg_targets)), Cint(arg_dep_bounds))
        
        nVar = ccall(@dlsym(libblockSQP, "Condenser_nVar"), Cint, (Ptr{Cvoid},), Condenser_obj)
        nCon = ccall(@dlsym(libblockSQP, "Condenser_nCon"), Cint, (Ptr{Cvoid},), Condenser_obj)
        nBlocks = ccall(@dlsym(libblockSQP, "Condenser_nBlocks"), Cint, (Ptr{Cvoid},), Condenser_obj)
        
        new_grad_obj = ccall(@dlsym(libblockSQP, "create_Matrix"), Ptr{Cvoid}, (Cint, Cint), nVar, Cint(1))
        new_constr_jac = ccall(@dlsym(libblockSQP, "create_Sparse_Matrix_default"), Ptr{Cvoid}, ())
        new_hess = ccall(@dlsym(libblockSQP, "create_SymMatrix_array"), Ptr{Cvoid}, (Cint,), nBlocks)
        for i = Cint(1):nBlocks
            ccall(@dlsym(libblockSQP, "SymMatrix_array_index_resize"), Cvoid, (Ptr{Cvoid}, Cint, Cint), new_hess, Cint(i - 1), Cint(arg_hsizes[i]))
        end
        new_lb_var = ccall(@dlsym(libblockSQP, "create_Matrix"), Ptr{Cvoid}, (Cint, Cint), nVar, Cint(1));
        new_ub_var = ccall(@dlsym(libblockSQP, "create_Matrix"), Ptr{Cvoid}, (Cint, Cint), nVar, Cint(1));
        new_lb_con = ccall(@dlsym(libblockSQP, "create_Matrix"), Ptr{Cvoid}, (Cint, Cint), nCon, Cint(1));
        new_ub_con = ccall(@dlsym(libblockSQP, "create_Matrix"), Ptr{Cvoid}, (Cint, Cint), nCon, Cint(1));
        
        condensed_nVar = ccall(@dlsym(libblockSQP, "Condenser_condensed_nVar"), Cint, (Ptr{Cvoid},), Condenser_obj)
        condensed_nCon = ccall(@dlsym(libblockSQP, "Condenser_condensed_nCon"), Cint, (Ptr{Cvoid},), Condenser_obj)
        condensed_nBlocks = ccall(@dlsym(libblockSQP, "Condenser_condensed_nBlocks"), Cint, (Ptr{Cvoid},), Condenser_obj)
        
        new_condensed_grad_obj = ccall(@dlsym(libblockSQP, "create_Matrix"), Ptr{Cvoid}, (Cint, Cint), condensed_nVar, Cint(1))
        new_condensed_constr_jac = ccall(@dlsym(libblockSQP, "create_Sparse_Matrix_default"), Ptr{Cvoid}, ())
        new_condensed_hess = ccall(@dlsym(libblockSQP, "create_SymMatrix_array"), Ptr{Cvoid}, (Cint,), condensed_nBlocks)
        new_condensed_lb_var = ccall(@dlsym(libblockSQP, "create_Matrix_default"), Ptr{Cvoid}, ());
        new_condensed_ub_var = ccall(@dlsym(libblockSQP, "create_Matrix_default"), Ptr{Cvoid}, ());
        new_condensed_lb_con = ccall(@dlsym(libblockSQP, "create_Matrix_default"), Ptr{Cvoid}, ());
        new_condensed_ub_con = ccall(@dlsym(libblockSQP, "create_Matrix_default"), Ptr{Cvoid}, ());
        
        new_xi_cond = ccall(@dlsym(libblockSQP, "create_Matrix"), Ptr{Cvoid}, (Cint, Cint), condensed_nVar, Cint(1))
        new_lambda_cond = ccall(@dlsym(libblockSQP, "create_Matrix"), Ptr{Cvoid}, (Cint, Cint), condensed_nVar + condensed_nCon, Cint(1))
        
        new_xi_rest = ccall(@dlsym(libblockSQP, "create_Matrix"), Ptr{Cvoid}, (Cint, Cint), nVar, Cint(1))
        new_lambda_rest = ccall(@dlsym(libblockSQP, "create_Matrix"), Ptr{Cvoid}, (Cint, Cint), nVar + nCon, Cint(1))
        
        new_Condenser = new(new_vblock_array_obj, new_cblock_array_obj, new_hsize_array_obj, new_target_array_obj,
                            Condenser_obj,
                            new_grad_obj, new_constr_jac, new_hess, new_lb_var, new_ub_var, new_lb_con, new_ub_con,
                            new_condensed_grad_obj, new_condensed_constr_jac, new_condensed_hess, new_condensed_lb_var, new_condensed_ub_var, new_condensed_lb_con, new_condensed_ub_con,
                            new_xi_cond, new_lambda_cond,
                            new_xi_rest, new_lambda_rest
                            )
        
        function Condenser_finalizer(J_cond::Condenser)
            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_lambda_rest)
            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_xi_rest)
            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_lambda_cond)
            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_xi_cond)

            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_condensed_ub_con)
            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_condensed_lb_con)
            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_condensed_ub_var)
            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_condensed_lb_var)
            ccall(@dlsym(libblockSQP, "delete_SymMatrix_array"), Cvoid, (Ptr{Cvoid},), J_cond.SymMatrix_array_condensed_hess)
            ccall(@dlsym(libblockSQP, "delete_Sparse_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Sparse_Matrix_condensed_constr_jac)
            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_condensed_grad_obj)

            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_ub_con)
            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_lb_con)
            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_ub_var)
            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_lb_var)
            ccall(@dlsym(libblockSQP, "delete_SymMatrix_array"), Cvoid, (Ptr{Cvoid},), J_cond.SymMatrix_array_hess)
            ccall(@dlsym(libblockSQP, "delete_Sparse_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Sparse_Matrix_constr_jac)
            ccall(@dlsym(libblockSQP, "delete_Matrix"), Cvoid, (Ptr{Cvoid},), J_cond.Matrix_grad_obj)

            ccall(@dlsym(libblockSQP, "delete_Condenser"), Cvoid, (Ptr{Cvoid},), J_cond.Condenser_obj)

            ccall(@dlsym(libblockSQP, "delete_target_array"), Cvoid, (Ptr{Cvoid},), J_cond.target_array_obj)
            ccall(@dlsym(libblockSQP, "delete_hsize_array"), Cvoid, (Ptr{Cvoid},), J_cond.hsize_array_obj)
            ccall(@dlsym(libblockSQP, "delete_cblock_array"), Cvoid, (Ptr{Cvoid},), J_cond.cblock_array_obj)
            ccall(@dlsym(libblockSQP, "delete_vblock_array"), Cvoid, (Ptr{Cvoid},), J_cond.vblock_array_obj)
        end
        finalizer(Condenser_finalizer, new_Condenser)
    end
end

function print_info(J_cond::Condenser)
    ccall(@dlsym(libblockSQP, "Condenser_print_debug"), Cvoid, (Ptr{Cvoid},), J_cond.Condenser_obj)
end

function condensed_nBlocks(J_cond::Condenser)
    return ccall(@dlsym(libblockSQP, "Condenser_condensed_nBlocks"), Cint, (Ptr{Cvoid},), J_cond.Condenser_obj)
end


mutable struct Sparse_Matrix
    m::Int32
    n::Int32
    nz::Vector{Cdouble}
    row::Vector{Cint}
    colind::Vector{Cint}
    Sparse_Matrix(arg_m::Integer, arg_n::Integer, arg_nz::Vector{FLOAT_T}, arg_row::Vector{INT_T}, arg_colind::Vector{INT_T}) where {FLOAT_T <: AbstractFloat, INT_T <: Integer} = begin
        new(arg_m, arg_n, Cdouble[Cdouble(v) for v in arg_nz], Cint[Cint(v) for v in arg_row], Cint[Cint(v) for v in arg_colind])
    end
end

Sparse_Matrix() = Sparse_Matrix(Cint(0), Cint(0), Cdouble[], Cint[], Cint[])
Sparse_Matrix(arg_m::Integer, arg_n::Integer, arg_nnz::Integer) = begin
    colind = Vector{Cint}(undef, arg_n + 1)
    colind[arg_n + 1] = Cint(arg_nnz)
    return Sparse_Matrix(Cint(arg_m), Cint(arg_n), Vector{Cdouble}(undef, arg_nnz), Vector{Cint}(undef, arg_nnz), colind)
end


function full_condense!(J_cond::Condenser, grad_obj::Vector{Float64}, constr_jac::Sparse_Matrix, hess::Vector{Matrix{Float64}}, lb_var::Vector{Float64}, ub_var::Vector{Float64}, lb_con::Vector{Float64}, ub_con::Vector{Float64})
    cond = J_cond.Condenser_obj
    
    nVar = ccall(@dlsym(libblockSQP, "Condenser_nVar"), Cint, (Ptr{Cvoid},), cond)
    nCon = ccall(@dlsym(libblockSQP, "Condenser_nCon"), Cint, (Ptr{Cvoid},), cond)
    nBlocks = ccall(@dlsym(libblockSQP, "Condenser_nBlocks"), Cint, (Ptr{Cvoid},), cond)
    
    nnz = Cint(length(constr_jac.nz))
    @assert nnz == constr_jac.colind[constr_jac.n + 1]
    
    unsafe_wrap(Vector{Cdouble}, 
                ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_grad_obj), 
                Cint(nVar);
                own = false
                )[:] = grad_obj
    
    ccall(@dlsym(libblockSQP, "Sparse_Matrix_set_structure"), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint), J_cond.Sparse_Matrix_constr_jac, nCon, nVar, nnz)
    
    unsafe_wrap(Vector{Cdouble}, 
                ccall(@dlsym(libblockSQP, "Sparse_Matrix_nz"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Sparse_Matrix_constr_jac), 
                nnz; 
                own = false
                )[:] = constr_jac.nz
    
    unsafe_wrap(Vector{Cint}, 
                ccall(@dlsym(libblockSQP, "Sparse_Matrix_row"), Ptr{Cint}, (Ptr{Cvoid},), J_cond.Sparse_Matrix_constr_jac), 
                nnz; 
                own = false
                )[:] = constr_jac.row
    
    unsafe_wrap(Vector{Cint}, 
                ccall(@dlsym(libblockSQP, "Sparse_Matrix_colind"), Ptr{Cint}, (Ptr{Cvoid},), J_cond.Sparse_Matrix_constr_jac), 
                Cint(constr_jac.n + 1);
                own = false
                )[:] = constr_jac.colind
    
    nBlocks = ccall(@dlsym(libblockSQP, "Condenser_nBlocks"), Cint, (Ptr{Cvoid},), cond)
    hsizes = unsafe_wrap(Vector{Cint}, 
                         ccall(@dlsym(libblockSQP, "Condenser_hsizes"), Ptr{Cint}, (Ptr{Cvoid},), cond), 
                         nBlocks; 
                         own = false)
    for i = Cint(1):nBlocks
        hsize = hsizes[i]
        hessblock_data = unsafe_wrap(Vector{Cdouble}, 
                                     ccall(@dlsym(libblockSQP, "SymMatrix_array_index_array"), Ptr{Cdouble}, (Ptr{Cvoid}, Cint), J_cond.SymMatrix_array_hess, Cint(i - 1)), 
                                     Cint((hsize*(hsize + 1))//2);
                                     own = false)
        full_to_lower!(hessblock_data, reshape(hess[i], Int64(hsize^2)), hsize)
    end
    
    unsafe_wrap(Vector{Cdouble},
                ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_lb_var),
                nVar;
                own = false
                )[:] = lb_var
    
    unsafe_wrap(Vector{Cdouble},
                ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_ub_var),
                nVar;
                own = false
                )[:] = ub_var
    
    unsafe_wrap(Vector{Cdouble},
                ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_lb_con),
                nCon;
                own = false
                )[:] = lb_con
    
    unsafe_wrap(Vector{Cdouble},
                ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_ub_con),
                nCon;
                own = false
                )[:] = ub_con
    
    ccall(@dlsym(libblockSQP, "Condenser_full_condense"), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
                cond, 
                J_cond.Matrix_grad_obj, J_cond.Sparse_Matrix_constr_jac, J_cond.SymMatrix_array_hess, J_cond.Matrix_lb_var, J_cond.Matrix_ub_var, J_cond.Matrix_lb_con, J_cond.Matrix_ub_con,
                J_cond.Matrix_condensed_grad_obj, J_cond.Sparse_Matrix_condensed_constr_jac, J_cond.SymMatrix_array_condensed_hess, J_cond.Matrix_condensed_lb_var, J_cond.Matrix_condensed_ub_var, J_cond.Matrix_condensed_lb_con, J_cond.Matrix_condensed_ub_con
          )
    
    condensed_nVar = ccall(@dlsym(libblockSQP, "Condenser_condensed_nVar"), Cint, (Ptr{Cvoid},), cond)
    condensed_nCon = ccall(@dlsym(libblockSQP, "Condenser_condensed_nCon"), Cint, (Ptr{Cvoid},), cond)
    condensed_nnz = ccall(@dlsym(libblockSQP, "Sparse_Matrix_nnz"), Cint, (Ptr{Cvoid},), J_cond.Sparse_Matrix_condensed_constr_jac)
    
    condensed_grad_obj = Vector{Cdouble}(undef, condensed_nVar)
    condensed_grad_obj[:] = unsafe_wrap(Vector{Cdouble},
                                        ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_condensed_grad_obj),
                                        condensed_nVar;
                                        own = false
                                        )
    
    condensed_constr_jac = Sparse_Matrix(condensed_nCon, condensed_nVar, condensed_nnz)
    condensed_constr_jac.nz[:] = unsafe_wrap(Vector{Cdouble},
                                             ccall(@dlsym(libblockSQP, "Sparse_Matrix_nz"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Sparse_Matrix_condensed_constr_jac),
                                             condensed_nnz;
                                             own = false
                                             )
    condensed_constr_jac.row[:] = unsafe_wrap(Vector{Cint},
                                              ccall(@dlsym(libblockSQP, "Sparse_Matrix_row"), Ptr{Cint}, (Ptr{Cvoid},), J_cond.Sparse_Matrix_condensed_constr_jac),
                                              condensed_nnz;
                                              own = false
                                              )
    condensed_constr_jac.colind[:] = unsafe_wrap(Vector{Cint},
                                                 ccall(@dlsym(libblockSQP, "Sparse_Matrix_colind"), Ptr{Cint}, (Ptr{Cvoid},), J_cond.Sparse_Matrix_condensed_constr_jac),
                                                 condensed_nVar + 1;
                                                 own = false
                                                 )
    
    condensed_nBlocks = ccall(@dlsym(libblockSQP, "Condenser_condensed_nBlocks"), Cint, (Ptr{Cvoid},), cond)
    condensed_hsizes = unsafe_wrap(Vector{Cint}, 
                                   ccall(@dlsym(libblockSQP, "Condenser_condensed_hsizes"), Ptr{Cint}, (Ptr{Cvoid},), cond),
                                   condensed_nBlocks;
                                   own = false)
    
    condensed_hess = Vector{Matrix{Cdouble}}(undef, condensed_nBlocks)
    for i = 1:condensed_nBlocks
        condensed_hsize = condensed_hsizes[i]
        condensed_hess[i] = Matrix{Float64}(undef, condensed_hsize, condensed_hsize)
        condensed_hessblock_data = unsafe_wrap(Vector{Cdouble},
                                               ccall(@dlsym(libblockSQP, "SymMatrix_array_index_array"), Ptr{Cdouble}, (Ptr{Cvoid}, Cint), J_cond.SymMatrix_array_condensed_hess, Cint(i - 1)),
                                               condensed_hsize^2;
                                               own = false
                                               )
        lower_to_full!(reshape(condensed_hess[i], Int64(condensed_hsize^2)), condensed_hessblock_data, condensed_hsize)
    end
    
    condensed_lb_var = Vector{Cdouble}(undef, condensed_nVar)
    condensed_ub_var = Vector{Cdouble}(undef, condensed_nVar)
    condensed_lb_con = Vector{Cdouble}(undef, condensed_nCon)
    condensed_ub_con = Vector{Cdouble}(undef, condensed_nCon)
    condensed_lb_var[:] = unsafe_wrap(Vector{Cdouble},
                                      ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_condensed_lb_var),
                                      condensed_nVar;
                                      own = false
                                      )
    condensed_ub_var[:] = unsafe_wrap(Vector{Cdouble},
                                      ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_condensed_ub_var),
                                      condensed_nVar;
                                      own = false
                                      )
    condensed_lb_con[:] = unsafe_wrap(Vector{Cdouble},
                                      ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_condensed_lb_con),
                                      condensed_nCon;
                                      own = false
                                      )
    condensed_ub_con[:] = unsafe_wrap(Vector{Cdouble},
                                      ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_condensed_ub_con),
                                      condensed_nCon;
                                      own = false
                                      )
    return condensed_grad_obj, condensed_constr_jac, condensed_hess, condensed_lb_var, condensed_ub_var, condensed_lb_con, condensed_ub_con
end

function recover_var_mult(J_cond::Condenser, xi_cond::Array{Float64, 1}, lambda_cond::Array{Float64, 1})
    cond = J_cond.Condenser_obj
    
    nVar = ccall(@dlsym(libblockSQP, "Condenser_nVar"), Cint, (Ptr{Cvoid},), cond)
    nCon = ccall(@dlsym(libblockSQP, "Condenser_nCon"), Cint, (Ptr{Cvoid},), cond)
    condensed_nVar = ccall(@dlsym(libblockSQP, "Condenser_condensed_nVar"), Cint, (Ptr{Cvoid},), cond)
    condensed_nCon = ccall(@dlsym(libblockSQP, "Condenser_condensed_nCon"), Cint, (Ptr{Cvoid},), cond)
    
    unsafe_wrap(Vector{Cdouble},
                ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_xi_cond),
                condensed_nVar;
                own = false
                )[:] = xi_cond
    unsafe_wrap(Vector{Cdouble},
                ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_lambda_cond),
                condensed_nVar + condensed_nCon;
                own = false
                )[:] = lambda_cond
    
    ccall(@dlsym(libblockSQP, "Condenser_recover_var_mult"), Cvoid, 
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), 
        cond, J_cond.Matrix_xi_cond, J_cond.Matrix_lambda_cond, J_cond.Matrix_xi_rest, J_cond.Matrix_lambda_rest)
    
    xi_rest = Vector{Cdouble}(undef, nVar)
    lambda_rest = Vector{Cdouble}(undef, nVar + nCon)
    xi_rest[:] = unsafe_wrap(Vector{Cdouble},
                             ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_xi_rest),
                             nVar;
                             own = false
                             )
    lambda_rest[:] = unsafe_wrap(Vector{Cdouble},
                                 ccall(@dlsym(libblockSQP, "Matrix_array"), Ptr{Cdouble}, (Ptr{Cvoid},), J_cond.Matrix_lambda_rest),
                                 nVar + nCon;
                                 own = false
                                 )
    return xi_rest, lambda_rest
end



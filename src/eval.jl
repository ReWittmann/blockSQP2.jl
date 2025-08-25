function initialize_dense(Prob::Ptr{Nothing}, xi::Ptr{Cdouble}, lam::Ptr{Cdouble}, Jac::Ptr{Cdouble})
    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    xi_arr = unsafe_wrap(Array{Cdouble, 1}, xi, Jprob.nVar, own = false)
    lam_arr = unsafe_wrap(Array{Cdouble, 1}, lam, Jprob.nVar + Jprob.nCon, own = false)

    xi_arr[:] = Jprob.x0
    lam_arr[:] = Jprob.lambda0
    return
end

function evaluate_dense(Prob::Ptr{Nothing}, xi::Ptr{Cdouble}, lam::Ptr{Cdouble}, objval::Ptr{Cdouble}, constr::Ptr{Cdouble}, gradObj::Ptr{Cdouble}, constrJac::Ptr{Cdouble}, hess::Ptr{Ptr{Cdouble}}, dmode::Cint, info::Ptr{Cint})
    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    xi_arr = unsafe_wrap(Array{Cdouble, 1}, xi, Jprob.nVar, own = false)
    lam_arr = unsafe_wrap(Array{Cdouble, 1}, lam, Jprob.nVar + Jprob.nCon, own = false)
    constr_arr = unsafe_wrap(Array{Cdouble, 1}, constr, Jprob.nCon, own = false)
    
    unsafe_store!(objval, Cdouble(Jprob.f(xi_arr)))
    constr_arr[:] .= Jprob.g(xi_arr)

    if dmode > 0
        gradObj_arr = unsafe_wrap(Array{Cdouble, 1}, gradObj, Jprob.nVar, own = false)
        constrJac_arr = unsafe_wrap(Array{Cdouble, 2}, constrJac, (Jprob.nCon, Jprob.nVar), own = false)
        gradObj_arr[:] = Jprob.grad_f(xi_arr)
        constrJac_arr[:,:] .= Jprob.jac_g(xi_arr)
        if dmode == 2
            hess_arr = unsafe_wrap(Array{CxxPtr{Cdouble}, 1}, hess, Jprob.n_hessblocks, own = false)
            
            s = Jprob.blockIdx[Jprob.n_hessblocks + 1] - Jprob.blockIdx[Jprob.n_hessblocks]
            hess_last = unsafe_wrap(Array{Cdouble,1}, hess_arr[Jprob.n_hessblocks], Cint((s*(s+Cint(1)))//(Cint(2))), own = false)
            hess_last[:] = Jprob.last_hessBlock(xi_arr, lam_arr[Jprob.nVar + 1 : Jprob.nVar + Jprob.nCon])
        elseif dmode == 3
            hessPTR_arr = unsafe_wrap(Array{CxxPtr{Cdouble}, 1}, hess, Jprob.n_hessblocks, own = false)
            hess_arr = Array{Array{Cdouble, 1}, 1}(undef, Jprob.n_hessblocks)
            for i = 1:Jprob.n_hessblocks
                Bsize = Jprob.blockIdx[i+1] - Jprob.blockIdx[i]
                hess_arr[i] = unsafe_wrap(Array{Cdouble,1}, hessPTR_arr[i], Cint((Bsize*(Bsize + Cint(1)))//Cint(2)), own = false)
            end
            
            hess_eval = Jprob.hess(xi_arr, lam_arr[Jprob.nVar + 1 : Jprob.nVar + Jprob.nCon])
            for i = 1:Jprob.n_hessblocks
                hess_arr[i][:] = hess_eval[i]
            end
        end
    end
    unsafe_store!(info, Cint(0))
    return
end

function evaluate_simple(Prob::Ptr{Nothing}, xi::Ptr{Cdouble}, objval::Ptr{Cdouble}, constr::Ptr{Cdouble}, info::Ptr{Cint})
    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    xi_arr = unsafe_wrap(Array{Cdouble, 1}, xi, Jprob.nVar, own = false)
    constr_arr = unsafe_wrap(Array{Cdouble, 1}, constr, Jprob.nCon, own = false)
    
    unsafe_store!(objval, Cdouble(Jprob.f(xi_arr)))
    constr_arr[:] .= Jprob.g(xi_arr)
    
    unsafe_store!(info, Cint(0))
    return
end


function initialize_sparse(Prob::Ptr{Nothing}, xi::Ptr{Cdouble}, lam::Ptr{Cdouble}, jac_nz::Ptr{Cdouble}, jac_row::Ptr{Cint}, jac_colind::Ptr{Cint})
    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    xi_arr = unsafe_wrap(Array{Cdouble, 1}, xi, Jprob.nVar, own = false)
    lam_arr = unsafe_wrap(Array{Cdouble, 1}, lam, Jprob.nVar + Jprob.nCon, own = false)
    jac_row_arr = unsafe_wrap(Array{Cint, 1}, jac_row, Jprob.nnz, own = false)
    jac_colind_arr = unsafe_wrap(Array{Cint, 1}, jac_colind, Jprob.nVar + 1, own = false)

    xi_arr[:] = Jprob.x0
    lam_arr[:] = Jprob.lambda0
    jac_row_arr[:] = Jprob.jac_g_row
    jac_colind_arr[:] = Jprob.jac_g_colind
    return
end


function evaluate_sparse(Prob::Ptr{Nothing}, xi::Ptr{Cdouble}, lam::Ptr{Cdouble}, objval::Ptr{Cdouble}, constr::Ptr{Cdouble}, gradObj::Ptr{Cdouble}, jac_nz::Ptr{Cdouble}, jac_row::Ptr{Cint}, jac_colind::Ptr{Cint}, hess::Ptr{Ptr{Cdouble}}, dmode::Cint, info::Ptr{Cint})
    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    xi_arr = unsafe_wrap(Array{Cdouble, 1}, xi, Jprob.nVar, own = false)
    lam_arr = unsafe_wrap(Array{Cdouble, 1}, lam, Jprob.nVar + Jprob.nCon, own = false)
    constr_arr = unsafe_wrap(Array{Cdouble, 1}, constr, Jprob.nCon, own = false)
    jac_nz_arr = unsafe_wrap(Array{Cdouble, 1}, jac_nz, Jprob.nnz, own = false)

    unsafe_store!(objval, Cdouble(Jprob.f(xi_arr)))
    constr_arr[:] .= Jprob.g(xi_arr)

    if dmode > 0
        gradObj_arr = unsafe_wrap(Array{Cdouble, 1}, gradObj, Jprob.nVar, own = false)
        gradObj_arr[:] = Jprob.grad_f(xi_arr)
        # jac_nz_arr[:] = Jprob.jac_g_nz(xi_arr)
        
        jac_g_nz_eval = Jprob.jac_g_nz(xi_arr)
        @assert length(jac_g_nz_eval) == Jprob.nnz
        jac_nz_arr[:] = jac_g_nz_eval
        if dmode == 2
            hess_arr = unsafe_wrap(Array{CxxPtr{Cdouble}, 1}, hess, Jprob.n_hessblocks, own = false)

            s_last = Jprob.blockIdx[Jprob.n_hessblocks + 1] - Jprob.blockIdx[Jprob.n_hessblocks]
            hess_last = unsafe_wrap(Array{Cdouble,1}, hess_arr[Jprob.n_hessblocks], Cint((s_last*(s_last + Cint(1)))//(Cint(2))), own = false)
            hess_last[:] = Jprob.last_hessBlock(xi_arr, lam_arr[Jprob.nVar + 1 : Jprob.nVar + Jprob.nCon])
        end

        if dmode == 3
            hessPTR_arr = unsafe_wrap(Array{CxxPtr{Cdouble}, 1}, hess, Jprob.n_hessblocks, own = false)
            hess_arr = Array{Array{Cdouble, 1}, 1}(undef, Jprob.n_hessblocks)
            for i = 1:Jprob.n_hessblocks
                Bsize = Jprob.blockIdx[i+1] - Jprob.blockIdx[i]
                hess_arr[i] = unsafe_wrap(Array{Cdouble,1}, hessPTR_arr[i], Cint((Bsize*(Bsize + Cint(1)))//Cint(2)), own = false)
            end

            hess_eval = Jprob.hess(xi_arr, lam_arr[Jprob.nVar + 1 : Jprob.nVar + Jprob.nCon])
            for i = 1:Jprob.n_hessblocks
                hess_arr[i][:] = hess_eval[i]
            end
        end
    end
    unsafe_store!(info, Cint(0))
    return
end


function reduceConstrVio(Prob::Ptr{Cvoid}, xi::Ptr{Cdouble}, info::Ptr{Cint})
    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    if Jprob.continuity_restoration == fnothing
        unsafe_store!(info, Cint(1))
    else
        xi_arr = unsafe_wrap(Array{Cdouble, 1}, xi.cpp_object, Jprob.nVar, own = false)
        xi_arr[:] = Jprob.continuity_restoration(xi_arr)
        unsafe_store!(info, Cint(0))
    end
    return
end

function initialize_dense(Prob::Ptr{Nothing}, xi::CxxPtr{Float64}, lam::CxxPtr{Float64}, Jac::CxxPtr{Float64})
    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    xi_arr = unsafe_wrap(Array{Float64, 1}, xi.cpp_object, Jprob.nVar, own = false)
    lam_arr = unsafe_wrap(Array{Float64, 1}, lam.cpp_object, Jprob.nVar + Jprob.nCon, own = false)

    xi_arr[:] = Jprob.x0
    lam_arr[:] = Jprob.lambda0
    return
end

function evaluate_dense(Prob::Ptr{Nothing}, xi::ConstCxxPtr{Float64}, lam::ConstCxxPtr{Float64}, objval::CxxPtr{Float64}, constr::CxxPtr{Float64}, gradObj::CxxPtr{Float64}, constrJac::CxxPtr{Float64}, hess::CxxPtr{CxxPtr{Float64}}, dmode::Int32, info::CxxPtr{Int32})

    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    xi_arr = unsafe_wrap(Array{Float64, 1}, xi.cpp_object, Jprob.nVar, own = false)
    lam_arr = unsafe_wrap(Array{Float64, 1}, lam.cpp_object, Jprob.nVar + Jprob.nCon, own = false)
    constr_arr = unsafe_wrap(Array{Float64, 1}, constr.cpp_object, Jprob.nCon, own = false)

    objval[] = Jprob.f(xi_arr)
    constr_arr[:] .= Jprob.g(xi_arr)

    if dmode > 0
        gradObj_arr = unsafe_wrap(Array{Float64, 1}, gradObj.cpp_object, Jprob.nVar, own = false)
        constrJac_arr = unsafe_wrap(Array{Float64, 2}, constrJac.cpp_object, (Jprob.nCon, Jprob.nVar), own = false)
        gradObj_arr[:] = Jprob.grad_f(xi_arr)
        constrJac_arr[:,:] .= Jprob.jac_g(xi_arr)

        if dmode == 2
            hess_arr = unsafe_wrap(Array{CxxPtr{Float64}, 1}, hess.cpp_object, Jprob.n_hessblocks, own = true)

            s = Jprob.blockIdx[Jprob.n_hessblocks + 1] - Jprob.blockIdx[Jprob.n_hessblocks]
            hess_last = unsafe_wrap(Array{Float64,1}, hess_arr[Jprob.n_hessblocks].cpp_object, Int32((s*(s+Int32(1)))//(Int32(2))), own = false)
            hess_last[:] = Jprob.last_hessBlock(xi_arr, lam_arr[Jprob.nVar + 1 : Jprob.nVar + Jprob.nCon])
        elseif dmode == 3
            hessPTR_arr = unsafe_wrap(Array{CxxPtr{Float64}, 1}, hess.cpp_object, Jprob.n_hessblocks, own = true)
            hess_arr = Array{Array{Float64, 1}, 1}(undef, Jprob.n_hessblocks)
            for i = 1:Jprob.n_hessblocks
                Bsize = Jprob.blockIdx[i+1] - Jprob.blockIdx[i]
                hess_arr[i] = unsafe_wrap(Array{Float64,1}, hessPTR_arr[i].cpp_object, Int32((Bsize*(Bsize + Int32(1)))//Int32(2)), own = false)
            end

            hess_eval = Jprob.hess(xi_arr, lam_arr[Jprob.nVar + 1 : Jprob.nVar + Jprob.nCon])
            for i = 1:Jprob.n_hessblocks
                hess_arr[i][:] = hess_eval[i]
            end
        end
    end
    info[] = Int32(0);
    return
end

function evaluate_simple(Prob::Ptr{Nothing}, xi::ConstCxxPtr{Float64}, objval::CxxPtr{Float64}, constr::CxxPtr{Float64}, info::CxxPtr{Int32})
    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    xi_arr = unsafe_wrap(Array{Float64, 1}, xi.cpp_object, Jprob.nVar, own = false)
    constr_arr = unsafe_wrap(Array{Float64, 1}, constr.cpp_object, Jprob.nCon, own = false)

    objval[] = Jprob.f(xi_arr)
    constr_arr[:] .= Jprob.g(xi_arr)

    info[] = Int32(0);
    return
end


function initialize_sparse(Prob::Ptr{Nothing}, xi::CxxPtr{Float64}, lam::CxxPtr{Float64}, jac_nz::CxxPtr{Float64}, jac_row::CxxPtr{Int32}, jac_colind::CxxPtr{Int32})
    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    xi_arr = unsafe_wrap(Array{Float64, 1}, xi.cpp_object, Jprob.nVar, own = false)
    lam_arr = unsafe_wrap(Array{Float64, 1}, lam.cpp_object, Jprob.nVar + Jprob.nCon, own = false)
    jac_row_arr = unsafe_wrap(Array{Int32, 1}, jac_row.cpp_object, Jprob.nnz, own = false)
    jac_colind_arr = unsafe_wrap(Array{Int32, 1}, jac_colind.cpp_object, Jprob.nVar + 1, own = false)

    xi_arr[:] = Jprob.x0
    lam_arr[:] = Jprob.lambda0
    jac_row_arr[:] = Jprob.jac_g_row
    jac_colind_arr[:] = Jprob.jac_g_colind
    return
end


function evaluate_sparse(Prob::Ptr{Nothing}, xi::ConstCxxPtr{Float64}, lam::ConstCxxPtr{Float64}, objval::CxxPtr{Float64}, constr::CxxPtr{Float64}, gradObj::CxxPtr{Float64}, jac_nz::CxxPtr{Float64}, jac_row::CxxPtr{Int32}, jac_colind::CxxPtr{Int32}, hess::CxxPtr{CxxPtr{Float64}}, dmode::Int32, info::CxxPtr{Int32})
    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    xi_arr = unsafe_wrap(Array{Float64, 1}, xi.cpp_object, Jprob.nVar, own = false)
    lam_arr = unsafe_wrap(Array{Float64, 1}, lam.cpp_object, Jprob.nVar + Jprob.nCon, own = false)
    constr_arr = unsafe_wrap(Array{Float64, 1}, constr.cpp_object, Jprob.nCon, own = false)
    jac_nz_arr = unsafe_wrap(Array{Float64, 1}, jac_nz.cpp_object, Jprob.nnz, own = false)

    objval[] = Jprob.f(xi_arr)
    constr_arr[:] = Jprob.g(xi_arr)

    if dmode > 0
        gradObj_arr = unsafe_wrap(Array{Float64, 1}, gradObj.cpp_object, Jprob.nVar, own = false)
        gradObj_arr[:] = Jprob.grad_f(xi_arr)
        jac_nz_arr[:] = Jprob.jac_g_nz(xi_arr)

        if dmode == 2
            hess_arr = unsafe_wrap(Array{CxxPtr{Float64}, 1}, hess.cpp_object, Jprob.n_hessblocks, own = true)

            s_last = Jprob.blockIdx[Jprob.n_hessblocks + 1] - Jprob.blockIdx[Jprob.n_hesslbocks]
            hess_last = unsafe_wrap(Array{Float64,1}, hess_arr[Jprob.n_hessblocks].cpp_object, Int32((s*(s+Int32(1)))//(Int32(2))), own = false)
            hess_last[:] = Jprob.last_hessBlock(xi_arr, lam_arr[Jprob.nVar + 1 : Jprob.nVar + Jprob.nCon])
        end

        if dmode == 3
            hessPTR_arr = unsafe_wrap(Array{CxxPtr{Float64}, 1}, hess.cpp_object, Jprob.n_hessblocks, own = true)
            hess_arr = Array{Array{Float64, 1}, 1}(undef, Jprob.n_hessblocks)
            for i = 1:Jprob.n_hessblocks
                Bsize = Jprob.blockIdx[i+1] - Jprob.blockIdx[i]
                hess_arr[i] = unsafe_wrap(Array{Float64,1}, hessPTR_arr[i].cpp_object, Int32((Bsize*(Bsize + Int32(1)))//Int32(2)), own = false)
            end

            hess_eval = Jprob.hess(xi_arr, lam_arr[Jprob.nVar + 1 : Jprob.nVar + Jprob.nCon])
            for i = 1:Jprob.n_hessblocks
                hess_arr[i][:] = hess_eval[i]
            end
        end
    end
    info[] = Int32(0);
    return
end

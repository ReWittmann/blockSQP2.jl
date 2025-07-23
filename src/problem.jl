mutable struct blockSQPProblem

    nVar::Cint
    nCon::Cint
    nnz::Cint
    blockIdx::Vector{Cint}
    
    vblocks::Vector{vblock}
    cond::Union{Condenser, Nothing}
    
    lb_var::Vector{Cdouble}
    ub_var::Vector{Cdouble}
    lb_con::Vector{Cdouble}
    ub_con::Vector{Cdouble}
    lb_obj::Cdouble
    ub_obj::Cdouble
    
    f::Function #Vector{Cdouble}[nVar] -> Cdouble
    g::Function #Vector{Cdouble}[nVar] -> Vector{Cdouble}[nCon]
    grad_f::Function #Vector{Cdouble}[nVar] -> Vector{Cdouble}[nVar]
    jac_g::Function #Vector{Cdouble}[nVar] -> Matrix{Cdouble, 2}[nCon, nVar]
    last_hessBlock::Function #Vector{Cdouble}[nVar] -> Vector{Cdouble}[lastBlocksize*(lastBlocksize + 1)/2] #lower diagonal elements
    hess::Function #Vector{Cdouble}[nVar] -> Vector{Vector{Cdouble}}[blocksize*(blocksize + 1)/2][length(blockIdx) - 1] #See utils.jl lower_to_full!, full_to_lower!
    
    continuity_restoration::Function #Vector{Cdouble}[nVar] -> Vector{Cdouble}[nVar]
    
    
    jac_g_nz::Function
    jac_g_row::Vector{Cint}
    jac_g_colind::Vector{Cint}

    x0::Vector{Cdouble}
    lambda0::Vector{Cdouble}

    blockSQPProblem(f::Function,
                    g::Function,
                    grad_f::Function,
                    jac_g::Function,
                    lb_var::Vector{Cdouble},
                    ub_var::Vector{Cdouble},
                    lb_con::Vector{Cdouble},
                    ub_con::Vector{Cdouble},
                    x0::Vector{Cdouble},
                    lambda0::Vector{Cdouble};
                    lb_obj::Cdouble = Cdouble(-Inf), 
                    ub_obj::Cdouble = Cdouble(Inf),
                    nnz::Cint = Cint(-1),
                    blockIdx::Vector{Cint} = [0, Cint(length(lb_var))],
                    vblocks::Vector{vblock} = vblock[],
                    cond::Union{Condenser, Nothing} = nothing,
                    jac_g_row::Vector{Cint} = Cint[],
                    jac_g_colind::Vector{Cint} = Cint[],
                    jac_g_nz::Function = fnothing, continuity_restoration::Function = fnothing,
                    last_hessBlock::Function = fnothing, hess::Function = fnothing
                    ) = new(Cint(length(lb_var)), Cint(length(lb_con)), nnz, blockIdx, vblocks, cond,
                            lb_var, ub_var, lb_con, ub_con, lb_obj, ub_obj,
                            f, g, grad_f, jac_g, jac_g_nz, continuity_restoration, last_hessBlock,
                            hess, jac_g_row, jac_g_colind, x0, lambda0
                            )
end
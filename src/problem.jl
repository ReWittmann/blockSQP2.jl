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
                    lb_var::Vector{FLOAT_T},
                    ub_var::Vector{FLOAT_T},
                    lb_con::Vector{FLOAT_T},
                    ub_con::Vector{FLOAT_T},
                    x0::Vector{FLOAT_T},
                    lambda0::Vector{FLOAT_T};
                    lb_obj::AbstractFloat = -Inf, 
                    ub_obj::AbstractFloat = Inf,
                    nnz::Integer = -1,
                    blockIdx::Vector{INT_T_1} = Int64[0, length(lb_var)],
                    vblocks::Vector{vblock} = vblock[],
                    cond::Union{Condenser, Nothing} = nothing,
                    jac_g_row::Vector{INT_T_2} = Int64[],
                    jac_g_colind::Vector{INT_T_3} = Int64[],
                    jac_g_nz::Function = fnothing, 
                    continuity_restoration::Function = fnothing,
                    last_hessBlock::Function = fnothing, 
                    hess::Function = fnothing
                    ) where {FLOAT_T <: AbstractFloat, INT_T_1 <: Integer, INT_T_2 <: Integer, INT_T_3 <: Integer} = 
                        new(Cint(length(lb_var)), Cint(length(lb_con)), Cint(nnz), [Cint(x) for x in blockIdx], 
                            vblocks, cond,
                            [Cdouble(x) for x in lb_var], [Cdouble(x) for x in ub_var], [Cdouble(x) for x in lb_con], [Cdouble(x) for x in ub_con], Cdouble(lb_obj), Cdouble(ub_obj),
                            f, g, grad_f, jac_g, last_hessBlock, hess,
                            continuity_restoration,
                            jac_g_nz, [Cint(x) for x in jac_g_row], [Cint(x) for x in jac_g_colind], 
                            x0, lambda0
                            )
end

function make_sparse!(B_prob::blockSQPProblem, nnz::Integer, jac_nz::Function, jac_row::Vector{T}, jac_col::Vector{T}) where T <: Integer
    B_prob.jac_g_nz = jac_nz
    B_prob.jac_g_row = [Cint(x) for x in jac_row]
    B_prob.jac_g_colind = [Cint(x) for x in jac_col]
    B_prob.nnz = Cint(nnz)
end
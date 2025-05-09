mutable struct blockSQPProblem

    nVar::Int32
    nCon::Int32
    nnz::Int32
    blockIdx::Array{Int32, 1}
    vblocks::Array{vblock, 1}
    
    lb_var::Array{Float64, 1}
    ub_var::Array{Float64, 1}
    lb_con::Array{Float64, 1}
    ub_con::Array{Float64, 1}
    lb_obj::Float64
    ub_obj::Float64

    f::Function
    g::Function
    grad_f::Function
    jac_g::Function
    continuity_restoration::Function
    last_hessBlock::Function
    hess::Function
    
    jac_g_nz::Function
    jac_g_row::Array{Int32, 1}
    jac_g_colind::Array{Int32, 1}

    x0::Array{Float64, 1}
    lambda0::Array{Float64, 1}

    blockSQPProblem(f::Function,
                    g::Function,
                    grad_f::Function,
                    jac_g::Function,
                    lb_var::Array{Float64, 1},
                    ub_var::Array{Float64, 1},
                    lb_con::Array{Float64, 1},
                    ub_con::Array{Float64, 1},
                    x0::Array{Float64, 1},
                    lambda0::Array{Float64,1};
                    lb_obj = -Inf, 
                    ub_obj = Inf,
                    nnz = Int32(-1),
                    blockIdx = [0, length(lb_var)],
                    vblocks::Array{vblock, 1} = vblock[],
                    jac_g_row = Int32[],
                    jac_g_colind = Int32[],
                    jac_g_nz = fnothing, continuity_restoration = fnothing,
                    last_hessBlock = fnothing, hess = fnothing
                    ) = new(Int32(length(lb_var)), Int32(length(lb_con)), nnz, blockIdx, vblocks,
                    lb_var, ub_var, lb_con, ub_con, lb_obj, ub_obj,
                    f, g, grad_f, jac_g, jac_g_nz, continuity_restoration, last_hessBlock,
                    hess, jac_g_row, jac_g_colind, x0, lambda0
    )
end
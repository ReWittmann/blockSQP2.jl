module SymbolicsHessUtilities

# Some helper functions to obtain the block structure as needed by blockSQP,
# placed into an extension to avoid the heavy Symbolics dependency. 
using Symbolics
using OptimizationBase, OptimizationBase.SciMLBase

using SparseArrays
using blockSQP

function blockSQP.compute_hessian_blocks(A::AbstractMatrix)
    _A = SparseArrays.SparseMatrixCSC(A)
    blockSQP.compute_hessian_blocks(_A)
end

function blockSQP.compute_hessian_blocks(A::SparseArrays.SparseMatrixCSC)
    n = size(A,1) # assume quadratic matrix here, since it is a hessian
    blockIdx = [0]
    max_row = zeros(n)
    for i=1:n
        col_start, col_end = A.colptr[i], A.colptr[i+1]-1
        if col_end >= col_start
            corresponding_rows = A.rowval[col_start:col_end]
            max_row[i] = maximum(corresponding_rows)
        else
            # TODO: How to handle empty column? For now: block continues
            max_row[i] = i+1
        end
    end

    for i=1:n-1
        if max_row[i] <= i && max_row[i+1] > i
            append!(blockIdx, i)
        end
    end
    append!(blockIdx, n)
    return blockIdx
end

function blockSQP.compute_hessian_blocks(f, g, num_x::Integer,
                 num_cons::Integer; parameters=[])
    lag(x, mu) = begin
        fx = f(x, parameters)
        cons = zeros(eltype(x), num_cons)
        g(cons, x)
        return fx + sum(mu .* cons)
    end
    input = Vector{Float64}(undef, num_x);
    sparse_hess = Symbolics.hessian_sparsity(x -> lag(x, ones(num_cons)), input)
    @info "Hessian sparsity structure:"
    display(sparse_hess)
    blockSQP.compute_hessian_blocks(sparse_hess)
end

function blockSQP.compute_hessian_blocks(prob::SciMLBase.OptimizationProblem)
    num_x = length(prob.u0)
    num_cons = prob.ucons === nothing ? 0 : length(prob.ucons)
    function cons_ip(cons,x)
        let _p = prob.p
            prob.f.cons(cons, x, _p)
            return cons
        end
    end
    blockSQP.compute_hessian_blocks(prob.f.f, cons_ip, num_x, num_cons; parameters=prob.p)
end


end #SymbolicsHelperUtilities

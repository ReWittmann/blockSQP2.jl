function lower_to_full!(arr1::Array{Float64, 1}, arr2::Array{Float64, 1}, n::Int32)
    for i = 1:n
        for j = 0:(i-1)
            arr2[i + j*n] = arr1[i + j*n - Int64(j*(j+1)//2)]
        end
    end

    for i = 1:n
        for j = i:(n-1)
            arr2[i+j*n] = arr1[(j+1) + (i-1)*n - Int64(i*(i-1)//2)]
        end
    end
end

function full_to_lower!(arr1::Array{Float64, 1}, arr2::Array{Float64, 1}, n::Int32)
    for i = 1:n
        for j = 0:(i-1)
            arr2[i + j*n - Int64(j*(j+1)//2)] = arr1[i + j*n]
        end
    end
end

function make_sparse(B_prob::blockSQPProblem, nnz::Integer, jac_nz::Function, jac_row::Array{Int32, 1}, jac_col::Array{Int32, 1})
    B_prob.jac_g_nz = jac_nz
    B_prob.jac_g_row = jac_row
    B_prob.jac_g_colind = jac_col
    B_prob.nnz = Int32(nnz)
end


function reduceConstrVio(Prob::Ptr{Nothing}, xi::CxxPtr{Float64}, info::CxxPtr{Int32})
    Jprob = unsafe_pointer_to_objref(Prob)::blockSQPProblem
    if Jprob.continuity_restoration == fnothing
        info[] = Int32(1)
    else
        xi_arr = unsafe_wrap(Array{Float64, 1}, xi.cpp_object, Jprob.nVar, own = false)
        xi_arr[:] = Jprob.continuity_restoration(xi_arr)
    end
    return
end

# Some functions to deal with sparsity
function compute_hessian_blocks(A::AbstractMatrix)
    _A = SparseArrays.SparseMatrixCSC(A)
    compute_hessian_blocks(_A)
end

function compute_hessian_blocks(A::SparseArrays.SparseMatrixCSC)
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

function compute_hessian_blocks(f, g, num_x::Integer,
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
    compute_hessian_blocks(sparse_hess)
end
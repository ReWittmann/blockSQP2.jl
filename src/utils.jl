function lower_to_full!(arr1::Vector{FLOAT_T_1}, arr2::Vector{FLOAT_T_2}, n::Integer) where {FLOAT_T_1 <: AbstractFloat, FLOAT_T_2 <: AbstractFloat}
    for i = 1:n
        for j = 0:(i-1)
            arr1[i + j*n] = arr2[i + j*n - Int64((j*(j+1))//2)]
        end
    end

    for i = 1:n
        for j = i:(n-1)
            arr1[i+j*n] = arr2[(j+1) + (i-1)*n - Int64(((i-1)*i)//2)]
        end
    end
end

function full_to_lower!(arr1::Vector{FLOAT_T_1}, arr2::Vector{FLOAT_T_2}, n::Integer) where {FLOAT_T_1 <: AbstractFloat, FLOAT_T_2 <: AbstractFloat}
    for i = 1:n
        for j = 0:(i-1)
            arr1[i + j*n - Int64((j*(j+1))//2)] = arr2[i + j*n]
        end
    end
end

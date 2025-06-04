module ComponentArraysExtension
using ComponentArrays
using blockSQP
@info "Loading ComponentArrays.jl extension for blockSQP..."

blockSQP.__lowerbounds(x::ComponentArray) = x[:]
blockSQP.__upperbounds(x::ComponentArray) = x[:]
blockSQP.__initial_values(x::ComponentArray) = x[:]

end
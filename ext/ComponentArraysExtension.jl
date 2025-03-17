module ComponentArraysExtension
using ComponentArrays
using blockSQP
@info "Loading ComponentArrays Extension..."

blockSQP.__lowerbounds(x::ComponentArray) = x[:]
blockSQP.__upperbounds(x::ComponentArray) = x[:]
blockSQP.__initial_values(x::ComponentArray) = x[:]

end
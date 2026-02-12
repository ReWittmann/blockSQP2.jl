# module ComponentArraysExtension
using ComponentArrays
# using blockSQP
# @info "Loading ComponentArrays.jl extension for blockSQP..."

blockSQP2.__lowerbounds(x::ComponentArray) = x[:]
blockSQP2.__upperbounds(x::ComponentArray) = x[:]
blockSQP2.__initial_values(x::ComponentArray) = x[:]

# end

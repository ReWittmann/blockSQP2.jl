module NLPstructures
using Reexport
@reexport using ComponentArrays
using DocStringExtensions

include("BlockDescriptor.jl")
export BlockDescriptor, blocktypeof

export Btype

"""Export the given names with a prefix, 
e.g. > @prefixexport pref SomeType
     becomes
     > prefSomeType = SomeType
     > export prefSomeType
"""
macro prefixexport(_prefix, args...)
    ret_exprs = Expr[]
    for arg in args
        push!(ret_exprs, Expr(:(=), Symbol("$(_prefix)$(arg)"), Symbol("$(arg)")))
        push!(ret_exprs, Expr(:export, Symbol("$(_prefix)$(arg)")))
    end
    return esc(Expr(:block, ret_exprs...))
end

macro nlpexport(args...)
    return esc(:(@prefixexport(nlp, $(args...))))
end

@nlpexport Block Variables Constraints Hess Matching Matchings

#Note: The AbstractVector should be of type AbstractVector{TupleBD{S}}
TupleBD{S} = Tuple{BlockDescriptor, Union{AbstractVector, S}} where {S <: Integer}
export TupleBD


include("NLPstructure.jl")

export NLPstructure, tagmap


include("LayoutUtils.jl")
export axlength, axsubrange, to_NamedTuple, to_ComponentArray, to_UR, to_Axis, axsubindex, axsubkeys,
       simple_vBlocks, simple_cBlocks, hessBlocks, hessBlockSizes, hessBlockZeroBasedIndex, hessBlockOneBasedIndex,
       has_parent, parent_of, has_parent_type, get_BlockDescriptors, tagmap


include("MultipleShootingSystemSC.jl")

@nlpexport StateMatching ParameterMatching ControlMatching
@nlpexport StateMatchings ParameterMatchings ControlMatchings
@nlpexport MultipleShootingSystemSC


end #module
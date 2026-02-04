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
       has_parent, parent_of, has_parent_type, get_BlockDescriptors, tagmap, subBlocks


include("MultipleShootingDF.jl")

@nlpexport MultipleShootingDF MSfree MSdependent
# export msSystemSC, msFree, msDependent, msMatching, msMatchings
# @nlpexport StateMatching ParameterMatching ControlMatching
# @nlpexport StateMatchings ParameterMatchings ControlMatchings
# @nlpexport MultipleShootingSystemSC

"""
    Extract structure information from objects of a modeling framework, 
    returning the variable layout and constraint layout as two TupleBD[]'s.
     
    This function requires extensions to add methods for specific frameworks.
    As these are highly specific, see the specific extension for required inputs.
"""
function extract_preLayouts()
    error("No extension method available to extract layout information from the given objects")
end

"""
    Extract structure information from objects of a modeling framework,
    returning a complete NLPlayout. See `extract_preLayouts` for specifics.
    
    This function should only be used if the complete layout information
    is available from the passed arguments and does not need to be modified 
    afterwards.
"""
function extract_NLPstructure(args...; kwargs...)
    prevLayout, precLayout = extract_preLayouts(args...; kwargs...)
    return NLPstructure((get_BlockDescriptors(prevLayout)...,), 
                     to_Axis(prevLayout), 
                     (get_BlockDescriptors(precLayout)...,), 
                     to_Axis(precLayout))
end

# export extract_preLayouts, extract_nlpLayout


end #module
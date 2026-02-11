module NLPstructures
using Reexport
@reexport using ComponentArrays
using DocStringExtensions

include("BlockDescriptor.jl")
export BlockDescriptor, blocktypeof


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

#Note: The AbstractVector is expected to be an AbstractVector{TupleBD{S}} or AbstractUnitRange
TupleBD = Tuple{BlockDescriptor, Union{AbstractVector, Integer}}
# TupleBD = Tuple{BlockDescriptor, T} where T <: Union{AbstractVector, Integer, Union{AbstractVector, Integer}}
export TupleBD

include("NLPlayout.jl")
export NLPlayout, tagmap

include("LayoutUtils.jl")
export to_NamedTuple, 
       to_ComponentArray, 
       to_UR, 
       to_Axis, 
       axlength, 
       axsubkeys, 
       axsublength,
       axsubrange, 
       simple_vBlocks, 
       simple_cBlocks, 
       has_parent, 
       has_parent_type, 
       has_parent_subtype,
       hessBlocks, 
       hessBlockSizes, 
       hessBlockIndexZeroBased, 
       hessBlockIndexOneBased,
       blockDescriptors, 
       tagmap, 
       subBlocks


include("MultipleShootingDF.jl")
@nlpexport MultipleShootingDF MSfree MSdependent

"""
    Extract structure information from objects of a modeling framework, 
    returning the variable layout and constraint layout as two TupleBD[]'s.
     
    This function requires extensions to add methods for specific frameworks.
    As these are highly specific, see the specific extension for required inputs.
"""
function get_preLayouts()
    error("No extension method available to extract layout information from the given objects")
end

"""
    Extract structure information from objects of a modeling framework,
    returning a complete NLPlayout. See `extract_preLayouts` for specifics.
    
    This function should only be used if the complete layout information
    is available from the passed arguments and does not need to be modified 
    afterwards.
"""
function get_layout(args...; kwargs...)
    prevLayout, precLayout = get_preLayouts(args...; kwargs...)
    return NLPlayout((blockDescriptors(prevLayout)...,), 
                     to_Axis(prevLayout), 
                     (blockDescriptors(precLayout)...,), 
                     to_Axis(precLayout))
end

export get_preLayouts, get_layout


end #module
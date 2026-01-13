module NLPstructures
using Reexport
@reexport using ComponentArrays
using DocStringExtensions

include("BlockDescriptor.jl")
export BlockDescriptor, blocktypeof

nlpBlock, nlpVariables, nlpConstraints, nlpHess, nlpMatching, nlpCondensingTarget = Block, Variables, Constraints, Hess, Matching, CondensingTarget
export nlpBlock, nlpVariables, nlpConstraints, nlpHess, nlpMatching, nlpCondensingTarget

#Note: The AbstractVector should be of type AbstractVector{TupleBD{S}}
TupleBD{S} = Tuple{BlockDescriptor, Union{AbstractVector, S}} where {S <: Integer}
export TupleBD


"""
$(TYPEDEF)
A type containing structure information about a nonlinear program,
e.g. layout data required for condensing of quadratic subproblems

# Fields
$(FIELDS)
"""
struct NLPstructure{VB, VL <: ComponentArrays.AbstractAxis, CB, CL <: ComponentArrays.AbstractAxis}
    "Collection of BlockDescriptor{Variables} describing blocks of variables"
    vBlocks::VB
    "A mapping of the BlockDescriptors to sections of the variables"
    vLayout::VL
    
    "Collection of all BlockDescriptor{Constraints} describing blocks of constraints"
    cBlocks::CB
    "See vLayout"
    cLayout::CL
end

export NLPstructure


include("LayoutUtils.jl")
export axlength, axsubrange, to_NamedTuple, to_ComponentArray, to_UR, to_Axis, axsubindex, hessblockindex, simple_vBlocks, simple_cBlocks, hessBlocks, has_parent



# abstract type AbstractBlockDescriptor end

# abstract type Block end
# abstract type Variables <: Block end
# abstract type Constraints <: Block end

# abstract type Hess <: Variables end
# abstract type Matching <: Constraints end


# nlpBlock, nlpVariables, nlpConstraints, nlpHess, nlpMatching = Block, Variables, Constraints, Hess, Matching
# export nlpBlock, nlpVariables, nlpConstraints, nlpHess, nlpMatching

# mutable struct BlockDescriptor{T <: Block} <: AbstractBlockDescriptor
#     tag::Symbol
#     flags::Tuple
#     attr::NamedTuple
#     function BlockDescriptor{arg_T}(args...; tag = gensym(), kwargs...) where {arg_T <: Block}
#         new{arg_T}(tag, (args...,), (; kwargs...))
#     end
#     BlockDescriptor(args...; kwargs...) = BlockDescriptor{Block}(args...; kwargs...)
    
#     #TODO Enforce certain attributes for certain block types, e.g. input and output blocks for Matching constraint descriptors
#     #function BlockDescriptor{arg_T}(args...; tag = gensym(), *=..., ... *=..., kwargs...) where {arg_T <: ...}
# end

# function blocktypeof(::BlockDescriptor{T}) where T
#     return T
# end


# _tags(arg::AbstractVector{Tuple{B, T}}) where {B <: AbstractBlockDescriptor, T} = arg .|> first .|> Base.Fix2(getfield, :tag)

# function _to_NamedTuple(arg::T, start::Vector{S}) where {T<:Integer, S <: Integer}
#     return UnitRange{T}(start[1], (start[1] += arg) - 1)
# end
# function _to_NamedTuple(arg::UnitRange, start::Vector{S}) where S <: Integer
#     i = start[1]
#     start[1] += length(arg)
#     return arg .- first(arg) .+ i
# end
# function _to_NamedTuple(arg::Vector{Tuple{B, T}}, start::Vector{S}) where {B <: AbstractBlockDescriptor, T <: Union{Integer, AbstractVector, Union{Integer, AbstractVector}}, S <: Integer}
#     return _to_NamedTuple(arg, start)
# end

# function to_NamedTuple(arg::Vector{Tuple{B, T}}, start::Vector{S} = [1]) where {B <: AbstractBlockDescriptor, T <: Union{Integer, AbstractVector, Union{Integer, AbstractVector}}, S <: Integer}
#     vals = Union{UnitRange, NamedTuple}[]
#     for elem in arg
#         push!(vals, _to_NamedTuple(last(elem), start))
#     end
#     return NamedTuple{(_tags(arg)...,)}((vals...,))
# end




# function to_ComponentArray(arg::Vector{Tuple{B, T}}) where {B <: AbstractBlockDescriptor, T <: Union{Integer, AbstractVector, Union{Integer, AbstractVector}}}
#     return arg |> to_NamedTuple |> ComponentArray
# end


# _to_UR!(arg::S, start::Vector{T} = [1]) where {S <: Integer, T <: Integer} = 
#     UnitRange{S}(start[1], (start[1] += arg) - 1)
# _to_UR!(arg::UnitRange{S}, start::Vector{T} = [1]) where {S <: Integer, T <: Integer} = 
#     UnitRange{S}(start[1], (start[1] += length(arg)) - 1)
# function _to_UR!(arg::Vector{Tuple{B, T}}, start::Vector{S} = [1]) where {B <: AbstractBlockDescriptor, T <: Union{AbstractVector, Integer}, S <: Integer}
#     ret = Tuple{B, AbstractVector}[]
#     for elem in arg
#         push!(ret, (first(elem), _to_UR!(last(elem), start)))
#     end
#     return ret
# end

# to_UR(arg::Vector{Tuple{B, T}}, start::S = 1) where {B <: AbstractBlockDescriptor, T <: Union{AbstractVector, Integer}, S <: Integer} = 
#     _to_UR!(arg, [start])


# axlength(arg::ComponentArrays.Axis) = sum([length(arg[key]) for key in keys(arg)])

# function to_ViewAxis(arg::UnitRange, start::Vector{S}) where S <: Integer
#     i = start[1]
#     start[1] += length(arg)
#     return arg .- first(arg) .+ i 
# end
# function to_ViewAxis(arg::T, start::Vector{S}) where {T <: Integer, S <: Integer}
#     return UnitRange(start[1], (start[1] += arg) - 1)
# end
# function to_ViewAxis(arg::AbstractVector, start::Vector{S}) where S <: Integer
#     ax = to_Axis(arg)
#     return ViewAxis(UnitRange(start[1], (start[1] += axlength(ax)) - 1), ax)
# end

# function to_Axis(arg::Vector{Tuple{B, T}}) where {B <: AbstractBlockDescriptor, T <: Union{Union{Integer,AbstractVector}, Integer, AbstractVector}}
#     names = Symbol[]
#     vals = Union{UnitRange,ViewAxis}[]
#     start = [1]
#     for elem in arg
#         push!(vals, to_ViewAxis(last(elem), start))
#     end
#     return Axis(NamedTuple{(_tags(arg)...,)}((vals...,)))
# end


# function Base.getindex(collection::T, ind::B) where {B <: AbstractBlockDescriptor, T <: ComponentArrays.ComponentArray}
#     return :parent in keys(ind.attr) ? collection[ind.attr[:parent]][ind.tag] : collection[ind.tag]
# end
# function Base.getindex(collection::T, @nospecialize(ind::BlockDescriptor)) where {T <: ComponentArrays.ComponentArray}
#     return :parent in keys(ind.attr) ? collection[ind.attr[:parent]][ind.tag] : collection[ind.tag]
# end


# function Base.view(collection::T, ind::B) where {B <: AbstractBlockDescriptor, T <: ComponentArrays.ComponentArray}
#     return :parent in keys(ind.attr) ? view(view(collection, ind.attr[:parent]), ind.tag) : view(collection, ind.tag)
# end
# function Base.view(collection::T, @nospecialize(ind::BlockDescriptor)) where {T <: ComponentArrays.ComponentArray}
#     return :parent in keys(ind.attr) ? view(view(collection, ind.attr[:parent]), ind.tag) : view(collection, ind.tag)
# end

# #Note: The AbstractVector should be of type AbstractVector{TupleBD{S}}
# TupleBD{S} = Tuple{BlockDescriptor, Union{AbstractVector, S}} where {S <: Integer}



# export BlockDescriptor, NLPstructure, to_Axis, to_UR, to_ComponentArray, to_NamedTuple, TupleBD, axlength, blocktypeof




#=
List of possible flags:
    :hess - designates block as a Hessian block
    :param - a block of parameter values
    :control - a block of control values
    :dstate - a block of differential state values
    :astate - a block of algebraic state values
    :dinit - a block of differential initial states
    :ainit - a block of algebraic initial states
    :cond_target - target of condensing, will walk its dependencies and eliminate all dependent variables
    :cond_barrier - will prevent condensing from eliminating this block and dependent block that the pass would reach.
    ...
=#

#=
    :parent - If the block is a subblock of another block, designates that parent block. Enables subblock-indexing
    :sub - The block's subblocks
    :([*]group, ::Any) - allows tagging several blocks as belonging to the same group, e.g. paramter/control values of the same parameter/control
    :(dep, Vector{Tuple{BlockDescriptor, BlockDescriptor}}) - designates variable blocks that the block depends on through constraints
    :(eq, Vector{BlockDescriptor}) - designtates blocks that linked via equality constraints, e.g. other variables of the same parameter
    ...
=#

end #module
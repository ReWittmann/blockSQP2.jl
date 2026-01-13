abstract type AbstractBlockDescriptor end

abstract type Block end
abstract type Variables <: Block end
abstract type Constraints <: Block end

abstract type Hess <: Variables end
abstract type Matching <: Constraints end

abstract type CondensingTarget <: Variables end

# nlpBlock, nlpVariables, nlpConstraints, nlpHess, nlpMatching = Block, Variables, Constraints, Hess, Matching
# export nlpBlock, nlpVariables, nlpConstraints, nlpHess, nlpMatching

mutable struct BlockDescriptor{T <: Block} <: AbstractBlockDescriptor
    tag::Symbol
    flags::Tuple
    attr::NamedTuple
    function BlockDescriptor{arg_T}(args...; tag = gensym(), kwargs...) where {arg_T <: Block}
        new{arg_T}(tag, (args...,), (; kwargs...))
    end
    BlockDescriptor(args...; kwargs...) = BlockDescriptor{Block}(args...; kwargs...)
    
    #TODO Enforce certain attributes for certain block types, e.g. input and output blocks for Matching constraint descriptors
    #function BlockDescriptor{arg_T}(args...; tag = gensym(), *=..., ... *=..., kwargs...) where {arg_T <: ...}
end

function blocktypeof(::BlockDescriptor{T}) where T
    return T
end

_tags(arg::AbstractVector{Tuple{B, T}}) where {B <: AbstractBlockDescriptor, T} = arg .|> first .|> Base.Fix2(getfield, :tag)
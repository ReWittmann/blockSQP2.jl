abstract type AbstractBlockDescriptor end

abstract type Block end

abstract type Variables <: Block end
abstract type Constraints <: Block end

abstract type Hess <: Variables end

abstract type Matching <: Constraints end

"""
Wrapper type for inner constructor selection
"""
abstract type Btype{T <: Block} <: Block end

mutable struct BlockDescriptor{T <: Block} <: AbstractBlockDescriptor
    tag::Symbol
    flags::Tuple
    attr::NamedTuple    
    function BlockDescriptor{arg_bT}(args...; tag = gensym(), kwargs...) where arg_bT <: Btype{arg_T} where arg_T
        new{arg_T}(tag, (args...,), (; kwargs...))
    end 
end

BlockDescriptor(args...; kwargs...) = BlockDescriptor{Btype{Block}}(args...; kwargs...)

function BlockDescriptor{arg_T}(args...; kwargs...) where arg_T <: Block
    return BlockDescriptor{Btype{arg_T}}(args...; kwargs...)
end

function blocktypeof(::BlockDescriptor{T}) where T
    return T
end


_tags(arg::AbstractVector{Tuple{B, T}}) where {B <: AbstractBlockDescriptor, T} = arg .|> first .|> Base.Fix2(getfield, :tag)

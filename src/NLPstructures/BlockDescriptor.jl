abstract type AbstractBlockDescriptor end

"""generic BlockDescriptor type"""
abstract type Block end

"""BlockDescriptor type for describing a section of variables"""
abstract type Variables <: Block end

"""BlockDescriptor type for describing a section of constraints"""
abstract type Constraints <: Block end

"""BlockDescriptor type for describing a section of variables 
corresponding to a block of a block-diagonal Hessian"""
abstract type Hess <: Variables end

"""BlockDescriptor type for describing a section of constraints
that demands a set of variables be some function of other variables"""
abstract type Matching <: Constraints end

"""BlockDescriptor type for describing a section of constraints
that is made up of several Matching blocks"""
abstract type Matchings <: Constraints end

# """BlockDescriptor type for describing a section of variables
# subject to matching conditions
# """
# abstract type Dependent <: Variables end

"""
Wrapper type for inner BlockDescriptor constructor selection
"""
abstract type Btype{T <: Block} <: Block end

mutable struct BlockDescriptor{T <: Block} <: AbstractBlockDescriptor
    tag::Symbol
    parent::Union{BlockDescriptor, Nothing}
    flags::Tuple
    attr::NamedTuple
    function BlockDescriptor{arg_bT}(args...; tag = gensym(), parent = nothing, kwargs...) where arg_bT <: Btype{arg_T} where arg_T
        new{arg_T}(tag, parent, (args...,), (; kwargs...))
    end
end


"""
    Methods for specific block types should provide sanity checks for the
    flags and attributes of the BlockDescriptor.
"""
function BlockDescriptor{arg_T}(args...; kwargs...) where arg_T <: Block
    return BlockDescriptor{Btype{arg_T}}(args...; kwargs...)
end

BlockDescriptor(args...; kwargs...) = BlockDescriptor{Btype{Block}}(args...; kwargs...)


function blocktypeof(::BlockDescriptor{T}) where T
    return T
end


_tags(arg::AbstractVector{Tuple{B, T}}) where {B <: AbstractBlockDescriptor, T} = arg .|> first .|> Base.Fix2(getfield, :tag)


function Base.show(io::IO, ::MIME"text/plain", B::BlockDescriptor{T}) where T <: Block
    _taketag(B::TB) where TB = B
    _taketag(B::TB) where TB <: BlockDescriptor = B.tag
    attr = NamedTuple((key, _taketag(B.attr[key])) for key in keys(B.attr))
    parent = isnothing(B.parent) ? "_" : B.parent.tag
    print(io, "BlockDescriptor{", T, "}(tag = ", B.tag, ", parent = ", parent, ", flags = ", B.flags, ", attr = ", attr, ")")
end

Base.show(io::IO, B::BlockDescriptor{T}) where T <: Block = Base.show(io, MIME"text/plain"(), B)
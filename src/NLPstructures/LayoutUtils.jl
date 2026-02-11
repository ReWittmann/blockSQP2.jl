

# function _to_NamedTuple(arg::T, start::Vector{S}) where {T<:Integer, S <: Integer}
function _to_NamedTuple(arg::Integer, start::Ref{I}) where I <: Integer
    return UnitRange(start[], (start[] += I(arg)) - I(1))
end
function _to_NamedTuple(arg::UnitRange, start::Ref{I}) where I <: Integer
    i = start[]
    start[] += length(arg)
    return arg .- I(first(arg)) .+ i
end
# function _to_NamedTuple(arg::Vector{Tuple{B, T}}, start::Ref{Int64}) where {B <: AbstractBlockDescriptor, T <: Union{Integer, AbstractVector, Union{Integer, AbstractVector}}}
function _to_NamedTuple(arg::Vector{Tuple{BD, T}}, start::Ref{I}) where {BD <: AbstractBlockDescriptor, T <: Union{AbstractVector, Integer}, I <: Integer}
    return to_NamedTuple(arg, start)
end

function to_NamedTuple(arg::Vector{Tuple{BD, T}}, start::Ref{I} = Ref{Int64}(1)) where {BD <: AbstractBlockDescriptor, T <: Union{AbstractVector, Integer}, I <: Integer}
    vals = Union{NamedTuple, UnitRange}[]
    for elem in arg
        push!(vals, _to_NamedTuple(last(elem), start))
    end
    return NamedTuple{(_tags(arg)...,)}((vals...,))
end


function to_ComponentArray(arg::Vector{Tuple{BD, T}}) where {BD <: AbstractBlockDescriptor, T <: Union{AbstractVector, Integer}}
    return arg |> to_NamedTuple |> ComponentArray
end


_to_UR!(arg::Integer, start::Ref{I} = Ref{Int64}(1)) where {I <: Integer} = 
    UnitRange(start[], (start[] += I(arg)) - I(1))
_to_UR!(arg::UnitRange, start::Ref{I} = Ref{Int64}(1)) where {I <: Integer} = 
    UnitRange(start[], (start[] += I(length(arg))) - I(1))
function _to_UR!(arg::Vector{Tuple{BD, T}}, start::Ref{I} = Ref{Int64}(1)) where {BD <: AbstractBlockDescriptor, T <: Union{AbstractVector, Integer}, I <: Integer}
    ret = Tuple{B, AbstractVector}[]
    for elem in arg
        push!(ret, (first(elem), _to_UR!(last(elem), start)))
    end
    return ret
end

to_UR(arg::Vector{Tuple{BD, T}}, start::Integer = 1) where {BD <: AbstractBlockDescriptor, T <: Union{AbstractVector, Integer}} = 
    _to_UR!(arg, Ref(start))


axlength(arg::ComponentArrays.Axis) = sum([length(arg[key]) for key in keys(arg)])

function to_ViewAxis(arg::UnitRange, start::Ref{I}) where I <: Integer
    i = start[]
    start[] += length(arg)
    return arg .- I(first(arg)) .+ i 
end
function to_ViewAxis(arg::Integer, start::Ref{I}) where {I <: Integer}
    return UnitRange(start[], (start[] += I(arg)) - I(1))
end
function to_ViewAxis(arg::AbstractVector, start::Ref{I}) where I <: Integer
    ax = to_Axis(arg)
    return ViewAxis(UnitRange(start[], (start[] += I(axlength(ax))) - I(1)), ax)
end

function to_Axis(arg::Vector{Tuple{BD, T}}) where {BD <: AbstractBlockDescriptor, T <: Union{AbstractVector, Integer}}
    names = Symbol[]
    vals = Union{UnitRange,ViewAxis}[]
    start = Ref(1)
    for elem in arg
        push!(vals, to_ViewAxis(last(elem), start))
    end
    return Axis(NamedTuple{(_tags(arg)...,)}((vals...,)))
end


function blockDescriptors(arg::T) where {T <: Union{Integer, UnitRange}}
    return BlockDescriptor[]
end

function blockDescriptors(arg::Tuple{BD, T}) where {BD <: AbstractBlockDescriptor, T <: Union{AbstractVector, Integer}}
    return vcat([first(arg)], blockDescriptors(last(arg)))
end

function blockDescriptors(arg::AbstractVector{Tuple{BD, T}}) where {BD <: AbstractBlockDescriptor, T <: Union{AbstractVector, Integer}}
    return arg .|> blockDescriptors |> splat(vcat)
end

function Base.getindex(collection::ComponentArray, @nospecialize(ind::BlockDescriptor))
    return !isnothing(ind.parent) ? collection[ind.parent][ind.tag] : collection[ind.tag]
end

function Base.getindex(ax::AbstractAxis, @nospecialize(ind::BlockDescriptor))
    return !isnothing(ind.parent) ? ax[ind.parent].ax[ind.tag] : ax[ind.tag]
end

function Base.view(collection::ComponentArray, @nospecialize(ind::BlockDescriptor))
    return !isnothing(ind.parent) ? view(view(collection, ind.parent), ind.tag) : view(collection, ind.tag)
end

function axsubkeys(ax::AbstractAxis, @nospecialize(ind::BlockDescriptor))
    return ax[ind].ax |> keys |> collect
end

function axsubrange(ax::AbstractAxis, @nospecialize(ind::BlockDescriptor))
    return !isnothing(ind.parent) ? axsubrange(ax, ind.parent)[ax[ind].idx] : ax[ind].idx
end

function axsublength(ax::AbstractAxis, @nospecialize(blk::BlockDescriptor))
    return length(axsubrange(ax, blk))
end

# function simple_blocks(blocks::BLKS, layout::NLPlayout{VB,VL,CB,CL}) where {BLKS, VB, VL <: ComponentArrays.Axis, CB, CL}
#     blocks |> Base.Fix1(filter, x-> (typeof(layout.vLayout[x].ax) <: Shaped1DAxis)) |> collect |> (arr -> sort(arr, by = (x -> first(axsubrange(layout.vLayout, x)))))
# end

"""Extract \"simple\" variable blocks, i.e. blocks that have no subblocks"""
function simple_vBlocks(layout::NLPlayout{VB,VL,CB,CL}) where {VB, VL <: ComponentArrays.Axis, CB, CL}
    return layout.vBlocks |> Base.Fix1(filter, x-> (typeof(layout.vLayout[x].ax) <: Shaped1DAxis)) |> collect |> (arr -> sort(arr, by = (x -> first(axsubrange(layout.vLayout, x)))))
    # return simple_blocks(layout.vBlocks, layout)
end

"""Extract \"simple\" constraint blocks, i.e. blocks that have no subblocks"""
function simple_cBlocks(layout::NLPlayout{VB,VL,CB,CL}) where {VB, VL, CB, CL <: ComponentArrays.Axis}
    return layout.cBlocks |> Base.Fix1(filter, x-> (typeof(layout.cLayout[x].ax) <: Shaped1DAxis)) |> collect |> (arr -> sort(arr, by = (x -> first(axsubrange(layout.cLayout, x)))))
    # return simple_blocks(layout.vBlocks, layout)
end


function has_parent(@nospecialize(b::BlockDescriptor), @nospecialize(p::BlockDescriptor))
    return !isnothing(b.parent) && (b.parent == p || has_parent(b.parent, p))
end

"""
Check whether the BlockDescriptor or any of its parents are of a certain block type.
"""
function has_parent_type(@nospecialize(blk::BlockDescriptor), TP::Type{T}) where T <: Block
    return blocktypeof(blk) isa Type{T} || (!isnothing(blk.parent) && has_parent_type(blk.parent, TP))
end

function has_parent_subtype(@nospecialize(blk::BlockDescriptor), TP::Type{T}) where T <: Block
    return blocktypeof(blk) <: T || (!isnothing(blk.parent) && has_parent_subtype(blk.parent, TP))
end

function hessBlocks(layout::NLPlayout{VB,VL,CB,CL}) where {VB, VL <: ComponentArrays.Axis, CB, CL}
    return layout.vBlocks |> Base.Fix1(filter, x -> blocktypeof(x) <: Hess) |> collect |> (arr -> sort(arr, by = (x -> first(axsubrange(layout.vLayout, x)))))
end

function hessBlockSizes(layout::NLPlayout{VB,VL,CB,CL}) where {VB, VL <: ComponentArrays.Axis, CB, CL}
    hBlocks = hessBlocks(layout)
    length(hBlocks) == 0 && return [axlength(layout.vLayout)]
    return hBlocks .|> Base.Fix1(axsublength, layout.vLayout)
end

function hessBlockIndexZeroBased(layout::NLPlayout{VB,VL,CB,CL}) where {VB, VL <: ComponentArrays.Axis, CB, CL}
    return cumsum(Int64[0, hessBlockSizes(layout)...])
end

function hessBlockIndexOneBased(layout::NLPlayout{VB,VL,CB,CL}) where {VB, VL <: ComponentArrays.Axis, CB, CL}
    return cumsum(Int64[1, hessBlockSizes(layout)...])
end

function subBlocks(layout::NLPlayout{VB,VL,CB,CL}, @nospecialize(blk::BlockDescriptor)) where {VB, VL <: ComponentArrays.Axis, CB, CL <: ComponentArrays.Axis}
    let __MP = tagmap(layout)
        if blk in layout.vBlocks
            return axsubkeys(layout.vLayout, blk) .|> x->__MP[x]
        elseif blk in layout.cBlocks
            return axsubkeys(layout.cLayout, blk) .|> x->__MP[x]
        end
    end
    error("Block is not part of the given NLPlayout")
end
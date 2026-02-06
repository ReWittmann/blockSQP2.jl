

function _to_NamedTuple(arg::T, start::Vector{S}) where {T<:Integer, S <: Integer}
    return UnitRange{T}(start[1], (start[1] += arg) - 1)
end
function _to_NamedTuple(arg::UnitRange, start::Vector{S}) where S <: Integer
    i = start[1]
    start[1] += length(arg)
    return arg .- first(arg) .+ i
end
function _to_NamedTuple(arg::Vector{Tuple{B, T}}, start::Vector{S}) where {B <: AbstractBlockDescriptor, T <: Union{Integer, AbstractVector, Union{Integer, AbstractVector}}, S <: Integer}
    return to_NamedTuple(arg, start)
end

function to_NamedTuple(arg::Vector{Tuple{B, T}}, start::Vector{S} = [1]) where {B <: AbstractBlockDescriptor, T <: Union{Integer, AbstractVector, Union{Integer, AbstractVector}}, S <: Integer}
    vals = Union{UnitRange, NamedTuple}[]
    for elem in arg
        push!(vals, _to_NamedTuple(last(elem), start))
    end
    return NamedTuple{(_tags(arg)...,)}((vals...,))
end


function to_ComponentArray(arg::Vector{Tuple{B, T}}) where {B <: AbstractBlockDescriptor, T <: Union{Integer, AbstractVector, Union{Integer, AbstractVector}}}
    return arg |> to_NamedTuple |> ComponentArray
end


_to_UR!(arg::S, start::Vector{T} = [1]) where {S <: Integer, T <: Integer} = 
    UnitRange{S}(start[1], (start[1] += arg) - 1)
_to_UR!(arg::UnitRange{S}, start::Vector{T} = [1]) where {S <: Integer, T <: Integer} = 
    UnitRange{S}(start[1], (start[1] += length(arg)) - 1)
function _to_UR!(arg::Vector{Tuple{B, T}}, start::Vector{S} = [1]) where {B <: AbstractBlockDescriptor, T <: Union{AbstractVector, Integer}, S <: Integer}
    ret = Tuple{B, AbstractVector}[]
    for elem in arg
        push!(ret, (first(elem), _to_UR!(last(elem), start)))
    end
    return ret
end

to_UR(arg::Vector{Tuple{B, T}}, start::S = 1) where {B <: AbstractBlockDescriptor, T <: Union{AbstractVector, Integer}, S <: Integer} = 
    _to_UR!(arg, [start])


axlength(arg::ComponentArrays.Axis) = sum([length(arg[key]) for key in keys(arg)])

function to_ViewAxis(arg::UnitRange, start::Vector{S}) where S <: Integer
    i = start[1]
    start[1] += length(arg)
    return arg .- first(arg) .+ i 
end
function to_ViewAxis(arg::T, start::Vector{S}) where {T <: Integer, S <: Integer}
    return UnitRange(start[1], (start[1] += arg) - 1)
end
function to_ViewAxis(arg::AbstractVector, start::Vector{S}) where S <: Integer
    ax = to_Axis(arg)
    return ViewAxis(UnitRange(start[1], (start[1] += axlength(ax)) - 1), ax)
end

function to_Axis(arg::Vector{Tuple{B, T}}) where {B <: AbstractBlockDescriptor, T <: Union{Union{Integer,AbstractVector}, Integer, AbstractVector}}
    names = Symbol[]
    vals = Union{UnitRange,ViewAxis}[]
    start = [1]
    for elem in arg
        push!(vals, to_ViewAxis(last(elem), start))
    end
    return Axis(NamedTuple{(_tags(arg)...,)}((vals...,)))
end


function blockDescriptors(arg::T) where {T <: Union{Integer, UnitRange}}
    return BlockDescriptor[]
end

function blockDescriptors(arg::Tuple{B, T}) where {B <: AbstractBlockDescriptor, T <: Union{Union{Integer,AbstractVector}, Integer, AbstractVector}}
    return vcat([first(arg)], blockDescriptors(last(arg)))
end

function blockDescriptors(arg::Vector{Tuple{B, T}}) where {B <: AbstractBlockDescriptor, T <: Union{Union{Integer,AbstractVector}, Integer, AbstractVector}}
    return arg .|> blockDescriptors |> splat(vcat)
end

    
# function Base.getindex(collection::T, ind::B) where {B <: AbstractBlockDescriptor, T <: ComponentArrays.ComponentArray}
#     return !isnothing(ind.parent) ? collection[ind.parent][ind.tag] : collection[ind.tag]
# end
function Base.getindex(collection::T, @nospecialize(ind::BlockDescriptor)) where {T <: ComponentArrays.ComponentArray}
    return !isnothing(ind.parent) ? collection[ind.parent][ind.tag] : collection[ind.tag]
end

# function Base.getindex(ax::AX, ind::B) where {AX <: ComponentArrays.AbstractAxis, B <: AbstractBlockDescriptor}
#     return !isnothing(ind.parent) ? ax[ind.parent].ax[ind.tag] : ax[ind.tag]
# end
function Base.getindex(ax::AX, @nospecialize(ind::BlockDescriptor)) where {AX <: ComponentArrays.AbstractAxis}
    return !isnothing(ind.parent) ? ax[ind.parent].ax[ind.tag] : ax[ind.tag]
end

# function axsubrange(ax::AX, ind::B)  where {AX <: ComponentArrays.AbstractAxis, B <: AbstractBlockDescriptor}
#     return !isnothing(ind.parent) ? axsubrange(ax, ind.parent)[ax[ind].idx] : ax[ind].idx
# end
function axsubrange(ax::AX, @nospecialize(ind::BlockDescriptor))  where {AX <: ComponentArrays.AbstractAxis}
    return !isnothing(ind.parent) ? axsubrange(ax, ind.parent)[ax[ind].idx] : ax[ind].idx
end

function axsublength(ax::AX, @nospecialize(blk::BlockDescriptor)) where {AX <: ComponentArrays.AbstractAxis}
    return length(axsubrange(ax, blk))
end

# function Base.view(collection::T, ind::B) where {B <: AbstractBlockDescriptor, T <: ComponentArrays.ComponentArray}
#     return !isnothing(ind.parent) ? view(view(collection, ind.parent), ind.tag) : view(collection, ind.tag)
# end
function Base.view(collection::T, @nospecialize(ind::BlockDescriptor)) where {T <: ComponentArrays.ComponentArray}
    return !isnothing(ind.parent) ? view(view(collection, ind.parent), ind.tag) : view(collection, ind.tag)
end


function axsubkeys(ax::AX, @nospecialize(ind::BlockDescriptor)) where {AX <: ComponentArrays.AbstractAxis}
    return ax[ind].ax |> keys |> collect
end


"""Extract \"simple\" variable blocks, i.e. blocks that have no subblocks"""
function simple_vBlocks(layout::NLPlayout{VB,VL,CB,CL}) where {VB, VL <: ComponentArrays.Axis, CB, CL}
    return layout.vBlocks |> Base.Fix1(filter, x-> (typeof(layout.vLayout[x].ax) <: Shaped1DAxis)) |> collect |> (arr -> sort(arr, by = (x -> first(axsubrange(layout.vLayout, x)))))
end

"""Extract \"simple\" constraint blocks, i.e. blocks that have no subblocks"""
function simple_cBlocks(layout::NLPlayout{VB,VL,CB,CL}) where {VB, VL, CB, CL <: ComponentArrays.Axis}
    return layout.cBlocks |> Base.Fix1(filter, x-> (typeof(layout.cLayout[x].ax) <: Shaped1DAxis)) |> collect |> (arr -> sort(arr, by = (x -> first(axsubrange(layout.cLayout, x)))))
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


parent_of(@nospecialize(b::BlockDescriptor), @nospecialize(p::BlockDescriptor)) = !isnothing(b.parent) && b.parent == p

function hessBlocks(layout::NLPlayout{VB,VL,CB,CL}) where {VB, VL <: ComponentArrays.Axis, CB, CL}
    return layout.vBlocks |> Base.Fix1(filter, x -> blocktypeof(x) <: Hess) |> collect |> (arr -> sort(arr, by = (x -> first(axsubrange(layout.vLayout, x)))))
end

function hessBlockSizes(layout::NLPlayout{VB,VL,CB,CL}) where {VB, VL <: ComponentArrays.Axis, CB, CL}
    hBlocks = hessBlocks(layout)
    length(hBlocks) == 0 && return [axlength(layout.vLayout)]
    # return hBlocks .|> Base.Fix1(getindex, layout.vLayout) .|> Base.Fix2(getfield, :idx) |> collect |> sort! .|> length
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
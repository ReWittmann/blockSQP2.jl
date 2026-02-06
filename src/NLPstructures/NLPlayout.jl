"""
$(TYPEDEF)
A generic type to store structure information of a nonlinear program,
e.g. Hessian block structure, dependency structure, etc.

# Fields
$(FIELDS)
"""
struct NLPlayout{VB, VL <: ComponentArrays.AbstractAxis, CB, CL <: ComponentArrays.AbstractAxis}
    "Collection of BlockDescriptor{Variables} describing blocks of variables"
    vBlocks::VB
    "A mapping of the BlockDescriptors to sections of the variables"
    vLayout::VL
    
    "Collection of all BlockDescriptor{Constraints} describing blocks of constraints"
    cBlocks::CB
    "See vLayout"
    cLayout::CL
    
    function NLPlayout(vB::VB, vL::VL, cB::CB, cL::CL) where {VB, VL <: ComponentArrays.AbstractAxis, CB, CL <: ComponentArrays.AbstractAxis}
        layout = new{VB, VL, CB, CL}(vB, vL, cB, cL)
        for vBlock in layout.vBlocks
            assert_layout(vBlock, layout)
        end
        #If HessBlocks are present, must cover the full range exactly.
        hBlocks = hessBlocks(layout)
        hlength = sum(hBlocks .|> Base.Fix1(axsubrange, layout.vLayout) .|> length)
        @assert hlength in (0, axlength(layout.vLayout))
        
        return layout
    end
end

function NLPlayout(pre_vL::AbstractVector{TupleBD}, pre_cL::AbstractVector{TupleBD})
    return NLPlayout((get_BlockDescriptors(pre_vL)...,), to_Axis(pre_vL), (get_BlockDescriptors(pre_cL)...,), to_Axis(pre_cL))
end

"""
    Methods for specific block types should provide sanity checks for the 
    layout the BlockDescriptor is embedded in.
"""
function assert_layout(::BlockDescriptor{B}, ::NLPlayout) where B <: Block 
    return nothing
end


function assert_layout(blk::BlockDescriptor{B}, layout::NLPlayout) where B <: Matchings
    @assert all(blocktypeof(subblk) <: Matching for subblk in subBlocks(layout, blk)) "All direct subblocks of a Matchings block must be of blocktype Matching"
end


tagmap(layout::NLPlayout) = Dict((BD.tag => BD) for BD in union(layout.vBlocks, layout.cBlocks))

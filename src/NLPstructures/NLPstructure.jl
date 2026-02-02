"""
$(TYPEDEF)
A generic type to store structure information of a nonlinear program,
e.g. Hessian block structure, dependency structure, etc.

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
    
    function NLPstructure(vB::VB, vL::VL, cB::CB, cL::CL) where {VB, VL <: ComponentArrays.AbstractAxis, CB, CL <: ComponentArrays.AbstractAxis}
        struc = new{VB, VL, CB, CL}(vB, vL, cB, cL)
        for vBlock in struc.vBlocks
            assert_layout(vBlock, struc)
        end
        return struc
    end
end

function NLPstucture(pre_vL::AbstractVector{TupleBD}, pre_cL::AbstractVector{TupleBD})
    return NLPstructure((get_BlockDescriptors(pre_vL)...,), to_Axis(pre_vL), (get_BlockDescriptors(pre_cL)...,), to_Axis(pre_cL))
end

"""
    Methods for specific block types should provide sanity checks for the 
    layout the BlockDescriptor is embedded in.
"""
function assert_layout(::BlockDescriptor{B}, ::NLPstructure) where B <: Block 
    return nothing
end


tagmap(S::NLPstructure) = Dict((BD.tag => BD) for BD in union(S.vBlocks, S.cBlocks))

"""
    System of NLP variables resulting from a multiple shooting,
    with the variables being ordered as [States,Controls] in each
    shooting interval. All direct children must be of blocktype Hess,
    which should have two children, the first being the states.
"""
abstract type MultipleShootingSystemSC <: Variables end


abstract type StateMatching <: Matching end
abstract type ParameterMatching <: Matching end
abstract type ControlMatching <: Matching end

"""
    System of NLP constraints that denotes a collection of
    state matching conditions for a MultipleShootingSystemSC.
    It must have one child for each shooting stage that is
    a matching for the dependent variables of that stage.
"""
abstract type StateMatchings <: Matchings end
abstract type ParameterMatchings <: Matchings end
abstract type ControlMatchings <: Matchings end


function BlockDescriptor{arg_T}(args...; kwargs...) where arg_T <: MultipleShootingSystemSC
    @assert :matchings in keys(kwargs) && blocktypeof(kwargs[:matchings]) == StateMatchings
    return BlockDescriptor{Btype{MultipleShootingSystemSC}}(args...; kwargs...)
end


function assert_layout(BD::BlockDescriptor{B}, struc::NLPstructure) where B <: MultipleShootingSystemSC
    MP = tagmap(struc)
    axsubBD(AX, ind) = let __MP = MP
        axsubkeys(AX, ind) .|> x->__MP[x]
    end
    
    #a) Direct subblocks must be Hessian blocks
    Hchildren = axsubBD(struc.vLayout, BD)
    @assert all(isa.(blocktypeof.(Hchildren), Type{Hess})) "All direct sub-BlockDescriptors of a MultipleShootingSystemSC must have BlockType \"Hess\""
    
    #b) One matching block for each shooting stage
    MTC = BD.attr[:matchings]
    Mchildren = axsubBD(struc.cLayout, MTC)
    @assert (length(Mchildren) == length(Hchildren) - 1) "Number matching conditions does not match Hessian block structure"
    
    #c) First Hessian block only has one control C subblock
    @assert (length(axsubkeys(struc.vLayout, first(Hchildren))) == 1) "First Hessian block of a MultipleShootingSystemSC must have exactly one subblock"
    
    #d) Hessian blocks except the first and the last must have one control/parameter C subblock and one state S subblock
    for i in eachindex(Hchildren)[2:end-1]
        subtags = axsubkeys(struc.vLayout, Hchildren[i])
        @assert length(subtags) == 2 "Every Hessian block of a MultipleShootingSystemSC except the first and last must have exactly two subblocks (free section and dependent section)"
        @assert (length(axsubrange(struc.vLayout, MP[first(subtags)])) == length(axsubrange(struc.cLayout, Mchildren[i-1]))) "Output dimension of matching $(i-1) does not match size of associated dependent variable block (block $i)"
    end
    
    return nothing
end


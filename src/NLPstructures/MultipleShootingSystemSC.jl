"""
    System of NLP variables resulting from a multiple shooting,
    with the variables being ordered as [States,Controls] in each
    shooting interval. All direct children must be of blocktype Hess,
    and have subblocks indicating dependent and free sections.
    
    msSystemSC
    |...Hess
        |...msFree
    |...Hess
        |...msDependent
        |...msFree
    |...Hess
    ...
    |...Hess
        |...msDependent
        |...msFree
    
    The dependent a free sections may have further sublayouts, marking
    e.g. controls and parameters.
    In addition, it must have an attribute :matchings referring to a 
    MultipleShootingMatchings constraint block indicating the state matching
    constraints of the system.
"""
abstract type msSystemSC <: Variables end

abstract type msFree <: Variables end
abstract type msDependent <: Variables end


"""
Matching condition for a the dependent states of a multiple shooting system.
"""
abstract type msMatching <: Matching end

"""
    System of NLP constraints that denotes a collection of
    state matching conditions for a MultipleShootingSystemSC.
    It must have one child for each shooting stage that is
    a matching for the dependent variables of that stage.
"""
abstract type msMatchings <: Matchings end



# abstract type ParameterMatching <: Matching end
# abstract type ParameterMatchings <: Matchings end
# abstract type ControlMatchings <: Matchings end


function BlockDescriptor{arg_T}(args...; kwargs...) where arg_T <: msSystemSC
    @assert :matchings in keys(kwargs) && blocktypeof(kwargs[:matchings]) == msMatchings
    return BlockDescriptor{Btype{msSystemSC}}(args...; kwargs...)
end


function assert_layout(BD::BlockDescriptor{B}, struc::NLPstructure) where B <: msSystemSC
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
    
    #c) First Hessian block only has one subblock of type msFree
    DFchildren = axsubBD(struc.vLayout, first(Hchildren))
    # @assert (length(axsubkeys(struc.vLayout, first(Hchildren))) == 1) "First Hessian block of a MultipleShootingSystemSC must have exactly one subblock"
    @assert (length(DFchildren) == 1 && blocktypeof(first(DFchildren)) == msFree) "First Hessian block of a MultipleShootingSystemSC must have exactly one subblock"

    
    #d) Hessian blocks except the first and the last must have one control/parameter C subblock and one state S subblock
    for i in eachindex(Hchildren)[2:end-1]
        # subtags = axsubkeys(struc.vLayout, Hchildren[i])
        DFchildren = axsubBD(struc.vLayout, Hchildren[i])
        # @assert length(subtags) == 2 "Every Hessian block of a MultipleShootingSystemSC except the first and last must have exactly two subblocks (free section and dependent section)"
        @assert (length(DFchildren) == 2 && blocktypeof(first(DFchildren)) == msDependent && blocktypeof(last(DFchildren)) == msFree) "Every Hessian block of a MultipleShootingSystemSC except the first and last must have exactly two subblocks (msFree and msDependent)"

        # @assert (length(axsubrange(struc.vLayout, MP[first(subtags)])) == length(axsubrange(struc.cLayout, Mchildren[i-1]))) "Output dimension of matching $(i-1) does not match size of associated dependent variable block (block $i)"
        @assert (length(axsubrange(struc.vLayout, first(DFchildren))) == length(axsubrange(struc.cLayout, Mchildren[i-1]))) "Output dimension of matching $(i-1) does not match size of associated dependent variable block (block $i)"
    end
    
    #e) Last Hessian block may have one dependent or one dependent and one free subblock.
    i = lastindex(Hchildren)
    DFchildren = axsubBD(struc.vLayout, Hchildren[i])
    @assert (1 <= length(DFchildren) <= 2 && blocktypeof(first(DFchildren)) == msDependent && (length(DFchildren) == 1 || blocktypeof(last(DFchildren)) == msFree)) "Every Hessian block of a MultipleShootingSystemSC except the first and last must have exactly two subblocks (msFree and msDependent)"
    @assert (length(axsubrange(struc.vLayout, first(DFchildren))) == length(axsubrange(struc.cLayout, Mchildren[i-1]))) "Output dimension of matching $(i-1) does not match size of associated dependent variable block (block $i)"
    
    return nothing
end

function assert_layout(BD::BlockDescriptor{B}, struc::NLPstructure) where B <: msMatchings
    MP = tagmap(struc)
    axsubBD(AX, ind) = let __MP = MP
        axsubkeys(AX, ind) .|> x->__MP[x]
    end
    
    for subBD in axsubBD(struc.cLayout, BD)
        @assert (blocktypeof(subBD) == msMatching) "Direct subblocks of a msMatchings must be of blocktype msMatching"
    end
end

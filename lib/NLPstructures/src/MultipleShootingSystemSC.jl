"""
    System of NLP variables resulting from a multiple shooting,
    with the variables being ordered as [States,Controls] in each
    shooting interval. All direct children must be of blocktype Hess,
    which should have two children, the first being the states.
"""
abstract type MultipleShootingSystemSC <: Variables end

"""
    System of NLP constraints that denotes a collection of
    matching constraints for a MultipleShootingSystemSC.
    It must have one child for each shooting stage that is
    a matching for the dependent variables of that stage.
"""
abstract type MultipleShootingMatchings <: Constraints end


function BlockDescriptor{arg_T}(args...; kwargs...) where arg_T <: MultipleShootingSystemSC
    @assert :matchings in keys(kwargs)
    return BlockDescriptor{Btype{MultipleShootingSystemSC}}(args...; kwargs...)
end


function assert_layout(BD::BlockDescriptor{B}, struc::NLPstructure) where B <: MultipleShootingSystemSC
    MP = tagmap(struc)
    axsubBD(AX, ind) = let __MP = MP
        axsubkeys(AX, ind) .|> x->__MP[x]
    end
    
    #a) Direct subblocks must be Hessian blocks
    Hchildren = axsubBD(struc.vLayout, BD)
    if !all(isa.(blocktypeof.(Hchildren), Type{Hess}))
        error("All direct sub-BlockDescriptors of a MultipleShootingSystemSC must have BlockType \"Hess\"")
    end
    
    #b) One matching block for each shooting stage
    MTC = BD.attr[:matchings]
    Mchildren = axsubBD(struc.cLayout, MTC)
    if !(length(Mchildren) == length(Hchildren) - 1)
        error("Number matching conditions does not match Hessian block structure")
    end
    
    #c) First Hessian block only has one control C subblock
    if !(length(axsubkeys(struc.vLayout, first(Hchildren))) == 1)
        error("First Hessian block of a MultipleShootingSystemSC must have exactly one subblock")
    end
    
    #d) Hessian blocks after the first have one control C subblock and one state S subblock
    for i in Iterators.drop(eachindex(Hchildren), 1)
        subtags = axsubkeys(struc.vLayout, Hchildren[i])
        if !(length(subtags) == 2)
            error("Every Hessian block of a MultipleShootingSystemSC after the first must have exactly two subblocks")
        elseif !(length(axsubrange(struc.vLayout, MP[first(subtags)])) == length(axsubrange(struc.cLayout, Mchildren[i-1])))
            error("Output dimension of matching $(i-1) does not match size of associated dependent variable block (block $i)")
        end
    end
    return nothing
end


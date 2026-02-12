using .NLPlayouts

"""
Create vblocks (variable block information) for blockSQP, not to be confused with NLPlayout.vBlocks, NLPlayouts.simple_vBlockes etc.
"""
function create_vblocks(struc::NLPlayout)
    return simple_vBlocks(struc) .|> x -> vblock(length(struc.vLayout[x]), has_parent_type(x, nlpMSdependent))
end

function create_condenser_args(struc::NLPlayout, add_dep_bounds = :all) #:none, :inactive, :all
    Dtargets = filter(x->(blocktypeof(x) <: nlpMultipleShootingDF), struc.vBlocks)
    if length(Dtargets) == 0
        return (nothing,)
    end
    
    Dmatchings = [DT.attr[:matchings] for DT in Dtargets]
    
    cBlocks = simple_cBlocks(struc)
    vBlocks = simple_vBlocks(struc)
    hBlocks = hessBlocks(struc)
    
    vblocks_args = [(size = length(struc.vLayout[x].idx), dependent = [false]) for x in vBlocks]
    
    cblocks = cBlocks .|> x -> cblock(length(struc.cLayout[x].idx))
    
    hsizes = hBlocks .|> x-> length(struc.vLayout[x].idx)
    
    targets = condensing_target[]
    for DT in Dtargets
        h0 = hBlocks[findfirst(Base.Fix2(has_parent, DT), hBlocks)]
        
        i0 = findfirst(Base.Fix2(has_parent,DT), vBlocks)
        i1 = findlast(Base.Fix2(has_parent,DT), vBlocks)
        for i in i0:i1
            vblocks_args[i].dependent[1] = has_parent_type(vBlocks[i], nlpMSdependent)
        end
        j0 = findfirst(Base.Fix2(has_parent, DT.attr[:matchings]), cBlocks)
        j1 = findlast(Base.Fix2(has_parent, DT.attr[:matchings]), cBlocks)
        N = count(Base.Fix2(has_parent, DT), hBlocks) - 1
        
        push!(targets, condensing_target(N, i0-1, i1, j0-1, j1)) 
    end
    
    vblocks = vblocks_args .|> x->vblock(x.size, x.dependent[1])
    return vblocks, cblocks, hsizes, targets, add_dep_bounds
end

function Condenser(struc::NLPlayout, add_dep_bounds = :all)
    @warn "Creating a Condenser from an NLPlayout is highly experimental for now"
    return Condenser(create_condenser_args(struc, add_dep_bounds)...)
end

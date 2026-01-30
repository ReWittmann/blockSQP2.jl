function create_vBlocks(struc::NLPstructure)
    Dtargets = filter(x->(blocktypeof(x) <: NLPstructures.MultipleShootingSystemSC), struc.vBlocks)
    return simple_vBlocks(struc) .|> x -> vblock(length(struc.vLayout[x]), :matching in keys(x.attr) && any(has_parent(x,T) for T in Dtargets))
end


function create_condenser_args(struc::NLPstructure, add_dep_bounds = :all) #:none, :inactive, :all
    Dtargets = filter(x->(blocktypeof(x) <: NLPstructures.MultipleShootingSystemSC), struc.vBlocks)
    if length(Dtargets) == 0
        return nothing
    end
    
    Dmatchings = [DT.attr[:matchings] for DT in Dtargets]
    
    Dcblocks = simple_cBlocks(struc)
    Dvblocks = simple_vBlocks(struc)
    Dhblocks = hessBlocks(struc)
    
    vblocks_args = [(size = length(struc.vLayout[x].idx), dependent = [false]) for x in Dvblocks]
    
    cblocks = Dcblocks .|> x -> cblock(length(struc.cLayout[x].idx))
    hsizes = Dhblocks .|> x-> length(struc.vLayout[x].idx)
    
    targets = condensing_target[]
    for T in Dtargets
        h0 = Dhblocks[findfirst(Base.Fix2(has_parent, T), Dhblocks)]
        
        i0 = findfirst(Base.Fix2(has_parent,T), Dvblocks)
        i1 = findlast(Base.Fix2(has_parent,T), Dvblocks)
        for i in i0:i1
            vblocks_args[i].dependent[1] = (:matching in keys(Dvblocks[i].attr) && !has_parent(Dvblocks[i], h0))
        end
        j0 = findfirst(Base.Fix2(has_parent, T.attr[:matchings]), Dcblocks)
        j1 = findlast(Base.Fix2(has_parent, T.attr[:matchings]), Dcblocks)
        N = count(Base.Fix2(has_parent, T), Dhblocks) - 1
        
        push!(targets, condensing_target(N, i0-1, i1, j0-1, j1)) 
    end
    
    vblocks = vblocks_args .|> x->blockSQP.vblock(x.size, x.dependent[1])
    return vblocks, cblocks, hsizes, targets, add_dep_bounds
    # return Condenser(vblocks, cblocks, hsizes, targets, add_dep_bounds)
end

create_condenser(struc::NLPstructure, add_dep_bounds = :all) = Condenser(create_condenser_args(struc, add_dep_bounds)...)

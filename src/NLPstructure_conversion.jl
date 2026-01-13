function create_vBlocks(struc::NLPstructure)
    Dtargets = filter(struc.vBlocks, x->(blocktypeof(x) <: nlpCondensingTarget))
    return simple_vBlocks(struc) .|> x -> vblock(length(struc.vLayout[x]), :matching in keys(x.attr) && any(has_parent(x,T) for T in Dtargets))
end

function create_condenser(struc::NLPstructure, add_dep_bounds = :all) #:None, :inactive, :all
    Dtargets = filter(struc.vBlocks, x->(blocktypeof(x) <: nlpCondensingTarget))
    Dmatchings = [Dtargets.attr[:matching]]
    
    Dcblocks = simple_cBlocks(struc)
    Dvblocks = simple_vBlocks(struc)
    Dhblocks = hessblocks(struc)
    
    vblocks = Dvblocks .|> x -> vblock(length(struc.vLayout[x].idx), false)
    cblocks = Dcblocks .|> x -> cblock(length(struc.cLayout[x].idx))
    hsizes = Dhblocks .|> x-> length(struc.vLayout[x].idx)
    
    targets = condensing_target[]
    for T in Dtargets
        i0 = findfirst(Base.Fix2(has_parent,T), Dvblocks)
        i1 = findlast(Base.Fix2(has_parent,T), Dvblocks)
        for i in i0:i1
            vblocks[i].dependent = (:matching in keys(Dvblocks[i]))
        end
        j0 = findfirst(Base.Fix2(has_parent, T.attr[:matching]), Dcblocks)
        j1 = findlast(Base.Fix2(has_parent, T.attr[:matching]), Dcblocks)
        N = count(Base.Fix2(has_parent, T), Dhblocks)
        
        push!(targets, condensing_target(N, i0-1, i1, j0-1, j1)) 
    end
    
    return Condenser(vblocks, cblocks, hsizes, targets, add_dep_bounds)
end
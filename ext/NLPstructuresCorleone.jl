module NLPstructuresCorleone
using LuxCore
using Corleone
using blockSQP.NLPstructures

@info "Loading Corleone.jl extension for blockSQP.NLPstructures..." 

#Prefixes: s - states, p - paramters, c - controls
#Suffixes: B - BlockDescriptor, L - Layout (::TupleBD[]), SUB - Sublayout (::TupleBD)
function NLPstructures.extract_preLayouts(shooting::SingleShootingLayer,
                    ps=LuxCore.initialparameters(Random.default_rng(), shooting),
                    st=LuxCore.initialstates(Random.default_rng(), shooting);
                    name = Symbol(gensym(), :_SingleShooting))
    # shootingSystemB = BlockDescriptor{nlpVariables}(tag = name)
    hessB = BlockDescriptor{nlpHess}(tag = Symbol(name, "_hess"))#, parent = shootingSystemB)        
        statesB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_states"), parent = hessB)
        paramB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_param"), parent = hessB)
        controlB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_control"), parent = hessB)
    
    statesparamcontrolL = TupleBD[(statesB, length(ps[:u0])), (paramB, length(ps[:p])), (controlB, length(ps[:controls]))]
    return TupleBD[(hessB, statesparamcontrolL)], TupleBD[]
end


function NLPstructures.extract_preLayouts(shooting::MultipleShootingLayer,
                    ps=LuxCore.initialparameters(Random.default_rng(), shooting),
                    st=LuxCore.initialstates(Random.default_rng(), shooting);
                    name = Symbol(gensym(), :_MultipleShooting))
    sMatchingsL = TupleBD[]
    pMatchingsL = TupleBD[]
    cMatchingsL = TupleBD[]
    
    sMatchingsB = BlockDescriptor{nlpMatchings}(tag = Symbol(name, "_stateMatchings"))
    pMatchingsB = BlockDescriptor{nlpMatchings}(tag = Symbol(name, "_parameterMatchings"))
    #Control matching are not supported yet as they may be partial
    cMatchingsB = BlockDescriptor{nlpConstraints}(tag = Symbol(name, "_controlMatchings"))
    
    shootingSystemB = BlockDescriptor{nlpMultipleShootingDF}(tag = name, matchings = sMatchingsB)    
    hessL = TupleBD[]
    
    traj, st = shooting(nothing, ps, st)
    CRL_matching_layout = traj.shooting
    # Assume Matchings to be ordered according to states-parameters-controls as 
    # (hopefully) returned by Corleone.get_shooting_constraints
    
    #First interval: No dependent section, layout according to indentation
    hessB = BlockDescriptor{nlpHess}(tag = Symbol(name, "_hess_1"), parent = shootingSystemB)
        freeB = BlockDescriptor{nlpMSfree}(tag = Symbol(name, "_free_1"), parent = hessB)
            statesB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_u0_1"), parent = freeB)
            paramB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_p_1"), parent = freeB)
            controlB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_control_1"), parent = freeB)
    #
    p1 = first(ps)
    hessLsub = (hessB, TupleBD[(freeB, [(statesB, length(p1[:u0])), 
                                        (paramB, length(p1[:p])), 
                                        (controlB, length(p1[:controls]))] )] )
    push!(hessL, hessLsub)
    sMatching_input = [freeB]
    
    for (k, interval) in Iterators.drop(enumerate(keys(ps)), 1)
        # a) Matching conditions to previous shooting interval
        sMatchingB = BlockDescriptor{nlpMatching}(tag = Symbol(name, "_stateMatching_$(k)"), parent = sMatchingsB, input = sMatching_input)
        pMatchingB = BlockDescriptor{nlpMatching}(tag = Symbol(name, "_paramterMatching_$(k)"), parent = pMatchingsB, input = [paramB])
        # Control matchings are not supported yet
        cMatchingB = BlockDescriptor{nlpConstraints}(tag = Symbol(name, "_controlMatching_$(k)"), parent = cMatchingsB)
        
        CRL_mtc_k = CRL_matching_layout[interval]
        push!(sMatchingsL, (sMatchingB, length(CRL_mtc_k[:u0])))
        push!(pMatchingsL, (pMatchingB, length(CRL_mtc_k[:p])))
        push!(cMatchingsL, (cMatchingB, length(CRL_mtc_k[:controls])))
        
        #b) Stage-local variables
        hessB = BlockDescriptor{nlpHess}(tag = Symbol(name, "_hess_$(k)"), parent = shootingSystemB)        
            depB = BlockDescriptor{nlpMSdependent}(tag = Symbol(name, "_dep_$(k)"), parent = hessB, matching = sMatchingB)
                statesB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_states_$(k)"), parent = depB)
            freeB = BlockDescriptor{nlpMSfree}(tag = Symbol(name, "_free_$(k)"), parent = hessB)
                paramB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_param_$(k)"), parent = freeB, matching = pMatchingB)
                controlB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_control_$(k)"), parent = freeB)
        
        ps_k = ps[interval]
        statesL = TupleBD[(statesB, length(ps_k[:u0]))]
        paramcontrolL = TupleBD[(paramB, length(ps_k[:p])), (controlB, length(ps_k[:controls]))]
        hessLsub = (hessB, TupleBD[(depB, statesL), (freeB, paramcontrolL)])
        push!(hessL, hessLsub)
        sMatching_input = [depB, freeB]
    end
    shootingSystemL = TupleBD[(shootingSystemB, hessL)]
    matchingsL = TupleBD[(sMatchingsB, sMatchingsL), (pMatchingsB, pMatchingsL), (cMatchingsB, cMatchingsL)]
    
    return shootingSystemL, matchingsL
end



end #NLPstructuresCorleone
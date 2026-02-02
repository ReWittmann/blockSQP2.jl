module NLPstructuresCorleone
using Corleone
using LuxCore
using blockSQP.NLPstructures

@info "Loading Corleone.jl extension for blockSQP.NLPstructures..." 

#Prefixes: s - states, p - paramters, c - controls
#Suffixes: B - BlockDescriptor, L - Layout (::TupleBD[]), SUB - Sublayout (::TupleBD) 
function NLPstructures.extract_preLayouts(shooting::MultipleShootingLayer,
                    ps=LuxCore.initialparameters(Random.default_rng(), shooting),
                    st=LuxCore.initialstates(Random.default_rng(), shooting);
                    name = Symbol(gensym(), :_shooting))
    sMatchingsL = TupleBD[]
    pMatchingsL = TupleBD[]
    cMatchingsL = TupleBD[]
    
    sMatchingsB = BlockDescriptor{msMatchings}(tag = Symbol(name, "_state_matchings"))
    pMatchingsB = BlockDescriptor{nlpMatchings}(tag = Symbol(name, "_parameter_matchings"))
    cMatchingsB = BlockDescriptor{nlpMatchings}(tag = Symbol(name, "_control_matchings"))
    
    shootingSystemB = BlockDescriptor{msSystemSC}(tag = name, matchings = sMatchingsB)    
    hessL = TupleBD[]
    
    CRL_matching_layout = shooting(nothing, ps, st).shooting
    # Assume Matchings to be ordered according to states-parameters-controls as 
    # returned by Corleone.get_shooting_constraints
    
    #First interval: No dependent section, layout according to indentation
    hessB = BlockDescriptor{nlpHess}(tag = Symbol(name, "_hess_$(1)"), parent = shootingSystemB)
        freeB = BlockDescriptor{msFree}(tag = Symbol(name, "_free_$(i)"), parent = hessB)
            statesB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_u0", parent = freeB))
            paramB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_p_$(i)", parent = freeB))
            controlB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_control_$(i)", parent = freeB))
    #
    p1 = first(ps)
    hessLsub = (hessB, TupleBD[(freeB, [(statesB, length(p1[:u0])), 
                                        (paramB, length(p1[:p])), 
                                        (controlB, length(p1[:controls]))]
                              )] )
    push!(hessL, hessLsub)
                              
    for (k, interval) in Base.tail(enumerate(keys(ps)))
        # a) Matching conditions to previous shooting interval
        sMatchingB = BlockDescriptor{msMatching}(tag = Symbol(name, "_state_matching_$(k)"), parent = sMatchingsB, input = [freeB])
        pMatchingB = BlockDescriptor{nlpMatching}(tag = Symbol(name, "_paramter_matching_$(k)"), parent = pMatchingsB, input = [paramB])
        # Currently, control matchings are not fully supported as they may be partial
        cMatchingB = BlockDescriptor{nlpMatching}(tag = Symbol(name, "_control_matching_$(k)"), parent = cMatchingsB)
        
        CRL_mtc_k = CRL_matching_layout[interval]
        push!(sMatchingsL, (sMatchingB, length(CRL_mtc_k[:u0])))
        push!(pMatchingsL, (pMatchingB, length(CRL_mtc_k[:p])))
        push!(cMatchingsL, (cMatchingB, length(CRL_mtc_k[:controls])))
        
        #b) Stage-local variables
        hessB = BlockDescriptor{nlpHess}(tag = Symbol(name, "_hess_$(k)"), parent = shootingB)        
            depB = BlockDescriptor{nlpDependent}(tag = Symbol(name, "_dep_$(k)"), parent = hessB, matching = sMatchingB)
                statesB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_states_$(k)"), parent = depB)
            freeB = BlockDescriptor{nlpFree}(tag = Symbol(name, "_free_$(k)"), parent = hessB)
                paramB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_param_$(k)"), parent = freeB, matching = pMatchingB)
                controlB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_control_$(k)", parent = freeB))
        
        ps_k = ps[interval]
        paramcontrolL = TupleBD[(paramB, length(ps_k[:p])), (controlB, length(ps_k[:controls]))]
        statesL = TupleBD[(statesB, length(ps_k[:u0]))]
        hessLsub = (hessB, TupleBD[(freeB, paramcontrolL), (depB, statesL)])
        push!(hessL, hessLsub)
    end
    shootingSystemL = TupleBD[(shootingSystemB, hessL)]
    matchingsL = TupleBD[(sMatchingsB, sMatchingsL), (pMatchingsB, pMatchingsL), (cMatchingsB, cMatchingsL)]
    
    return shootingSystemL, matchingsL
end


# function NLPstructures.extract_NLPstructure(
#     shooting::MultipleShootingLayer,
#     ps=LuxCore.initialparameters(Random.default_rng(), shooting),
#     st=LuxCore.initialstates(Random.default_rng(), shooting);
#     name = Symbol(gensym(), :_shooting)
# )
#     prevLayout, precLayout = preLayout(shooting, ps, st, name)
#     return NLPlayout((get_BlockDescriptors(prevLayout)...,), 
#                      to_Axis(prevLayout), 
#                      (get_BlockDescriptors(precLayout)...,), 
#                      to_Axis(precLayout))
# end



end
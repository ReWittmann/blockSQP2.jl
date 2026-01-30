module NLPstructuresCorleone
using Corleone
using LuxCore
using NLPstructures

@info "Loading Corleone.jl extension for NLPstructures..." 


#Prefixes: s - states, p - paramters, c - controls
#Suffixes: B - BlockDescriptor, L - Layout (::TupleBD[]), SUB - Sublayout (::TupleBD) 
function preLayout(shooting::MultipleShootingLayer,
                   ps=LuxCore.initialparameters(Random.default_rng(), shooting),
                   st=LuxCore.initialstates(Random.default_rng(), shooting);
                   name = Symbol(gensym(), :_shooting))
    sMatchingsL = TupleBD[]    
    pMatchingsL = TupleBD[]
    cMatchingsL = TupleBD[]

    sMatchingsB = BlockDescriptor{nlpMultipleShootingMatchings}(tag = Symbol(name, "_state_matchings"))
    pMatchingsB = BlockDescriptor{nlpMatching}(tag = Symbol(name, "_parameter_matchings"))
    cMatchingsB = BlockDescriptor{nlpMatching}(tag = Symbol(name, "_control_matchings"))
    
    
    shootingSystemB = BlockDescriptor{nlpMultipleShootingSystemSC}(tag = name, matchings = matchingsB)
        
    hessL = TupleBD[]
    
    CRL_matching_layout = shooting(nothing, ps, st).shooting
    
    #First interval: No dependent section
    hessB = BlockDescriptor{nlpHess}(tag = Symbol(name, "_hess_$(1)"), parent = shootingB)
        freeB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_free_$(i)"), parent = hessB)
            statesB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_u0", parent = freeB))
            paramB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_p_$(i)", parent = freeB))
            controlB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_control_$(i)", parent = freeB))
    
    #Layout
    p1 = first(ps)
    hessLT = (hessB, TupleBD[(freeB, [(initialstatesB, length(p1[:u0])), 
                                      (paramB, length(p1[:p])), 
                                      (controlB, length(p1[:controls]))]
                              )] )
    push!(hessL, hessLT)
                              
    for k, interval in Base.tail(enumerate(keys(ps)))
        sMatchingB = BlockDescriptor{nlpMatching}(tag = Symbol(name, "_state_matching_$(k)"), parent = matchingsB, input = [freeB])
        pMatchingB = BlockDescriptor{nlpMatching}(tag = Symbol(name, "_paramter_matching_$(k)"), parent = )
        cMatchingB = BlockDescriptor{nlpMatching}(tag = Symbol(name, "_control_matching_$(k)"))
        CRLmk = CRL_matching_layout[interval]
        push!(sMatchingsL, (sMatchingB, length(CRLmk[:u0])))
        push!(pMatchingsL, (pMatchingB, length(CRLmk[:p])))
        push!(cMatchingsL, (cMatchingB, length(CRLmk[:controls])))
        
        # Variables
        hessB = BlockDescriptor{nlpHess}(tag = Symbol(name, "_hess_$(k)"), parent = shootingB)        
            depB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_dep_$(k)"), parent = hessB, matching = smatchB)
                statesB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_states_$(k)"), parent = depB)
            freeB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_free_$(k)"), parent = hessB)
                paramB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_param_$(k)", parent = freeB))
                controlB = BlockDescriptor{nlpVariables}(tag = Symbol(name, "_control_$(k)", parent = freeB))
        
        pk = ps[interval]
        pcL = TupleBD[(paramB, length(pk[:p])), (controlB, length(pk[:controls]))]
        stL = TupleBD[(statesB, length(pk[:u0]))]
        hessSUB = (hessB, TupleBD[(freeB, param_controlL), (depB, statesL)])
        push!(hessL, hessSUB)
    end
    shootingSystemL = TupleBD[(shootingSystemB, hessL)]
    matchingsL = TupleBD[(sMatchingsB, sMatchingsL), (pMatchingsB, pMatchingsL), (cMatchingsB, cMatchingsL)]
    
    return shootingSystemL, matchingsL
end


function Layout(
    shooting::MultipleShootingLayer,
    ps=LuxCore.initialparameters(Random.default_rng(), shooting),
    st=LuxCore.initialstates(Random.default_rng(), shooting);
    name = Symbol(gensym(), :_shooting)
)
    prevLayout, precLayout = preLayout(shooting, ps, st, name)
    return NLPlayout((get_BlockDescriptors(prevLayout)...,), 
                     to_Axis(prevLayout), 
                     (get_BlockDescriptors(precLayout)...,), 
                     to_Axis(precLayout))
end
module NLPstructuresCorleoneOED
using LuxCore
using Corleone
using CorleoneOED
using blockSQP.NLPstructures

@info "Loading CorleoneOED extension for blockSQP.NLPstructures..." 

function NLPstructures.get_preLayouts(OED::OEDLayer,
                    ps=LuxCore.initialparameters(Random.default_rng(), OED),
                    st=LuxCore.initialstates(Random.default_rng(), OED);
                    name = Symbol(gensym(), :_OED))
    return NLPstructures.get_preLayouts(OED.layer, ps, st; name = name)
end

function NLPstructures.get_preLayouts(MultiExp::MultiExperimentLayer,
                    ps=LuxCore.initialparameters(Random.default_rng(), MultiExp),
                    st=LuxCore.initialstates(Random.default_rng(), MultiExp);
                    name = Symbol(gensym(), :_MultiExperiment))
    subPreLayouts = ((NLPstructures.get_preLayouts(MultiExp.layers[i], ps[Symbol(:experiment_, i)], st[Symbol(:experiment_, i)]; name = Symbol(name, :_, i))
                        for i in Base.OneTo(MultiExp.n_exp))...,)
    pre_vLayout = subPreLayouts .|> first |> Base.splat(vcat)
    pre_cLayout = subPreLayouts .|> last |> Base.splat(vcat)
    return pre_vLayout, pre_cLayout
end

end #NLPstructuresCorleoneOED
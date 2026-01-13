using NLPstructures


#u0, x1, u1, x2, u2, x3

h0 = BlockDescriptor{nlpHess}(tag = :h0)
u0 = BlockDescriptor{nlpVariables}(:control, parent = h0, tag = :u0)
p0 = BlockDescriptor{nlpVariables}(:parameter, parent = h0, tag = :p0)

c1 = BlockDescriptor{nlpMatching}(input = [u0, p0], tag = :c1)
h1 = BlockDescriptor{nlpHess}(tag = :h1)
x1 = BlockDescriptor{nlpVariables}(:dstate, matching = [c1], parent = h1, tag = :x1)
u1 = BlockDescriptor{nlpVariables}(:control, parent = h1, tag = :u1)
p1 = BlockDescriptor{nlpVariables}(:parameter, parent = h1, tag = :p1)

c2 = BlockDescriptor{nlpMatching}(input = [u1, p1, x1], tag = :c2)
h2 = BlockDescriptor{nlpHess}(tag = :h2)
x2 = BlockDescriptor{nlpVariables}(:dstate, matching = [c2], parent = h2, tag = :x2)
u2 = BlockDescriptor{nlpVariables}(:control, parent = h2, tag = :u2)
p2 = BlockDescriptor{nlpVariables}(:parameter, parent = h2, tag = :p2)

c3 = BlockDescriptor{nlpMatching}(input = [u2, p2, x2], tag = :c3)
h3 = BlockDescriptor{nlpHess}(tag = :h3)
x3 = BlockDescriptor{nlpVariables}(:dstate, :condensing_target, matching = [c3], parent = h3, tag = :x3)


vec_vLayout = TupleBD[(h0, [(u0, 2), (p0, 1)]), (h1, [(x1, 2), (u1, 2), (p1, 1)]), (h2, [(x2, 2), (u2, 2), (p2, 1)]), (h3,[(x3, 2)])]
vec_cLayout = TupleBD[(c1, 2), (c2, 2), (c3, 2)]

vBlocks = [u0, p0, x1, u1, p1, x2, u2, p2, x3]
cBlocks = [c1,c2,c3]

NLPlayout = NLPstructure(vBlocks, to_Axis(vec_vLayout), cBlocks, to_Axis(vec_cLayout))

#Dict((x.tag => x) for x in vBlocks)

x = 1
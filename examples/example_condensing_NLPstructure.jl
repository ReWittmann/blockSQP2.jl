using blockSQP
using blockSQP.NLPstructures


matchings = BlockDescriptor{nlpMatchings}(tag = :matchings)
MSsys = BlockDescriptor{nlpMultipleShootingDF}(tag = :MSsys, matchings = matchings)

h0 = BlockDescriptor{nlpHess}(tag = :h0, parent = MSsys)
h1 = BlockDescriptor{nlpHess}(tag = :h1, parent = MSsys)
h2 = BlockDescriptor{nlpHess}(tag = :h2, parent = MSsys)
h3 = BlockDescriptor{nlpHess}(tag = :h3, parent = MSsys)

u0 = BlockDescriptor{nlpMSfree}(:control, parent = h0, tag = :u0)

m1 = BlockDescriptor{nlpMatching}(input = [u0], tag = :m1, parent = matchings)
x1 = BlockDescriptor{nlpMSdependent}(:dstate, parent = h1, matching = [m1], tag = :x1)
u1 = BlockDescriptor{nlpMSfree}(:control, parent = h1, tag = :u1)

m2 = BlockDescriptor{nlpMatching}(input = [u1, x1], tag = :m2, parent = matchings)
x2 = BlockDescriptor{nlpMSdependent}(:dstate, matching = [m2], parent = h2, tag = :x2)
u2 = BlockDescriptor{nlpMSfree}(:control, parent = h2, tag = :u2)

m3 = BlockDescriptor{nlpMatching}(input = [u2, x2], tag = :m3, parent = matchings)
x3 = BlockDescriptor{nlpMSdependent}(:dstate, matching = [m3], parent = h3, tag = :x3)
u3 = BlockDescriptor{nlpMSfree}(:control, parent = h3, tag = :u3)

constr = BlockDescriptor{nlpConstraints}(tag = :constr)


vLayout = TupleBD[(MSsys, [(h0, [(u0, 1)]), (h1, [(x1, 2), (u1, 1)]), (h2, [(x2, 2), (u2, 1)]), (h3,[(x3, 2), (u3, 1)])])]
cLayout = TupleBD[(matchings, [(m1, 2), (m2, 2), (m3, 2)]), (constr, 1)]

vBlocks = [MSsys, h0, h1, h2, h3, u0, x1, u1, x2, u2, x3, u3]
cBlocks = [matchings, m1, m2, m3, constr]

layout = NLPlayout((vBlocks...,), to_Axis(vLayout), (cBlocks...,), to_Axis(cLayout))



using LinearAlgebra

con_jac_nz = Float64[-1,-2,1,1,-2,1,1,-1,1,-1,-2,1,1,-2,1,1,-1,1,-1,-2,1,1,1,1,1,1]
con_jac_row = Int32[0,1,6,0,2,6,1,3,6,2,3,6,2,4,6,3,5,6,4,5,6,4,6,5,6,6]
con_jac_colind = Int32[0,3,6,9,12,15,18,21,23,25,26]
con_jac = blockSQP.Sparse_Matrix(7, 10, con_jac_nz, con_jac_row, con_jac_colind)


full_block = Float64[0.75;-0.25;-0.25;; -0.25; 0.75; 0.25;; -0.25; 0.25; 1.25]

hess = Array{Array{Float64, 2}, 1}(undef, 4)
hess[1] = Float64[1.0;;]
hess[2] = full_block
hess[3] = full_block
hess[4] = full_block

grad_obj = ones(10)

lb_var = Float64[-0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3]
ub_var = Float64[0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
lb_con = Float64[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, -1.9]
ub_con = Float64[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.9]

using blockSQP
condenser = blockSQP.Condenser(layout)

H = zeros(10,10)
H[1,1] = hess[1][1,1]
H[2:4,2:4] = hess[2]
H[5:7,5:7] = hess[3]
H[8:10,8:10] = hess[4]
function to_dense(M::Integer, N::Integer, NZ::Vector{Float64}, ROW::Vector{Int32}, COLIND::Vector{Int32}) 
    A = zeros(M, N)
    for j = 1:N
        for i = COLIND[j]+1:COLIND[j+1]
            A[ROW[i]+1,j] = NZ[i]
        end
    end
    return A
end
A_con = to_dense(7,10,con_jac_nz, con_jac_row, con_jac_colind)
A = vcat(Diagonal(ones(10)),A_con)


using QPALM

#Full QP
model = QPALM.Model()
QPALM.setup!(model, Q=H, q=grad_obj, A=A, bmin=vcat(lb_var,lb_con), bmax=vcat(ub_var,ub_con), 
             print_iter=5, eps_rel=1e-10, eps_abs=1e-10)
results = QPALM.solve!(model)

xi = results.x
lam = results.y


#Condensed QP
condensed_h, condensed_jacobian, condensed_hess,
    condensed_lb_var, condensed_ub_var, condensed_lb_con, condensed_ub_con = 
    blockSQP.full_condense!(condenser, grad_obj, con_jac, hess, lb_var, ub_var, lb_con, ub_con)

conDENSEd_hess = condensed_hess[1]
conDENSEd_jacobian = zeros(7,4)
for j = 1:4
    for i = condensed_jacobian.colind[j]:(condensed_jacobian.colind[j + 1] - 1)
        conDENSEd_jacobian[condensed_jacobian.row[i + 1] + 1, j] = condensed_jacobian.nz[i + 1]
    end
end

A_cond = vcat(Diagonal(ones(4)), conDENSEd_jacobian)
QPALM.setup!(model, Q=conDENSEd_hess, q=condensed_h, A=A_cond, 
             bmin=vcat(condensed_lb_var,condensed_lb_con), 
             bmax=vcat(condensed_ub_var,condensed_ub_con), 
             print_iter=5, eps_rel=1e-10, eps_abs=1e-10)
condensed_results = QPALM.solve!(model)

xi_cond = condensed_results.x
lam_cond = condensed_results.y
xi_rest, lam_rest = blockSQP.recover_var_mult(condenser, xi_cond, lam_cond)

print("\n||xi - xi_rest||_∞ = ", maximum(xi - xi_rest))
print("\n||lam - lam_rest||_∞ ", maximum(lam - lam_rest), "\n")


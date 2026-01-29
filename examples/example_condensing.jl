include("../src/blockSQP.jl")
using .blockSQP


NZ = Float64[-1,-2,1,1,-2,1,1,-1,1,-1,-2,1,1,-2,1,1,-1,1,-1,-2,1,1,1,1,1,1]
ROW = Int32[0,1,6,0,2,6,1,3,6,2,3,6,2,4,6,3,5,6,4,5,6,4,6,5,6,6]
COLIND = Int32[0,3,6,9,12,15,18,21,23,25,26]
con_jac = blockSQP.Sparse_Matrix(7, 10, NZ, ROW, COLIND)


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
A_con = to_dense(7,10,NZ,ROW,COLIND)


print("###########################################\nQP info:\nH = \n")
display(H)
print("\n h = \n")
display(grad_obj)
print("A = \n")
display(A_con)
print("\nlb_var, ub_var = \n")
display(lb_var)
display(ub_var)
print("\nlb_con, ub_con = \n")
display(lb_con)
display(ub_con)

ID = zeros(10,10)
for i = 1:10
    ID[i,i] = 1.0
end
A = vcat(ID,A_con)

#Create structure data
vblocks = Array{blockSQP.vblock, 1}(undef, 7)
vblocks[1] = blockSQP.vblock(Int32(1), false)

vblocks[2] = blockSQP.vblock(Int32(2), true)
vblocks[3] = blockSQP.vblock(Int32(1), false)

vblocks[4] = blockSQP.vblock(Int32(2), true)
vblocks[5] = blockSQP.vblock(Int32(1), false)

vblocks[6] = blockSQP.vblock(Int32(2), true)
vblocks[7] = blockSQP.vblock(Int32(1), false)

cblocks = Array{blockSQP.cblock, 1}(undef, 4)
cblocks[1] = blockSQP.cblock(Int32(2))
cblocks[2] = blockSQP.cblock(Int32(2))
cblocks[3] = blockSQP.cblock(Int32(2))
cblocks[4] = blockSQP.cblock(Int32(1))

hsizes = Int32[1, 3, 3, 3]

targets = Array{blockSQP.condensing_target, 1}(undef, 1)
#3 stages, index of first free vblock, index after last dependent vblock, index of first condition, index after last condition
targets[1] = blockSQP.condensing_target(Int32(3), Int32(0), Int32(7), Int32(0), Int32(3))

cond = blockSQP.Condenser(vblocks, cblocks, hsizes, targets, Int32(2))

print("Created condenser julia struct. Condensing info:\n")
blockSQP.print_info(cond)


#Call a QP solver
using QPALM
model = QPALM.Model()
QPALM.setup!(model, Q=H, q=grad_obj, A=A, bmin=vcat(lb_var,lb_con), bmax=vcat(ub_var,ub_con), 
             print_iter=5, eps_rel=1e-10, eps_abs=1e-10)
results = QPALM.solve!(model)

xi = results.x
lam = results.y

condensed_h, condensed_jacobian, condensed_hess,
    condensed_lb_var, condensed_ub_var, condensed_lb_con, condensed_ub_con = 
    blockSQP.full_condense!(cond, grad_obj, con_jac, hess, lb_var, ub_var, lb_con, ub_con)



conDENSEd_hess = condensed_hess[1]


conDENSEd_jacobian = zeros(7,4)
for j = 1:4
    for i = condensed_jacobian.colind[j]:(condensed_jacobian.colind[j + 1] - 1)
        conDENSEd_jacobian[condensed_jacobian.row[i + 1] + 1, j] = condensed_jacobian.nz[i + 1]
    end
end


print("\n###########################################\nCondensed QP info:\ncondensed_h=\n") 
display(condensed_h) 
print("\ncondensed_jacobian=\n")
display(conDENSEd_jacobian)
print("\ncondensed_hess=\n")
display(condensed_hess[1])
print("\ncondensed_lb_var=\n")
display(condensed_lb_var)
print("\ncondensed_lb_con=\n")
display(condensed_lb_con)
print("\ncondensed_ub_con=\n")
display(condensed_ub_con)
print("\n")


ID = zeros(4,4)
for i = 1:4
    ID[i,i] = 1.0
end
A_2 = vcat(ID, conDENSEd_jacobian)

QPALM.setup!(model, Q=conDENSEd_hess, q=condensed_h, A=A_2, 
             bmin=vcat(condensed_lb_var,condensed_lb_con), 
             bmax=vcat(condensed_ub_var,condensed_ub_con), 
             print_iter=5, eps_rel=1e-10, eps_abs=1e-10)
condensed_results = QPALM.solve!(model)

xi_cond = condensed_results.x
lam_cond = condensed_results.y


print("\nPrimal solution of uncondensed QP:\n")
display(xi)
print("\nDual solution of uncondensed QP:\n")
display(lam)


print("\nPrimal solution of condensed QP:\n")
display(xi_cond)
print("\nDual solution of condensed QP:\n")
display(lam_cond)


xi_rest, lam_rest = blockSQP.recover_var_mult(cond, xi_cond, lam_cond)

print("\nPrimal restored solution:\n")
display(xi_rest)
print("\nDual restored solution:\n")
display(lam_rest)

print("\n||xi - xi_rest||_∞ = ", maximum(xi - xi_rest))
print("\n||lam - lam_rest||_∞ ", maximum(lam - lam_rest), "\n")



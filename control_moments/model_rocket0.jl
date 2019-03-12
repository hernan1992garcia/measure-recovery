using MultivariatePolynomials
using JuMP
using PolyJuMP
using SumOfSquares
using Mosek
using DynamicPolynomials
using Cubature
using PyPlot
###############################################
##function L: Calculates L operator value over v
##input: v: a test function in mons_test
########space_vars: the variables of the trajectory x
########f: such that x' = f
##output: coefficients in basis "mons" of the function L(v)
function L(v,space_vars,f,mons)
    dv_dt = differentiate(v,t)
    grad_v = differentiate(v,space_vars)
    if length(space_vars)==1
        L_v = dv_dt[1] + dot(grad_v[1],f)
    else
        L_v = dv_dt[1] + dot(grad_v,f)
    end
    return coefficients(L_v,mons)
end

####times_g: Calculates the matrix whose columns are coefficients of g*m
###input: mons: the global monomials; mon: monomial such that mon;g: polinomial which multiplicates
###otput: matrix of multiplication of g

function times_g(mons,mon,g)
    g_prod = coefficients(g*mon,mons)
    return g_prod
end

########
function moment_relax(mons,mons_test_tr,l,L_matrix,test_init_vec,test_final_vec,domain_rest)

coef_l = coefficients(l,mons)

size_mom = length(mons)
size_test = length(mons_test_tr)
l_rest = length(domain_rest)

modelo = Model(solver = MosekSolver())
@variable(modelo, A[1:size_mom])
@variable(modelo, B[1:size_test*(l_rest+1),1:size_test])

mons_matrix = mons_test_tr*mons_test_tr'

for i in 1:size_test
    for j in 1:size_test
        @constraint(modelo,B[i,j] == A[findin(mons.==mons_matrix[i,j],1)[1]])
    end
end
@SDconstraint(modelo,B[1:size_test,1:size_test] >= 0)

for k in 1:length(domain_rest)
    g = domain_rest[k]
    for i in 1:size_test
        for j in 1:size_test
            coef = dot(times_g(mons,mons_matrix[i,j],g),A)
            @constraint(modelo, B[k*size_test + i,j] == coef)
        end
    end
    @SDconstraint(modelo,B[k*size_test + 1:(k+1)*size_test,1:size_test] >= 0)
end

moments_vec = A

@constraint(modelo,test_init_vec+L_matrix'*moments_vec.== test_final_vec)

@objective(modelo,Min,dot(coef_l,moments_vec))

status = solve(modelo)
val_A = getvalue(A)

getvalue(dot(coef_l,moments_vec))
val_err = getvalue(sum(producto))
endcoefficients(L_v,mons)
#################################################
l_x = 1
l_u = 1
l_t = 1
##!!the last variables corresponds with control variable
@polyvar z[1:l_x+l_u+l_t]
t = z[1:l_t]
x = z[l_t+1:l_x+l_t]
u = z[l_x+l_t+1:l_x+l_t+l_u]
##d= used test function degree
d = 1
##max_d: maximum degree in restrictions
max_d = 1
##t_0: initial time
t_0 = 0
##t_F: Final time
t_F = 1
##x_0: initial state
x_0 = [0]
##x_F: final state
x_F = [0]
##f: such that x' = f(x)
f = u[1]
##l: the functional such that its integral will be optimized
l = x[1]^4 + (u[1]^2-1)^2
##the degree of l
deg_l = 4

##domain_rest: equations defining semialgebraic of X,T,U
domain_rest = [1-x[1],x[1],1-t[1],t[1],1-u[1],u[1]]
##deg_rest: maximum degree of restrictions
deg_rest = 1
################################################################
################################################################
##big_d: The work degree
big_d = max(d,max_d) + deg_l+ deg_rest
##mom_d: degree of the moments used in MR
mom_d = max(d,max_d) + deg_l
##test_d
test_d = 1

mons = []
for i in 0:big_d
    append!(mons,monomials(z,i))
end

mons = convert(Array{typeof(mons[1])},mons)

###mons_test_lim:test functions for final and initial measures
mons_test_lim = []
for i in 0:test_d
    append!(mons_test_lim,monomials(z[1:l_t+l_x],i))
end

mons_test_lim = convert(Array{typeof(mons_test_lim[1])},mons_test_lim)

###mons_test_tr:test functions for trajectory measure
mons_test_tr = []
for m in mons_test_lim
append!(mons_test_tr,monomials(m*f))
mons_test_tr = unique(mons_test_tr)
end

mons_test_tr = convert(Array{typeof(mons_test_tr[1])},mons_test_tr)

test_final_vec = [convert(Float64,subs(v,z[1:l_t+l_x] =>[t_F[1],x_F[1]])) for v in mons_test_lim]
test_init_vec = [convert(Float64,subs(v,z[1:l_t+l_x] =>[t_0[1],x_0[1]])) for v in mons_test_lim]
################
L_matrix = zeros(length(mons),length(mons_test_lim))

for i in 1:length(mons_test_lim)
    L_matrix[:,i] = L(mons_test_lim[i],x,f,mons)
end
###########################

using SymPy

##monomios: produces all monomials of degree exactly de in n variables
##input: de: the degree of the monomials;n: the number of variables
##output: M: the list f monomials
function monomios(de,n)
  de=convert(Int64,de)
  n=convert(Int64,n)
  M=[]
  for S in subsets(collect(1:(n+de-1)),n-1)
    E=[]
    E=Array{Int64}(E)
    for k in 1:n
      if k==1
        push!(E,S[k]-1)
      elseif k==n
        push!(E,de+n-1-S[n-1])
      else
      push!(E,S[k]-S[k-1]-1)
      end
    end
    push!(M,E)
  end
  return M
end

##given the Sympy polynomial expression p1 calculates coefficient of monomial mon
##vars: the variables with respect derivative is calculated(equal length of mon)
function coef_pol_sym(p1,mon,vars)
    partial = p1
    for i in 1:length(mon)
        partial = diff(partial,vars[i],mon[i])
    end

    prod_factoriales = prod([factorial(j) for j in mon])

    for i in 1:length(vars)
        partial = partial(vars[i]=>0)
    end
    return convert(Float64,(1/prod_factoriales)*partial)
end

####################################
function coef_pol_sym_tio(p0,mon0,vars0)
    partial = p0
    for i in 1:length(mon0)
        partial = diff(partial,vars0[i],mon0[i])
    end

    prod_factoriales = prod([factorial(j) for j in mon])

    for i in 1:length(vars)
        partial = partial(vars[i]=>0)
    end
    return convert(Float64,(1/prod_factoriales)*partial)
end



#####grad: gradient of a function
###input: f: the function; v:vector of variables
###output: grad(f): gradient of f respect v
function grad(f,v)
    deriv = [diff(f,v[i]) for i in 1:length(v)]
    return deriv
end

###D: calculates multivariables derivative of vector function f
###input: f: multivariate function;v: variables
###output D(f)
function D(f,v)
    Der = Array{typeof(f[1])}(length(f),length(v))
    for i in 1:length(f)
        Der[i,:] = grad(f[i],v)
    end
    return Der
end
##sustitucion: subtitute expr with values of "values"
##input: expr: the symbolic expression; var: variables, values: the values, in order corresponding to var
##output: numerical value corresponding to substitute var by values in expr
function sustitucion(expr, var,values)
    long = length(var)
    expr1 = expr
    for i in 1:long
        expr1 = subs(expr1,var[i]=> values[i])
    end
    expr2 = convert(Float64,expr1)
    return expr2
end


#////////////////////////////////////////////////////////////////
#---inverted pendulum equations
##1-arm inverted pendulum
N = 1
g = -9.8

###Constants
##Numeric values of parameters

l = [1.0]
#m= [m0;m1,..,]
m = [1;1.0]

#W: external forces
W = zeros(N+1)
#//////

#/////
u = symbols("u",real = true)
Rayleigh = zeros(N+1,N+1)
for i in 1:(N+1)
    Rayleigh[i,i] = 0.00
end
###Variables
theta = Array{SymPy.Sym}(N)
Dtheta = Array{SymPy.Sym}(N)
D2theta = Array{SymPy.Sym}(N)

q = symbols("q",real = true)
dq = symbols("dq",real = true)
d2q = symbols("d2q",real = true)


for i in 1:N
    theta[i] = symbols("theta$i",real = true)

    Dtheta[i] = symbols("dtheta$i",real = true)
    D2theta[i] = symbols("d2theta$i",real = true)
end

qs = Array{SymPy.Sym}(2,N+1)
qs[:,1] = [q;0]
for i in 2:(N+1)
    qs[:,i] = [l[i-1]*sin(theta[i-1]); l[i-1]*cos(theta[i-1])] + qs[:,i-1]
end

###dqs/dt
##v: the variables
v = [theta;q]
##dv_dt: derivative ofArray{SymPy.Sym}(N+1,N+1) variables respect time
dv_dt = [Dtheta;dq]
d2v_dt2 = [D2theta;d2q]

Dqs = Array{SymPy.Sym}(2,N+1)
for i in 1:(N+1)
    Dqs[:,i] = D(qs[:,i],v)*dv_dt
end

####K: kinetic energy
K = sum([0.5*m[i]*norm(Dqs[:,i])^2 for i in 1:(N+1)])

####P : potential energy
P = g*sum([m[i]*l[i-1]*qs[2,i] for i in 2:(N+1)])

###L: the Lagrangian
Lagr = K-P

grad_DqsL = grad(Lagr,dv_dt)
Dt_grad_DqsL = D(grad_DqsL,[v;dv_dt])*[dv_dt;d2v_dt2]

grad_qsL = grad(Lagr,v)

ext = -Rayleigh*dv_dt
#ext[N+1] = ext[N+1]+u
control = Array{SymPy.Sym}(N+1)
control[1:N] = [0.0 for i in 1:N]
control[N+1] = u
### eqs equal to zero
ecuaciones = Dt_grad_DqsL -grad_qsL - ext - control
ecuaciones = [simplify(i) for i in ecuaciones]

M_d2t = Array{SymPy.Sym}(N+1,N+1)
for i in 1:(N+1)
    for j in 1:(N+1)
        M_d2t[i,j] = simplify(coeff(ecuaciones[i],d2v_dt2[j]))
    end
end

F = ecuaciones - M_d2t*d2v_dt2
F = [simplify(F[i]) for i in 1:(N+1)]

M_inv = M_d2t^(-1)
################
################Construction Taylor approximation for SINGLE inverted pendulum
#d_expansion: degree of the trncated expansion + 1(high degree is missed in the method)
d_expansion = 2


#f_dyn: function such that y'' = f_dyn(t,y,y',u)
f_dyn = -M_inv*F

f_dyn_taylor = [removeO(series(f,v[1],0,d_expansion)) for f in f_dyn]
f_dyn_taylor = [simplify(f) for f in f_dyn_taylor]

######################################################
#####Converting problem in MultivariatePolynomials setting
########################################################

@polyvar z[1:2*length(v) + 2]
x = z[1:length(v)]
dx = z[length(v)+1:2*length(v)]
w = z[2*length(v) + 1]
t = z[2*length(v) + 2]
###d_expansion -1, because one degree is lost in Taylor expansion. +2 becacause dtheta gives two degrees more
###and the expansion was done just with respect theta1
d_f = d_expansion -1 + 2

monos = []
for i in 0:d_f
    append!(monos,monomios(i,length(z)-1))
end

vars = []
append!(vars,v)
append!(vars,dv_dt)
append!(vars,[u])
##f_pol: polynomials which defines de dynamics
f_pol = []
for f in f_dyn_taylor
    p = sum([coef_pol_sym(f,mon,vars)*prod(z[1:length(z)-1].^mon) for mon in monos])
    append!(f_pol,[p])
end

###################################
##################################
lim_inf_x = [0,-10]
lim_sup_x = [2*pi,10]
X_I = [[pi/17,0]]
X_F = [[0,0]]
t_F = 10
T = [0,t_F]
l = w^2 + x[1]^2
f = f_pol
domain_rest = [lim_sup_x[2]-x[2],x[2]-lim_inf_x[2],lim_sup_x[1]-x[1],x[1]-lim_inf_x[1],t,t_F-t]
d_test = 1
d_l = 2
d_rest = 1
d_f = 3
#############################
deg_mons = max(d_test + d_l,d_rest +d_f)
deg_mons_half = convert(Int64,floor(deg_mons/2))

mons =  []
for i in 0:deg_mons
    append!(mons,monomials(z,i))
end
mons = convert(Array{typeof(mons[1])},mons)

mons_half =  []
for i in 0:deg_mons_half
    append!(mons_half,monomials(z,i))
end
mons_half = convert(Array{typeof(mons[1])},mons_half)
size_half = length(mons_half)

#############################
(val_A,opt_value) = momentos0(X_I,X_F,T,l,f,domain_rest,d_test,d_l,d_rest,d_f)
integ_vec = val_A
matrizCoef = eye(length(integ_vec))
(A,errr) = error_SOS_aprox(integ_vec,mons,mons_half,matrizCoef)

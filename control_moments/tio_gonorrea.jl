using SymPy

x = symbols("x",real = true)
a = symbols("a",real = true)
b = symbols("b",real = true)
rho = symbols("rho",real = true)
mu = symbols("mu",real = true)

poly0 = rho*(b-a)
mon0 = rho*b
vars0 = [rho,b]
coef_pol_sym(poly0,mon0,vars0)

R = ((x-a)*rho+1)*(rho*(b-a)^2 +2*(b-a)) - (rho*(x-a)^2 +2*(x-a))*((b-a)*rho+1)
S = 2*mu*((b-a)*rho + 1)
T = R/S

gradiente = grad(T,[a,b])

gradiente2 = [simplify(u) for u in gradiente]

falso_del_tio = (a^2+b^2)/b^2

der = grad(falso_del_tio,[a,b])

#############################################
########POR FIN HIZO ALGO EL TIO
##############################################

x = symbols("x",real = true)
a = symbols("a",real = true)
b = symbols("b",real = true)
rho = symbols("rho",real = true)
mu = symbols("mu",real = true)
lam = symbols("lam",real = true)
m = symbols("m",real = true)
b_a = symbols("b_a",real = true)

a = z[1]
b = z[2]
x = z[3]
lam = z[4]
rho = z[5]
mu = z[6]
m = z[7]

poly1 = -2*mu*rho^2*(b-m) + lam*((b-a)*rho^2*(rho*(b-a)+1)+rho^2*(b-a))
poly0 = -2*mu*rho^2*(b-m) + lam*((b_a)*rho^2*(rho*(b_a)+1)+rho^2*(b_a))

deg_big = [coeff(expand(poly0),b_a^i) for i in 1:2]

coef_0 = simplify(poly0 - sum([deg_big[i]*b_a^i for i in 1:2]))

coefs = [coef_0,deg_big[1],deg_big[2]]
##########################
############Cuadratica
##########################
oj = -coefs[1] + ()

mons = convert(Array{typeof(mons[1])},mons)

#####################
#####################
@polyvar u
@polyvar v
@polyvar w

poly0 = u*(v-w)
poly0 = 1+3*u
mons0 = [1,u]

coefficients(poly0,mons0)

#################################################
########Particle pursuing a curve model
##################################################
@polyvar z[1:4]
t = z[1]
x = z[2:3]
u = z[4]
###w: the desired constant velocity
###n: the degree in series of sine and cosine
w = 4
n = 3
###Phi: the curve to be pursued(approximated as a polynomial)
Phi = [t;t^2]
X_I = [[0;0]]
X_F = [[1; 1]]
T = [0,1]

l = sum((x[1:2]-Phi).^2)

base_leg = [trans_legendre(u,n,T[2]) for n in 0:n]
sin_coefs = [quadgk(s->h(s)*sin(s),0,T[2])[1] for h in base_leg]
cos_coefs = [quadgk(s->h(s)*cos(s),0,T[2])[1] for h in base_leg]

approx_sin = sin_coefs'*base_leg
approx_cos = cos_coefs'*base_leg
#approx_sin = x-(1/6)*(u^3)
#approx_cos = 1-(1/2)*(u^2)

f1 = [w*approx_cos;w*approx_cos]
##domain_rest:For SOS version procedure.Must include everytime 1
#domain_rest = [1,-(X_F[1]-x)[1],-(X_F[1]-x)[2],-(x-X_I[1])[1],-(x-X_I[1])[2],t,1-t]
###domain_rest:For classical semidefinite procedure
domain_rest = [X_F[1][1]-x[1],X_F[1][2]-x[2],x[1]-X_I[1][1],x[2]-X_I[1][2],t-T[1],T[2]-t]#,u,2*pi-u]
d_test = 4
d_l = 2
d_rest = 1
d_f = 3

(val_A,opt_value) = momentos0(X_I,X_F,T,l,f1,domain_rest,d_test,d_l,d_rest,d_f)
##################################################################
#####Dirac delta's construction
###################################################################
deg_delta = 4

mons1 =  []
for i in 0:(1 + deg_delta)
    append!(mons1,monomials(z,i))
end
mons1 = convert(Array{typeof(mons1[1])},mons1)

expon_mons = [exponents(mons1[i]) for i in 1:length(mons1)]

############################################################
###Mixed Strategy
############################################################
#(val_A,opt_value) = momentosMix(X_I,X_F,T,l,f1,domain_rest,d_test,d_l,d_rest,d_f)
real_moments = val_A[1:length(mons1)]
####################################################
#######Recovering
####################################################
tiempos = collect(0.0:0.1:1)
x_delta = zeros(length(tiempos))
y_delta = zeros(length(tiempos))

for i in 1:length(tiempos)
    tiem = tiempos[i]
    delta_approx = dirac_legendre_n(t,T[2],tiem,deg_delta)
    coord_int_x = x[1]*delta_approx
    coord_int_y = x[2]*delta_approx
    coefs2_x = coefficients(coord_int_x,mons1)
    coefs2_y = coefficients(coord_int_y,mons1)
    x_delta[i] = real_moments'*coefs2_x
    y_delta[i] = real_moments'*coefs2_y
end

plot(x_delta,y_delta,marker = "*")
######
y2 = [convert(Float64,subs(Phi[2],t=>tiem)) for tiem in tiempos]
plot(tiempos,y2)

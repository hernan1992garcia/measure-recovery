#################################################
########Falling ball model
##################################################
@polyvar z[1:3]
t = z[1]
x = z[2:3]
u = z[3]

m = 1
g = -10

X_I = [[0;1]]
X_F = [[g/2+1; g+1]]
T = [0,1]
l = 0*x[1]
f1 = [u;g]
##domain_rest:For SOS version procedure.Must include everytime 1
#domain_rest = [1,-(X_F[1]-x)[1],-(X_F[1]-x)[2],-(x-X_I[1])[1],-(x-X_I[1])[2],t,1-t]
###domain_rest:For classical semidefinite procedure
domain_rest = [-(X_F[1]-x)[1],-(X_F[1]-x)[2],-(x-X_I[1])[1],-(x-X_I[1])[2],t,1-t]
d_test = 3
d_l = 0
d_rest = 1
d_f = 1
############################################################
###Construction of Exact Moments of the occupation measure
############################################################

##deg_delta:the degree of the truncated delta approximation
deg_delta = 2

mons1 =  []
for i in 0:(1 + deg_delta)
    append!(mons1,monomials(z,i))
end
mons1 = convert(Array{typeof(mons1[1])},mons1)

expon_mons = [exponents(mons1[i]) for i in 1:length(mons1)]

dyn_vec(s) = [s;(g/2)*s^2+s;g*s+1]

real_moments0 = [quadgk(prod(dyn_vec(t).^i),0,1)[1] for i in expon_mons]
noise = Normal(0,0.01)
real_moments = real_moments0 + rand(noise,length(real_moments0))
#############mix moments
real_moments = val_A[1:20]
############
####################################################
#######Recovering
####################################################

tiempos = collect(0.0:0.1:T[2])
y_delta = zeros(length(tiempos))

for i in 1:length(tiempos)
    tiem = tiempos[i]
    delta_approx = dirac_legendre_n(t,T[2],tiem,deg_delta)
    coord_int = x[1]*delta_approx
    coefs2 = coefficients(coord_int,mons1)
    y_delta[i] = real_moments'*coefs2
end

plot(tiempos,y_delta,marker="*")
#####################################################
#####Comparison Graph
######################################################
yy = [dyn_vec(i)[2] for i in tiempos]
plot(tiempos,yy)

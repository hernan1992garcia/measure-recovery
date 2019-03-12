########################################################
##########Spring model
########################################################
@polyvar z[1:3]
t = z[1]
x = z[2:3]
#u = z[3]

m = 1
k = 2

X_I = [[0;1]]
X_F = [[0; exp(-2*pi)]]
T = [0,2*pi]
l = 0*x[1]
f1 = [x[2];-k*x[1]-2*x[2]]
##domain_rest_DSOS:Must include everytime 1
domain_rest = [x[1]+1,1-x[1],x[2]+1,1-x[2],t-T[1],T[2]-t]
#domain_rest_DSOS = [1,x[1]+1,1-x[1],x[2]+1,1-x[2],t-T[1],T[2]-t]
d_test = 7
d_l = 0
d_rest = 1
d_f = 1

##################################################################
#####Dirac delta's approach
###################################################################

deg_delta = 7

mons1 =  []
for i in 0:(1 + deg_delta)
    append!(mons1,monomials(z,i))
end
mons1 = convert(Array{typeof(mons1[1])},mons1)

expon_mons = [exponents(mons1[i]) for i in 1:length(mons1)]

dyn_vec(s) = [s,exp(-s)*sin(s),exp(-s)*(cos(s)-sin(s))]

real_moments0 = [quadgk(w->prod(dyn_vec(w).^i),0,2*pi)[1] for i in expon_mons]

noise = Normal(0,1)
real_moments = real_moments0 + rand(noise,length(real_moments0))

############################################################
###Mixed Strategy
############################################################
deg_mons = d_test + d_l+d_rest +d_f
deg_mons_half = convert(Int64,floor(deg_mons/2))
mons0 =  []
for i in 0:deg_mons
    append!(mons0,monomials(z,i))
end
mons0 = convert(Array{typeof(mons0[1])},mons0)

#ort_basis = ort_random_matrix(length(mons0))
ort_basis = eye(length(mons0))

mons = ort_basis*mons0

size_half = convert(Int64,factorial(deg_mons_half + length(z))/(factorial(deg_mons_half)*factorial(length(z))))
mons_half = mons[1:size_half]

(val_A,opt_value) = momentosMix(X_I,X_F,T,l,f1,domain_rest,d_test,d_l,d_rest,d_f)
real_moments = val_A[1:length(mons1)]


####################################################
#######Recovering
####################################################
tiempos = collect(0.0:0.1:2*pi)
y_delta = zeros(length(tiempos))

for i in 1:length(tiempos)
    tiem = tiempos[i]
    delta_approx = dirac_legendre_n(t,2*pi,tiem,deg_delta)
    coord_int = x[1]*delta_approx
    coefs2 = coefficients(coord_int,mons1)
    y_delta[i] = real_moments'*coefs2
end

plot(tiempos,y_delta,marker = "*",label = "d_test = 7")

#####################################################
#####Comparison Graph
######################################################
yy = [dyn_vec(i)[2] for i in tiempos]
plot(tiempos,yy,label = "original")
title("legendre delta d = 7")
legend()

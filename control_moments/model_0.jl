
###############################################
#########Input Parameters
###############################################
@polyvar z[1:3]
t = z[1]
x = z[2]
u = z[3]

X_I = [0]
X_F = [1]
T = [0,1]
l = (u^2 - 1)^2#x^4 + (u^2 - 1)^2
f = u
domain_rest = [1-x,x,t,1-t]#[1-x,x,1-t,t,1-u,u]
d_test = 3
d_l = 4
d_rest = 1
d_f = 1
################################################
###########Pendulum Model
################################################
@polyvar z[1:3]
t = z[1]
x = z[2]
u = z[3]

m = 1
long = 1

X_I = [0]
X_F = [pi/2]
T = [0,pi/(2*(9.8)^(1/2))]
l = m*(0.5*(long*u)^2 + 9.8*long*(1-(x-x^3/6)))
f = u
domain_rest = [pi/2-x,x,t,pi/(2*(9.8)^(1/2))-t]#[1-x,x,1-t,t,1-u,u]
d_test = 1
d_l = 2
d_rest = 1
d_f = 1
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
f = [x[2];-k*x[1]-2*x[2]]
##domain_rest_DSOS:Must include everytime 1
domain_rest = [x[1]+1,1-x[1],x[2]+1,1-x[2],t-T[1],T[2]-t]
#domain_rest_DSOS = [1,x[1]+1,1-x[1],x[2]+1,1-x[2],t-T[1],T[2]-t]
d_test = 4
d_l = 0
d_rest = 1
d_f = 1
#######Calculating the matrix A of the trajectory SOS form m^t A m
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

##########################################
(val_A,opt_value) = momentosMix(X_I,X_F,T,l,f1,domain_rest,d_test,d_l,d_rest,d_f)
#(val_A,opt_value) = momentos0(X_I,X_F,T,l,f,domain_rest,d_test,d_l,d_rest,d_f)

#integ_vec = val_A[1:56]
matrizCoef = eye(length(integ_vec))
(A,errr) = error_SOS_aprox(integ_vec,mons0,mons_half,matrizCoef)
####################################################################################
#######################################################################################
####################################################################################
#eig_A = eig(A)
    diag_A = [A[i,i] for i in 1:size_half]
    mode_comp[j-1] = maximum(diag_A)
#    median_comp[j-1] = median(diag_A)
    mean_size[j-1] = size_half
end

plot(mean_size,mode_comp,color = :green,label = "maximum")
plot(mean_size,median_comp,color = :black,label = "median")
xlabel("# rows")
legend(loc="upper right",fancybox="true")
#plot(collect(2:1:11),mean_trunc,color = :red)
###################################################
###################################################
evaluation = [convert(Float64,MultivariatePolynomials.subs(f[1],z=>p)) for p in optimas]
tiempos = [p[1] for p in optimas]

ti = tiempos[tiempos.<=0.15]
ti = 0.1

dot(evaluation[tiempos.<=ti],coefes[tiempos.<=ti])

##-------------------Gives graphs in Figure 2

##epsilon: right limit of interval
epsilon = 0.1
##taus: the set of trial taus
taus = 0.0001:0.0003:epsilon
###Parameters of the recovery problem
k = 10
N = length(Puntos)
d = 30
n = 1
##Measure miu, with support Y
theta = pi/(20*N)
m = 120
indP = sample(1:N,k,replace = false)
#Y0-support of nu
Y0 = Puntos[indP]
##Y: support of original measure mu
Y = []
for i in 1:length(Y0)
  push!(Y,rotarEp(theta,Y0[i]))
end

measureCoefs =(1/k)*ones(k)

Y_approx = []
ind_approx =[]
for y in Y
  near = nearPoint(y)
  push!(Y_approx,near)
  push!(ind_approx,find_index(Puntos, near))
end

##-----
##-----"optimal" value
##-----

matrizMons = evaluamonomios(mons, Puntos)
matrizCoef = (1/(m)^(1/2))*coefMonomios(mons , d, m)
matrizMedicion = matrizCoef*matrizMons

weig_approx = measureCoefs
c_approx = zeros(N)
c_approx[ind_approx] = weig_approx
tol = tau(Y, measureCoefs, Y_approx, weig_approx)
tol = (1+tol/10)*tol
big_coef_aprox = solucion3(Y,measureCoefs,tol,matrizCoef,matrizMedicion)[2]

mons = []
for i in 0:d
  L = monomios(i,n+1)
  mons = append!(mons,L)
end

sparsities = []
for i in 1:length(taus)
  big_coef = solucion3(Y,measureCoefs,taus[i],matrizCoef,matrizMedicion)[2]
  push!(sparsities,length(big_coef))
end

##-----
##-----graph
##-----
plot(taus,sparsities,".",color = (0.0,0.0,0.0))
plot(tol,length(big_coef_aprox),"D",label = "optimal tau-sparsity",color = (0.0,0.0,0.0))
legend()
plot(taus,length(big_coef_aprox)*ones(length(taus)),"--",color=(0.0,0.0,0.0))
xlabel("tau")
ylabel("sparsity")

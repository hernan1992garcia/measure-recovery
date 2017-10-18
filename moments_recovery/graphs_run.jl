###---Calculates graphs in Figure 3

N = length(Puntos)
## d: degree of KSS polynomials used
d = 5
##n: dimension os the sphere
n = 1
##size_mean: the number of samples to be averaged
size_mean = 5
## m: number of measurements
m=120
## tol: the nu in approximated recovery- Theorem(A)
tol = 0.0003
#------------------------------
#------------------------------
#------------------------------

##mons: the monomials up to degree d
mons = []
for i in 0:d
  L = monomios(i,n+1)
  mons = append!(mons,L)
end

k = 3
epsilons = 0.0:0.001:pi/N
###Constructing support measures

errores = []
stand_dev = []
for epsilon in epsilons
  indP = sample(1:N,k,replace = false)
  #Y0-support of nu
  Y0 = Puntos[indP]
  ###Y-support of original measure mu
  Y = []
  for i in 1:length(Y0)
    push!(Y,rotarEp(epsilon,Y0[i]))
  end
  Y_approx = []
  ind_approx =[]
  for y in Y
    near = nearPoint(y)
    push!(Y_approx,near)
    push!(ind_approx,find_index(Puntos, near))
  end
  weig_mu = (1/k)*ones(k)
  weig_approx = weig_mu
  c_approx = zeros(N)
  c_approx[ind_approx] = weig_approx
  matrizMons = evaluamonomios(mons, Puntos)
  error_samples = []
  for cont in 1:size_mean
    matrizCoef = (1/(m)^(1/2))*coefMonomios(mons , d, m)
    matrizMedicion = matrizCoef*matrizMons
    c_recov = solucion2(Y,weig_mu,tol,matrizCoef,matrizMedicion)
    push!(error_samples,norm(c_recov-c_approx))
  end
  media = mean(error_samples)
  push!(errores,media)
  push!(stand_dev,norm(error_samples-media)/(length(error_samples)-1)^(1/2))
end

plot(epsilons,errores,"-" ,label = "error",color = (0.0,0.0,0.0))
plot(epsilons,errores-stand_dev,"--",label = "error± std. dev.",color = (0.0,0.0,0.0))
plot(epsilons,errores+stand_dev,"--",color = (0.0,0.0,0.0))
legend()
xlabel("angle")
ylabel("error")

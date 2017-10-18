##---

N = length(Puntos)
d = 10
n = 7
#Puntos = equidist(N)

##mons: the monomials up to degree d
mons = []
for i in 0:d
  L = monomios(i,n+1)
  mons = append!(mons,L)
end

ks = 1:50
M = [30,90,150,200]
style =  [":","-.","-","--"]
for t in 1:length(M)
  m= M[t]
  errores = []
  stand_dev = []
  for k in ks
    epsilon = 0.0
#####      m = Int(round(2*k*log(N)))
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
      tol = tau(Y, weig_mu, Y_approx, weig_approx)
      tol = (1+tol/10)*tol

      matrizMons = evaluamonomios(mons, Puntos)
        error_samples = []
        for cont in 1:20
          matrizCoef = (1/(m)^(1/2))*coefMonomios(mons , d, m)
          matrizMedicion = matrizCoef*matrizMons
          c_recov = solucion2(Y,weig_mu,tol,matrizCoef,matrizMedicion)
          push!(error_samples,norm(c_recov-c_approx))
        end
      media = mean(error_samples)
      push!(errores,media)
  end

  plot(ks,errores,style[t], label = string("m=",string(m)),color = (0.0,0.0,0.0))
  legend()

end
xlabel("sparsity")
ylabel("error")

savefig("d="*string(d)*"k's")
close()

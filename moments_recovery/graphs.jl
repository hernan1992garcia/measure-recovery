####---Functions used to calculated graphs in Figure 3

using StatsBase
using JuMP
using SCS
using Mosek
using PyPlot

##rotaEp: given a point, it is rotated in the direction of nearest point on the grid
##input: eps: angle of rotation; pto: point to be rotated
##output: Rotation point of pto eps degrees. Rotation is done in direction from
##pto to the nearest point to pto in the grid Puntos
function rotarEp(eps,pto)
  Puntos2 = []
  for p in Puntos
    if !(p==pto)
      push!(Puntos2,p-pto)
    end
  end
  normas = [norm(p) for p in Puntos2]
  minimos = find_index(normas,minimum(normas))
  direction = Puntos2[minimos[1]]+pto
  ortAxis = direction-dot(direction,pto)*pto
  ortAxis = ortAxis/norm(ortAxis)
  rotacion = cos(eps)*pto +sin(eps)*ortAxis
  return rotacion
end


###nearPoint: Gives the nearest point to a chosen point 'pto' over a grid 'Puntos'
###input:pto
###output: Near: nearest point to pto on Puntos
function nearPoint(pto)
  n = length(pto)
  Puntos2 = []
  for p in Puntos
    push!(Puntos2,p-pto)
  end
  normas = [norm(p) for p in Puntos2]
  ordNorm = sort(normas)
  min0 = find_index(normas,ordNorm[1])
  Near = Puntos[min0]
  return Near
end


##nearGridSupp:
##input:Y: support of the measure to be approximated
##output: supp: the points on Puntos that are nearest to points on Y
function nearGridSupp(Y)
  supp = []
  for y in Y
    push!(supp,nearPoint(y))
  end
  return supp
end
##tau: Calulates tau value of Lemma 4.12
##input:Y: support of measure mu; Y0: support of measure nu
##weig_mu: weights of mu; weig_nu: weights of nu
##output: tau value of Lemma 4.12
function tau(Y, weig_mu, Y0, weig_approx)
  N = length(Y0)
  k = length(Y)
  ##V matrix
  V = zeros(N,N)
  for s in 1:N
    for t in 1:N
      V[s,t] = ((1+dot(Y0[s],Y0[t]))/2)^d
    end
  end
  ## A matrix
  A = zeros(N,k)
  for s in 1:N
    for t in 1:k
      A[s,t] = ((1+dot(Y0[s],Y[t]))/2)^d
    end
  end
  ##D matrix
  D = zeros(k,k)
  for s in 1:k
    for t in 1:k
      D[s,t] = ((1+dot(Y[s],Y[t]))/2)^d
    end
  end
  t=weig_approx'*V*weig_approx-2*weig_mu'A*weig_approx + weig_approx'*D*weig_approx
  return t[1]
end

##find_index: find position of a point in a list
##input: list: list of points; item: point for which position is wanted
##output: ind: index in list of item or -1 if item is not in list
function find_index(list,item)
  ind = -1
  for l in 1:length(list)
      if norm(list[l]-item)<=0.0000001
        ind = l
        return ind
      end
  end
  return ind
end

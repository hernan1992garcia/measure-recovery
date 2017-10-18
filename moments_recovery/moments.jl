##-------This is a set of functions which willbe used in the develop of our experiments
##-------and also can be used to apply theory of our work in practice.

using Distributions
using MultiPoly
using JuMP
using SCS
using Mosek
using Iterators
using DataStructures

##Sum:sum "BigFloat" scalars
##input: L: a list of "BigFloat" numbers
## output: suma:the sum of elements in L

function Sum(L)
    suma=BigFloat(0.0)
    for i in L
        suma=suma+i
    end
    return suma
end

##expFact: calculates log(n!)
##input: n: a nonnegative integer
##output: sum(logs): the value of log(n!)
function expFact(n)
    if n==0
        return n

    elseif n>0
    logs=zeros(n)
      for i in 1:n
         logs[i]=log(i)
      end
         return sum(logs)
    end
end


## Comb(d,alpha):
##input: d: integer scalar; alpha: integer vector
##output: the result of compute C(d,alpha), as indicated in Definition 4.1 in the paper
function Comb(d,alpha)

  den1=expFact(d-sum(alpha))
  den2=zeros(length(alpha))
    for i in 1:length(alpha)
        den2[i]=expFact(alpha[i])
    end
      den3=BigFloat(sum(den2))

       return exp(expFact(d)-(den1+den3))
end

##coef: cofficient of monomial KSS with degree d and exponent alpha
##input: d: degree of the polynomials; alpha: exponent vector
##output: a: the coefficient A_alpha
function coef(d,alpha)
  variance=Comb(d,alpha)/2^d
  stDev=variance^(1/2)
  dist=Normal(0.0,stDev)
  a=BigFloat(rand(dist,1)[1])
  return a
end

##apoint: a point in S^(n+1) with spherical coordinates phi
##input: phi: the vector of n-1 angles in the S_n spherical parametrization. last one is the 2pi angle;
##output, the point in the sphere with parametric coordinates phi
function apoint(phi)

  n=length(phi)

  A=sin.(phi)
  B=cos.(phi)

  x=zeros(n+1)

  for i in 2:n
   x[i]=prod(A[1:(i-1)])*B[i]
  end
 x[1]=B[1]
 x[n+1]=prod(A)

 return x
end

##puntos: N points in S^n with spherical coordinates generated uniformly i.i.d
##inputs:n+1: number of entries in the points; N: number of points
##output: Puntos: the list of points in cartesian coordinates
function puntos(n,N)
  Angulos=zeros(N,n)
  DistAng0=Uniform(0,2*pi)
  DistAng1=Uniform(0,pi)
  for i in 1:N
      ang=zeros(n)
      ang0=rand(DistAng0,1)
      ang1=rand(DistAng1,n-1)
      ang[1:(n-1)]=ang1
      ang[n]=ang0[1]

      Angulos[i,:]=ang
   end
  Puntos=[]
  for j in 1:N
   Puntos=push!(Puntos,apoint(Angulos[j,:]))
  end
  return Puntos
end

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


##evaluamonomios: matrix with monomials evaluated at the points
##inputs: expMonom:list of monomial exponents; Puntos: points to be evaluated
##output: M: a matrix wich contains evaluation a monomial in expMon ovr all points in Puntos
##in each row
function evaluamonomios(expMonom,Puntos)
  M =zeros(length(expMonom),length(Puntos))
  for m in 1:length(expMonom)
      for p in 1:length(Puntos)
       M[m,p]=prod(Puntos[p].^expMonom[m])
      end
  end
  return M
end


##coefMonomios: Constructs coefficients of KSS polynomials given the list of exponents
##input: ExpMonom: list of monomial exponents; d: degree of monomials; m: number of polynomials
##output: N: matrix with coefficients of each KSS polynomial in each row. Monomials label columns
function coefMonomios(ExpMonom,d,m)
  N=zeros(m,length(ExpMonom))
  for i in 1:m
    for j in 1:length(ExpMonom)
      N[i,j]=coef(d,ExpMonom[j])
    end
  end
  return N
end


##evalPoly: evaluates a polynomial in a point
##input:d: the maximum degree of the polynomial;
##coefs:coefficients vector of the polynomial ; point: the point to be evaluated
##output:evaluation of polynomial with coefficients coefs in point
function evalPoly(d,coefs,point)
  monom=monomiosuptod(d,length(point))
  result=zeros(length(monom))
    for j in 1:length(monom)
        result[j]= prod([point[s]^monom[j][s] for s in 1:length(point)])
    end
  return coefs'*result
end


##solucion: Gives the solution of recovering a measure using aprroximate recovery in Theorem(A)
##input: Y0:Points to be rotated to obtain Y, the support of the measure;measureCoefs: coefficients of the measure
##tau: nu tolerance in Theorem(A);epsilon: rotation angle;
##matrizCoef: matrix of polynomial coefficients;matrizMedicion:measurement matrix
##output: the list [indBig,solBig], indBig:indices of points supports of biggest coefficients;}
##solBig: biggest coefficients
function solucion(Y0,measureCoefs,tau,epsilon,matrizCoef,matrizMedicion)
  Y = []
  for i in 1:length(Y0)
    push!(Y,rotarEp(epsilon,Y0[i]))
  end
  indiOp = []
  for i in 1:length(Y)
    push!(indiOp,puntAprox(Puntos,Y[i]))
  end
  ##b_miu:the moments of miu
  evY=evaluamonomios(mons, Y)
  b_miu =matrizCoef*evY*measureCoefs
  ##optimization problem
  modelo = Model(solver = MosekSolver())
  @variable(modelo, x[1:N])
  @variable(modelo, z[1:N])
  @objective(modelo, Min,sum(z))
  @constraint(modelo,norm(matrizMedicion*x-b_miu)<=tau)
  @constraint(modelo, x.<=z)
  @constraint(modelo, -x.<=z)
  status = solve(modelo)
  sol = getvalue(x)
  solBig = sol[sol.>=0.0001]
  indBig= findin(sol,solBig)
  return [indBig,solBig]
end



##solucion2:Solve approximate recovery in Theorem(A) for a measure supported mu on points in S^n
##input: Y: support of measure mu; measureCoefs: weights of measure mu;
##matrizCoef: matrix of polynomial coefficients; matrizMedicion: matrix of measures
##output: return vector of solution of the problem
function solucion2(Y,measureCoefs,tol,matrizCoef,matrizMedicion)
   ##b_miu:the moments of miu
   evY=evaluamonomios(mons, Y)
   b_miu =matrizCoef*evY*measureCoefs
   modelo = Model(solver = MosekSolver())

   @variable(modelo, x[1:N])
   @variable(modelo, z[1:N])

   @objective(modelo, Min,sum(z))

   @constraint(modelo,norm(matrizMedicion*x-b_miu)<=tol)
   @constraint(modelo, x.<=z)
   @constraint(modelo, -x.<=z)

   status = solve(modelo)
   return getvalue(x)
end

##solucion3:Solve approximate recovery in Theorem(A) for a measure supported mu on points in S^n
##input: Y: support of measure mu; measureCoefs: weights of measure mu;
##matrizCoef: matrix of polynomial coefficients; matrizMedicion: matrix of measures
##output: the list [indBig,solBig], indBig:indices of points supports of biggest coefficients;
##solBig: biggest coefficients
function solucion3(Y,measureCoefs,tol,matrizCoef,matrizMedicion)
   ##b_miu:the moments of miu
   evY=evaluamonomios(mons, Y)
   b_miu =matrizCoef*evY*measureCoefs
   modelo = Model(solver = MosekSolver())

   @variable(modelo, x[1:N])
   @variable(modelo, z[1:N])

   @objective(modelo, Min,sum(z))

   @constraint(modelo,norm(matrizMedicion*x-b_miu)<=tol)
   @constraint(modelo, x.<=z)
   @constraint(modelo, -x.<=z)

   status = solve(modelo)
   sol =  getvalue(x)
   solBig = sol[sol.>=0.0015]
   indBig= findin(sol,solBig)
   return [indBig,solBig]
end


##errorMed: Gives the momentum norm of difference of both measures
##input: Y, points original meadure;weig, weights original measure; Xsol,
##points approximation measure; Csol, coefficients Xsol; d, dregree of polynomials;m number of polynomials.
##output: the momentum norm of difference of both measures, for moments given by m KSS polynomials with degree d
function errorMed(Y, weig, Xsol,Csol, d, m)
  mons = []
     for i in 0:d
       L = monomios(i,n+1)
       mons = append!(mons,L)
     end
  polinomios = coefMonomios(mons,d,m)
  evY = evaluamonomios(mons,Y)
  evX = evaluamonomios(mons, Xsol)
  momOrig = polinomios*evY*weig
  momApr = polinomios*evX*Csol
  return norm(momOrig-momApr)
end

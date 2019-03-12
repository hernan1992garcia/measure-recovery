###used Packages
using Distributions
using MultivariatePolynomials
using JuMP
using PolyJuMP
using SumOfSquares
using Mosek
using Iterators
using DataStructures
using DynamicPolynomials
using Optim
using Cubature
using PyPlot
using Jacobi
#########################
#########################
###linear_system: Calculates coefficients given the support and moments of a measure
##imput:b:moments of the measure, soporte: the support of the measure,functions:
##functions which moments are taken
function linear_system(b,soporte,functions)
    ev_matriz = evalua_polys(soporte,functions)
    coeficientes = ev_matriz\b
    return coeficientes
end
#####################################
tau = 0.0
function CS(b,soporte,functions,tau)
    ev_matriz = evalua_polys(soporte,functions)
    N = length(soporte)
    modelo = Model(solver = MosekSolver())
  @variable(modelo, x[1:N])
  @variable(modelo, z[1:N])
  @constraint(modelo,ev_matriz*x.==b)
  @constraint(modelo, x.<=z)
  @constraint(modelo, -x.<=z)

  @objective(modelo, Min,sum(z))
 # @constraint(modelo,norm(ev_matriz*x-b)<=tau)

  status = solve(modelo)
  sol = getvalue(x)
  return sol
end
#####################################
##clean_measurements_line: clean the noised vector b_noise for a measure on the line
##mons: monomials of quadratic module
##S: domain of the measure to be recovered
##delta: parameter as in the Problem (2)
##b_noise: the noised vector of measurements
##dominio: list of polynomials defining the semialgebraic set
##return: the cleaned vector of measurements
function clean_measurements_line(mons,mons_ort,delta,b_noise,dominio)
    mons = convert(Array{typeof(mons[1])},mons)
    modelo = SOSModel(solver = MosekSolver())

    @variable(modelo, b)
    @variable(modelo, a[1:length(mons_ort)])

    p_minus = 1 - sum([a[i]*mons_ort[i] for i in 1:length(mons_ort)])
    p_plus = 1 + sum([a[i]*mons_ort[i] for i in 1:length(mons_ort)])

    @variable modelo p_m_0 SOSPoly(mons)
    @variable modelo p_p_0 SOSPoly(mons)
    @variable modelo p_m_1 SOSPoly(mons)
    @variable modelo p_p_1 SOSPoly(mons)

    @constraint(modelo,p_m_0 +p_m_1*dominio[1] == p_minus)
    @constraint(modelo,p_p_0 +p_p_1*dominio[1] == p_plus)


    @constraint(modelo, norm(a)<=b)

    @objective(modelo,Max,dot(a,b_noise)-delta*b)

    solve(modelo)
    a_opt = getvalue(a)

    #####################################
    #####################################
    modelo2 = Model(solver = MosekSolver())
    ##l: the values of the operator L over mons
    @variable(modelo2,l[1:length(mons_ort)])

    minus_vec = l - b_noise

    @constraint(modelo2,norm(minus_vec)<=delta)

    @objective(modelo2,Min,dot(a_opt,l))

    solve(modelo2)
    b_clean = getvalue(l)

    return (b_clean,a_opt)
end
#####################################
#####################################
##clean_measurements_square: clean the noised vector b_noise for a measure on the square
##mons: monomials of quadratic module
##S: domain of the measure to be recovered
##delta: parameter as in the Problem (2)
##b_noise: the noised vector of measurements
##dominio: list of polynomials defining the semialgebraic set
##return: the cleaned vector of measurements
function clean_measurements_square(mons,mons_ort,delta,b_noise,dominio)
    mons = convert(Array{typeof(mons[1])},mons)
    modelo = SOSModel(solver = MosekSolver())

    @variable(modelo, b)
    @variable(modelo, a[1:length(mons_ort)])

    p_minus = 1 - sum([a[i]*mons_ort[i] for i in 1:length(mons_ort)])
    p_plus = 1 + sum([a[i]*mons_ort[i] for i in 1:length(mons_ort)])

    @variable modelo p_m_0 SOSPoly(mons)
    @variable modelo p_p_0 SOSPoly(mons)
    @variable modelo p_m_1 SOSPoly(mons)
    @variable modelo p_p_1 SOSPoly(mons)
    @variable modelo p_m_2 SOSPoly(mons)
    @variable modelo p_p_2 SOSPoly(mons)


    @constraint(modelo,p_m_0 +p_m_1*dominio[1] + p_m_2*dominio[2] == p_minus)
    @constraint(modelo,p_p_0 +p_p_1*dominio[1] + p_p_2*dominio[2]== p_plus)


    @constraint(modelo, norm(a)<=b)

    @objective(modelo,Max,dot(a,b_noise)-delta*b)

    solve(modelo)
    a_opt = getvalue(a)

    #####################################
    #####################################
    modelo2 = Model(solver = MosekSolver())
    ##l: the values of the operator L over mons
    @variable(modelo2,l[1:length(mons_ort)])

    minus_vec = l - b_noise

    @constraint(modelo2,norm(minus_vec)<=delta)

    @objective(modelo2,Min,dot(a_opt,l))

    solve(modelo2)
    b_clean = getvalue(l)

    return (b_clean,a_opt)
end
#####################################
#####################################

##integ_v: vector of measurements
##mons:monomials of degree up to 2d
##mons_half:monomials of degree up to d
function error_SOS0(integ_v,mons,mons_half,mons_moments_half)
    prod_mons = Array{Any}(length(mons_half),length(mons_half))

    for i in 1:length(mons_half)
        for j in 1:length(mons_half)
            prod_mons[i,j] = mons_moments_half[i] + mons_moments_half[j]
        end
    end
    ################################
    ################################
    B = zeros(length(mons_half),length(mons_half))

    for i in 1:length(mons_half)
        for j in 1:length(mons_half)
            ind = findin(mons_moments,[prod_mons[i,j]])
            B[i,j] = integ_v[ind[1]]
        end
    end

    l = length(mons_half)
    n = length(integ_v)
    restr_matrix = matrizCoef[n-l+1:n,n-l+1:n]
    modelo = Model(solver = MosekSolver())
    @variable(modelo, A[1:l,1:l])
    @SDconstraint(modelo, A>=0)

    @constraint(modelo,trace(A)==1)
    @objective(modelo,Min,trace(A*B))

    status = solve(modelo)
    val_A = getvalue(A)
    val_err = getvalue(trace(A*B))

    return [val_A, val_err]
end

#####################################
#####################################
###evalua_polys_high: evaluate polynomials in dimension graeater than 1
##input:supp_pun: the list of points to be evaluated, mons_ort: the functions to be evaluated
##output: matrix with columns the evaluation of each point in all the functions mons_ort
function evalua_polys_high(supp_pun,mons_ort)
    ev_mat = zeros(length(mons_ort),length(supp_pun))
    vars = variables(mons_ort[length(mons_ort)])
    for i in 1:length(mons_ort)
        for j in 1:length(supp_pun)
            ev_mat[i,j] = convert(Float64,subs(mons_ort[i],vars => supp_pun[j]))
        end
    end
    return ev_mat
end
#####################################
#####################################
###evalua_polys: evaluate polynomials in dimension 1
##input:supp_pun: the list of points to be evaluated, mons_ort: the functions to be evaluated
##output: matrix with columns the evaluation of each point in all the functions mons_ort
function evalua_polys(supp_pun,mons_ort)
    ev_mat = zeros(length(mons_ort),length(supp_pun))
    vars = MultivariatePolynomials.variables(mons_ort[length(mons_ort)])
    for i in 1:length(mons_ort)
        for j in 1:length(supp_pun)
            if length(supp_pun[j])==1
                ev_mat[i,j] = convert(Float64,MultivariatePolynomials.subs(mons_ort[i],vars => [supp_pun[j]]))
            else
                ev_mat[i,j] = convert(Float64,MultivariatePolynomials.subs(mons_ort[i],vars => supp_pun[j]))
            end
        end
    end
    return ev_mat
end
#####################################
#####################################
##integ_v: vector of measurements
##mons:monomials of degree up to 2d
##mons_half:monomials of degree up to d
function error_SOS_aprox(integ_v,mons_ort,mons_half,matrizCoef)
    l = length(mons_half)
    n = length(integ_v)

    prod_ort = Array{Any}(length(mons_half),length(mons_half))
    mons_ort_half = mons_ort[1:l]

    for i in 1:length(mons_ort_half)
        for j in 1:length(mons_ort_half)
            prod_ij = mons_ort_half[i]*mons_ort_half[j]
            coefs_prod = coefficients(prod_ij,mons0)
            ort_basis_coef = (matrizCoef')\coefs_prod
            integral = dot(ort_basis_coef,integ_v)
            prod_ort[i,j] = integral
        end
    end
    ################################
    modelo = Model(solver = MosekSolver())
    @variable(modelo, A[1:l,1:l])
    @SDconstraint(modelo, A>=0)

    @constraint(modelo,trace(A)==1)
    producto = Array{Any}(l)
    for i in 1:l
        producto[i] = sum([A[i,k]*prod_ort[k,i] for k in 1:l])
    end

    @objective(modelo,Min,sum(producto))

    status = solve(modelo)
    val_A = getvalue(A)
    val_err = getvalue(sum(producto))

    return [val_A, val_err]
end

###############################
##integ_mon: integrates a monomial mon
##input: mon: the monomial, lims_inf: ythe inferior limits of the square;
##lims_sup:superior limits of square
##output: integral of mon over the square
function integ_mon(lims_inf,lims_sup,mon)
    exps_m = exponents(mon)
    if exps_m == zeros(length(exps_m))
        return prod(lims_sup-lims_inf)
    else
        exps_int = exps_m+ones(length(exps_m))
        coefs_int = [1/j for j in exps_int]

        integ = prod([coefs_int[i]*(lims_sup[i]^exps_int[i]- lims_inf[i]^exps_int[i]) for i in 1:length(exps_m)])
        return integ
    end
end

##################################
##integ_poly: integrates poly over square given for lims_inf,lims_sup
##input: lims_inf, lims_sup: inferior and superior limits of box
##poly, the polynomial
##output: the alue of the integral
function integ_poly(lims_inf,lims_sup,poly)
    monos = monomials(poly)
    coefis = coefficients(poly,monos)
    integ = sum([coefis[i]*integ_mon(lims_inf,lims_sup,monos[i]) for i in 1:length(monos)])
    return integ
end
###############################
#############
####projection: projects mon_2 over mon_1
##input:mon_1, mon_2 the monomials
##output: polynomial which is projection of mon_2 over mon_1
function projection(mon_1,mon_2, lim_inf, lim_sup)
    num = integ_poly(lim_inf,lim_sup,mon_1*mon_2)
    denom = integ_poly(lim_inf,lim_sup,mon_1^2)
    proy = (num/denom)*mon_1
    return proy
end
##########################################
##########################################
## ortonormal basis: makes grhm schmidt process to basis mons
##input: mons: te monomial basis,lim_inf: the inferior limit of integration;
##lim_sup: the superior limit of integration
###output: ortonormal basis
function ortonormalization(mons,lim_inf,lim_sup)
    ortonormal_basis = []
    append!(ortonormal_basis,[mons[1]])

    for i in 2:length(mons)
        mon = mons[i]
        restar = []
        for j in 1:(i-1)
            proy = projection(ortonormal_basis[j],mon,lim_inf,lim_sup)
            append!(restar, [proy])
        end
         append!(ortonormal_basis,[mon - sum(restar)])
    end

    for i in 1:length(ortonormal_basis)
        ortonormal_basis[i] = (1/(integ_poly(lim_inf,lim_sup,ortonormal_basis[i]^2))^(1/2))*ortonormal_basis[i]
    end
    return ortonormal_basis
end
##############################################
###############################################
##ev_mons_mas_cheb: evaluates z in the weight \frac{1}{\sqrt{1-x^2}} approximately
##from 0 to 1
function ev_mons_mas_cheb(z)
    ev_vec = [z^degree(m) for m in mons_ort]
    ev_vec = ((1.0000000001-z^2)^(1/2))*ev_vec
    return ev_vec
end

##ev_mons_mas_cheb: evaluates z in the weight \frac{1}{\sqrt{1-x^2}} approximately
##from -1 to 0
function ev_mons_menos_cheb(z)
    ev_vec = [(-z)^degree(m) for m in mons_ort]
    ev_vec = ((1.0000000001-(-z)^2)^(1/2))*ev_vec
    return ev_vec
end

######L:function which gives the coefs value of operator L in deduction of Henrion: [L(v)]_{mons}
function L(v,space_vars,f,mons)
    dv_dt = differentiate(v,t)
    grad_v = differentiate(v,space_vars)
    #if length(space_vars)==1
    #    L_v = dv_dt[1] + dot(grad_v[1],f)
    #else
        L_v = dv_dt + sum([grad_v[i]*f[i] for i in 1:length(space_vars)])
    #end
    return coefficients(L_v,mons)
end

##################################
##################################
###times_g: Gives the matrix of multiplication by g
function times_g(mons,mon,g)
    g_prod = coefficients(g*mon,mons)
    return g_prod
end
##########################################
##########################################
###calculates the moments of the occupation measure of the optimal trajectory of the problem, given its inputs.
###output: val_A the moments vector, opt_value, th optimal value of the functional l
function momentos0(X_I,X_F,T,l,f,domain_rest,d_test,d_l,d_rest,d_f)
    deg_mons = d_test + d_l+d_rest +d_f#max(d_test + d_l,d_rest +d_f)
    deg_mons_half = convert(Int64,floor(deg_mons/2))

    mons =  []
    for i in 0:deg_mons
        append!(mons,monomials(z,i))
    end
    mons = convert(Array{typeof(mons[1])},mons)

    mons_half =  []
    for i in 0:deg_mons_half
        append!(mons_half,monomials(z,i))
    end
    mons_half = convert(Array{typeof(mons[1])},mons_half)
    size_half = length(mons_half)
    mons_test = []
    for i in 0:d_test
        append!(mons_test,monomials([t;x],i))
    end
    mons_test = convert(Array{typeof(mons_test[1])},mons_test)
    ##################################################

    test_final_vec = [convert(Float64,vv([t;x] =>[T[2];X_F[1]])) for vv in mons_test]
    test_init_vec = [convert(Float64,vv([t;x] =>[T[1];X_I[1]])) for vv in mons_test]

    ######################################
    L_matrix = zeros(length(mons),length(mons_test))

    for i in 1:length(mons_test)
        L_matrix[:,i] = L(mons_test[i],x,f,mons)
    end

    ##################
    coef_l = coefficients(l,mons)

    size_mom = length(mons)
    size_test = length(mons_test)
    l_rest = length(domain_rest)
    mons_matrix = mons_half*mons_half'

    modelo = Model(solver = MosekSolver())
    @variable(modelo, A[1:size_mom])

    for m in domain_rest
        g = m
        rest_mat = Array{Any}(size_half,size_half)
        for i in 1:size_half
            for j in 1:size_half
                coef = dot(times_g(mons,mons_matrix[i,j],g),A)
                rest_mat[i,j] = coef
            end
        end
        @SDconstraint(modelo, rest_mat >= 0)
    end

    moments_vec = A

    @constraint(modelo,L_matrix'*moments_vec.== test_final_vec - test_init_vec)

    @objective(modelo,Min,dot(coef_l,moments_vec))

    status = solve(modelo)
    val_A = getvalue(A)
    opt_value = getvalue(dot(coef_l,moments_vec))
    return (val_A,opt_value)
end

#########################################
#########################################
function momentosMix(X_I,X_F,T,l,f,domain_rest,d_test,d_l,d_rest,d_f)
    deg_mons = d_test + d_l+d_rest +d_f#max(d_test + d_l,d_rest +d_f)
    deg_mons_half = convert(Int64,floor(deg_mons/2))

    mons =  []
    for i in 0:deg_mons
        append!(mons,MultivariatePolynomials.monomials(z,i))
    end
    mons = convert(Array{typeof(mons[1])},mons)

    mons_half =  []
    for i in 0:deg_mons_half
        append!(mons_half,monomials(z,i))
    end
    mons_half = convert(Array{typeof(mons[1])},mons_half)
    size_half = length(mons_half)
    mons_test = []
    for i in 0:d_test
        append!(mons_test,monomials([t;x],i))
    end
    mons_test = convert(Array{typeof(mons_test[1])},mons_test)
    ##################################################

    test_final_vec = [convert(Float64,vv([t;x] =>[T[2];X_F[1]])) for vv in mons_test]
    test_init_vec = [convert(Float64,vv([t;x] =>[T[1];X_I[1]])) for vv in mons_test]

    ######################################
    L_matrix = zeros(length(mons),length(mons_test))

    for i in 1:length(mons_test)
        L_matrix[:,i] = L(mons_test[i],x,f,mons)
    end

    ##################
    coef_l = coefficients(l,mons)

    size_mom = length(mons)
    size_test = length(mons_test)
    l_rest = length(domain_rest)
    mons_matrix = mons_half*mons_half'

    modelo = Model(solver = MosekSolver())
    @variable(modelo, A[1:size_mom])
    @variable(modelo, w[1:length(domain_rest)*size_half,1:size_half])

    for k in 1:length(domain_rest)
        g = domain_rest[k]
        rest_mat = Array{Any}(size_half,size_half)
        for i in 1:size_half
            for j in 1:size_half
                coef = dot(times_g(mons,mons_matrix[i,j],g),A)
                rest_mat[i,j] = coef
            end
        end
        abs_mat = w[((k-1)*size_half +1):k*size_half,1:size_half]
######Linear restrictions for diagonal dominant matrix
        @constraint(modelo, -abs_mat.<=rest_mat)
        @constraint(modelo, rest_mat.<=abs_mat)

        for l in 1:size_half
            @constraint(modelo, rest_mat[l,l]>=0)
            @constraint(modelo,abs_mat[l,l]>= sum(abs_mat[l,:]-abs_mat[l,l]))
        end

    end

   rest_mat0 = Array{Any}(size_half,size_half)
    g = 1
    for i in 1:size_half
        for j in 1:size_half
            coef = dot(times_g(mons,mons_matrix[i,j],g),A)
            rest_mat0[i,j] = coef
        end
    end
    @SDconstraint(modelo, rest_mat0 >= 0)

    moments_vec = A

    @constraint(modelo,L_matrix'*moments_vec.== test_final_vec - test_init_vec)

    @objective(modelo,Min,dot(coef_l,moments_vec)+sum(w))

    status = solve(modelo)
    val_A = getvalue(A)
    opt_value = getvalue(dot(coef_l,moments_vec))
    return (val_A,opt_value)
end

####################################
#####################################
#####################################
#ort_basis: Change of basis matrix

function momentosDSOS(X_I,X_F,T,l,f,domain_rest,d_test,d_l,d_rest,d_f,ort_basis)
    deg_mons = d_test + d_l+d_rest +d_f#max(d_test + d_l,d_rest +d_f)
    deg_mons_half = convert(Int64,floor(deg_mons/2))

    mons0 =  []
    for i in 0:deg_mons
        append!(mons0,MultivariatePolynomials.monomials(z,i))
    end
    mons0 = convert(Array{typeof(mons0[1])},mons0)
    mons = ort_basis*mons0

    size_half = convert(Int64,factorial(deg_mons_half + length(z))/(factorial(deg_mons_half)*factorial(length(z))))
    mons_half = mons[1:size_half]

    mons_test = []
    for i in 0:d_test
        append!(mons_test,monomials([t;x],i))
    end
    mons_test = convert(Array{typeof(mons_test[1])},mons_test)
    ##################################################

    test_final_vec = [convert(Float64,vv([t;x] =>[T[2];X_F[1]])) for vv in mons_test]
    test_init_vec = [convert(Float64,vv([t;x] =>[T[1];X_I[1]])) for vv in mons_test]

    ######################################
    L_matrix = zeros(length(mons),length(mons_test))

    for i in 1:length(mons_test)
        L_matrix[:,i] = ort_basis*L(mons_test[i],x,f,mons0)
    end

    ##################
    coef_l = ort_basis*coefficients(l,mons0)

    size_mom = length(mons)
    size_test = length(mons_test)
    l_rest = length(domain_rest)
    mons_matrix = mons_half*mons_half'

    modelo = Model(solver = MosekSolver())
    ###A: vector of integrals of polynomials in mons respect the desired measure
    @variable(modelo, A[1:size_mom])

    @variable(modelo, w[1:length(domain_rest)*size_half,1:size_half])

    for k in 1:length(domain_rest)
        g = domain_rest[k]
        rest_mat = Array{Any}(size_half,size_half)
        for i in 1:size_half
            for j in 1:size_half
                coef = dot(ort_basis*times_g(mons0,mons_matrix[i,j],g),A)
                rest_mat[i,j] = coef
            end
        end
#@constraint modelo mons_half'*rest_mat*mons_half in DSOSCone()
        abs_mat = w[((k-1)*size_half +1):k*size_half,1:size_half]
######Linear restrictions for diagonal dominant matrix
        @constraint(modelo, -abs_mat.<=rest_mat)
        @constraint(modelo, rest_mat.<=abs_mat)

        for l in 1:size_half
            @constraint(modelo, rest_mat[l,l]>=0)
            @constraint(modelo,abs_mat[l,l]>= sum(abs_mat[l,:]-abs_mat[l,l]))
        end

    end

    moments_vec = A

    @constraint(modelo,L_matrix'*moments_vec.== test_final_vec - test_init_vec)

    @objective(modelo,Min,dot(coef_l,moments_vec)+sum(w))

    status = solve(modelo)
    val_A = getvalue(A)
    opt_value = getvalue(dot(coef_l,moments_vec))
    return (val_A,opt_value)
end

################################################
################################################
function momentos_r_DSOS(X_I,X_F,T,l,f,domain_rest,d_test,d_l,d_rest,d_f,r)

    deg_mons = d_test + d_l+d_rest +d_f + 2*r
    deg_mons_half = convert(Int64,floor(deg_mons/2))

    mons =  []
    for i in 0:deg_mons
        append!(mons,MultivariatePolynomials.monomials(z,i))
    end
    mons = convert(Array{typeof(mons0[1])},mons0)

    size_half = convert(Int64,factorial(deg_mons_half + length(z))/(factorial(deg_mons_half)*factorial(length(z))))
    mons_half = mons[1:size_half]

    mons_test = []
    for i in 0:d_test
        append!(mons_test,monomials([t;x],i))
    end
    mons_test = convert(Array{typeof(mons_test[1])},mons_test)
    ##################################################

    test_final_vec = [convert(Float64,vv([t;x] =>[T[2];X_F[1]])) for vv in mons_test]
    test_init_vec = [convert(Float64,vv([t;x] =>[T[1];X_I[1]])) for vv in mons_test]

    ######################################
    L_matrix = zeros(length(mons),length(mons_test))

    for i in 1:length(mons_test)
        L_matrix[:,i] = L(mons_test[i],x,f,mons)
    end

    ##################
    coef_l = coefficients(l,mons)

    size_mom = length(mons)
    size_test = length(mons_test)
    l_rest = length(domain_rest)
    mons_matrix = mons_half*mons_half'

    modelo = Model(solver = MosekSolver())
    ###A: vector of integrals of polynomials in mons respect the desired measure
    @variable(modelo, A[1:size_mom])

    @variable(modelo, w[1:length(domain_rest)*size_half,1:size_half])

    norma = (sum([z[i]^2 for i in 1:length(z)]))^r

    for k in 1:length(domain_rest)
        g = norma*domain_rest[k]
        rest_mat = Array{Any}(size_half,size_half)
        for i in 1:size_half
            for j in 1:size_half
                coef = dot(times_g(mons,mons_matrix[i,j],g),A)
                rest_mat[i,j] = coef
            end
        end
#@constraint modelo mons_half'*rest_mat*mons_half in DSOSCone()
        abs_mat = w[((k-1)*size_half +1):k*size_half,1:size_half]
######Linear restrictions for diagonal dominant matrix
        @constraint(modelo, -abs_mat.<=rest_mat)
        @constraint(modelo, rest_mat.<=abs_mat)

        for l in 1:size_half
            @constraint(modelo, rest_mat[l,l]>=0)
            @constraint(modelo,abs_mat[l,l]>= sum(abs_mat[l,:]-abs_mat[l,l]))
        end

    end

    moments_vec = A

    @constraint(modelo,L_matrix'*moments_vec.== test_final_vec - test_init_vec)

    @objective(modelo,Min,dot(coef_l,moments_vec)+sum(w))

    status = solve(modelo)
    val_A = getvalue(A)
    opt_value = getvalue(dot(coef_l,moments_vec))
    return (val_A,opt_value)
end

################################################
################################################

##Generates a random matrix of size tam x tam with Haar measure("see what is a random matrix":Diaconis)
function ort_random_matrix(tam)
    dist=Normal(0.0,1)
    base0 = zeros((tam,tam))
    for i in 1:tam
        base0[i,:]=rand(dist,tam)
    end
    base_ort = gram_schmidt(base0)
    return base_ort
end

#################
##################
##################
##gram_schmidt: given an n x n matrix generates a matrix whose columns
##are gram schmidt ortonormalization process of original columns

function gram_schmidt(mat)
    ort_base = zeros(size(mat))
    base1 = [(1/norm(mat[:,1]))*mat[:,1]]
    ort_base[:,1] = base1[1]
    for i in 2:length(mat[1,:])
        w = mat[:,i]
        proy = zeros(length(mat[1,:]))
        for v in base1
            proy = proy + dot(v,w)*v
        end
        ort = w - proy
        append!(base1,[(1/norm(ort))*ort])
        ort_base[:,i] = base1[i]
    end
    return ort_base
end

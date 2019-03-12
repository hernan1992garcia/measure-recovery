function log_q_form(p)
    ev_basis_p = evalua_polys([p],mons[1:length(mons_half)])
    res = log.((ev_basis_p)'*A*ev_basis_p)
    return res[1]
end

####################################
function approx_grad(log_q_form,p,h)
    direc = eye(length(p))
    grad = zeros(length(p))
    for i in 1:length(p)
      val_mas = log_q_form(p+h*direc[i,:])
      val_menos = log_q_form(p-h*direc[i,:])
      #val_menos = log_q_form(p)
      grad[i] = ((2*h)^(-1))*(val_mas - val_menos)
      #grad[i] = ((h)^(-1))*(val_mas - val_menos)
    end
    return grad
end

####################################
function approx_grad_maurostyle(funcion,p,q,h)
    direc = eye(length(q))
    grad = zeros(length(q))
    for i in 1:length(q)
      val_mas = funcion(p,q+h*direc[i,:])
      val_menos = funcion(p,q-h*direc[i,:])
      #val_menos = log_q_form(p)
      grad[i] = ((2*h)^(-1))*(val_mas - val_menos)
      #grad[i] = ((h)^(-1))*(val_mas - val_menos)
    end
    return grad
end
#############################
#############################
p = [0.5,0.5,0.5]
function fun(p,q)
    return log_q_form(q)- (dot(approx_grad(log_q_form,p,h),q-p) + log_q_form(p))
end
##newton_step: calculates step of multivariate approximated newton method
##input: fun: function to fin zero, q: beggining point, p:the point defining the tangent plane(this is fixed!!)
##output: point wich intersects z = 0 and minimizes distance with p
q = [0.2,0.3,0.4]
function newton_step(fun,p,q,h)
    grad = approx_grad_maurostyle(fun,p,q,h)

    val_form = fun(p,q)
    inner = dot(grad,q)

    b_lin = -(val_form + inner)

    modelo = Model(solver = MosekSolver())
    @variable(modelo, w[1:length(grad)])
    @variable(modelo, z[1:length(grad)])

    @constraint(modelo, dot(grad,w) == -b_lin)
    @constraint(modelo, 0.<=z)
    @constraint(modelo, w-p.<=z)
    @constraint(modelo, -w+p.<=z)

    @objective(modelo, Min, sum(z))

    status = solve(modelo)
    val_x = getvalue(w)

    return val_x
end

######################
p = [0.3,0.1,0.2]
q = [0.1,0.1,0.1]


q = p + 0.01*approx_grad_maurostyle(fun,p,q,h)
p = newton_step(fun,p,q,h)


print(p)
print(fun(p,q))
print(log_q_form(p))

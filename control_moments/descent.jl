
###################################
#function log_q_form(p)
#    is_neg = 1*(0.>[convert(Float64,subs(m,z=>p)) for m in domain_rest])
#    if sum(is_neg) > 0
#        res = 1
#        return res
#    else
#        ev_basis_p = evalua_polys([p],mons[1:length(mons_half)])
#        res = log((ev_basis_p)'*A*ev_basis_p)
#        return res[1]
#    end
#end
function log_q_form(p)
    ev_basis_p = evalua_polys([p],mons[1:length(mons_half)])
    res = log.(1 + (ev_basis_p)'*A*ev_basis_p )#(abs.((ev_basis_p)'*A*ev_basis_p ))
    return res[1]
end

####################################
####################################
##avoid_dir: integer list [i_1,...,i_n]. derivative in direction e_{i_k} is not considered
function approx_grad(log_q_form,p,h,avoid_dir)
    direc = eye(length(p))
    grad = zeros(length(p))
    for i in 1:length(p)
        if i in avoid_dir
            grad[i] = 0
        else
          val_mas = log_q_form(p+h*direc[i,:])
          val_menos = log_q_form(p-h*direc[i,:])
          #val_menos = log_q_form(p)
          grad[i] = ((2*h)^(-1))*(val_mas - val_menos)
          #grad[i] = ((h)^(-1))*(val_mas - val_menos)
        end
    end
    return grad
end

####################################
####################################
function descent_step(log_q_form,p,h,nu,avoid_dir)
    grad = approx_grad(log_q_form,p,h,avoid_dir)
    return p-nu*grad
end

########################
########################
function choose_limit_pr(box_l,box_u,p,p0)
    proy = zeros(length(box_l))
    for i in 1:length(box_l)
        if p[i]<=box_l[i]
            proy[i] = box_l[i]
        elseif box_l[i] <= p[i] <= box_u[i]
            proy[i] = p[i]
        else
            proy[i] = box_u[i]
        end
    end
    return proy
end
####################################
####################################
##box_l: inferior limits of the box
##box_u: superior limits of the box
function descent_step_box(log_q_form,p,h,nu,box_l,box_u,avoid_dir)
    grad = approx_grad(log_q_form,p,h,avoid_dir)
    step = p-nu*grad
    p0 = p
    if sum(box_l.<= step .<= box_u)==length(step)
        return step
    else
        return choose_limit_pr(box_l,box_u,step,p0)
    end
end
#####################################
#####################################

function backtracking(log_q_form,p,grad_f,alph,bet)
    t = 1
    while log_q_form(p -t*grad_f) > log_q_form(p) + alph*t*norm(grad_f)^2
        t = bet*t
    end
    return t
end
#####################################
#####################################
function descenso(log_q_form,p,h,steps,avoid_dir)
    p_opt = p
    for i in 1:steps
        alph = 0.1
        bet = 0.5
        grad_f = approx_grad(log_q_form,p_opt,h,avoid_dir)
        nu = backtracking(log_q_form,p_opt,grad_f,alph,bet)
        p_opt = descent_step(log_q_form,p_opt,h,nu,avoid_dir)
    end
    return p_opt
end
#####################################
#####################################

function descenso_proy(log_q_form,p,h,box_l,box_u,steps,avoid_dir)
    p_opt = p
    i = 1
    while i <= steps
        alph = 0.4
        bet = 0.8
        grad_f = approx_grad(log_q_form,p_opt,h,avoid_dir)
        nu = backtracking(log_q_form,p_opt,grad_f,alph,bet)
        p_opt = descent_step_box(log_q_form,p_opt,h,nu,box_l,box_u,avoid_dir)
        i = i+1
    end
    return p_opt
end
###################################
#mons_ort = mons
#d = descenso(log_q_form,p,h,nu_0,steps)
#plot(d[1],d[2],marker = "*",color = :red)

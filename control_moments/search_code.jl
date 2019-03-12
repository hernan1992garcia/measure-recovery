
##################################
###########Construct the code
##################################

X_I = [[0;1]]
X_F = [[0; exp(-2*pi)]]
T = [0,2*pi]

dist_tr1 = Uniform(-1,1)
dist_tr2 = Uniform(-1,1)

supp = []
tiem = 0.0:0.1:2*pi
for t_p in tiem
    th_code2 = []
    for i in 1:15
        #append!(th_code2,[[rand(dist_t),rand(dist_tr1),rand(dist_tr2)]])
        append!(th_code2,[[t_p,rand(dist_tr1),rand(dist_tr2)]])
    end

    ###################################
    h = 0.000001
    steps = 200
    #box_l = [0,0,-2*pi]
    #box_u = [t_F,2*pi,2*pi]

    box_l = [0,-1,-1]
    box_u = [2*pi,1,1]

    #box_l = [0,-10,0,-10,-10,0]
    #box_u = [2*pi,10,10*pi,100,10,10]
    avoid_dir = [1]
    #th_code2 = optimas
    optimas = []
    for p in th_code2
        desc = descenso_proy(log_q_form,p,h,box_l,box_u,steps,avoid_dir)#descenso(log_q_form,p,h,steps,avoid_dir)#
        append!(optimas,[desc])
    end
    #####################################
    ###################################
evals = [log_q_form(p) for p in optimas]
optim = optimas[findin(evals,minimum(evals))][1]
append!(supp,[optim])

#    for p in th_code2
#       plot(p[1],p[2],marker = ".",color = :black)
    #    plot(p[1],p[2],marker = ".",color = :blue)
#    end

#    for p in optimas
        #plot(optim[1],optim[2],marker = "+",color = :green)
#        plot(p[1],p[2],marker = ".",color = :red)
    #    plot(p[1],p[2],marker = ".",color = :blue)
#    end
print(t_p,"##",optim[2])
end

x_c = [p[1] for p in supp]
y_c = [p[2] for p in supp]


plot(x_c,y_c,color =:purple)

plot(tiem,y_c)

for p in supp
    plot(p[1],p[2],marker = "o",color = :black)
end
#####################
posy = zeros(length(tiem))
for i in 1:length(tiem)
    t_p = tiem[i]
    pos = -5*t_p^2 + t_p
    posy[i] = exp(-2*t_p)*sin(2*t_p)
#    plot(t_p,pos,marker = ".",color = :yellow)
end

plot(tiem,posy,color =:black)
##############Coefficients of measure
coefes = linear_system(integ_vec,supp,mons)

ev_matriz = evalua_polys(supp,mons)
er = norm(ev_matriz*coefes-integ_vec)

a = [p[3] for p in supp]
##########################x_c = [p[1] for p in supp]
y_c = [p[2] for p in supp]


plot(x_c,y_c,color =:blue)

##########################
using Clustering

optimas_matrix = zeros(length(th_code2[1]),length(optimas))

for i in 1:length(optimas)
    optimas_matrix[:,i] = optimas[i]
end

centroides = kmeans(optimas_matrix,5)
centros = centroides.centers


vals_centros = [log_q_form(optimas[i]) for i in 1:length(optimas)]

centros_matrix = zeros(1,length(optimas))
for i in 1:length(optimas)
    centros_matrix[:,i] = vals_centros[i]
end
kmeans(centros_matrix,2)
##########################
##########################
for p in th_code2
   plot(p[1],p[2],marker = ".",color = :black)
#    plot(p[1],p[2],marker = ".",color = :blue)
end


t_p = 1
mons_ort = mons
value_m = zeros(length(th_tr),length(th_cont))
#value_m = zeros(length(th_tr),length(th_time))
#for u in 1:length(th_time)
   for w in 1:length(th_tr)
    for v in 1:length(th_cont)
        #punto = [th_time[v],th_tr[w],1]
        punto = [t_p,th_tr[w],th_cont[v]]
        value_m[w,v] = log_q_form(punto)
    end
    end
#end

pcolormesh(th_tr,th_cont,value_m')#log.(abs.((value_m)')))
colorbar()
xlabel("x")
ylabel("u")
#######################
mons_ort = mons
value_m = zeros(length(th_tr),length(th_time))
#value_m = zeros(length(th_tr),length(th_time))
#for u in 1:length(th_time)
   for w in 1:length(th_tr)
    for v in 1:length(th_time)
        #punto = [th_time[v],th_tr[w],1]
        punto = [th_tr[w],th_time[v],1]
        value_m[w,v] = log_q_form(punto)
    end
    end
pcolormesh(th_time,th_tr,value_m')#log.(abs.((value_m)')))
colorbar()
xlabel("t")
ylabel("x")
#######################


for p in th_code
    plot(p[1],p[2],marker = "o",color = :red)
end

for p in optimas
    plot(p[1],p[2],marker = "*",color = :orange)
end

for p in punt_supp
    plot(p[1],p[2],marker = ".",color = :white)
end

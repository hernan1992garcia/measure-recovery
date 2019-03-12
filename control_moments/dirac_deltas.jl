###trans_legendre:translate Legendre polynomials from [-1,1] to [0,b]
###input: x: the variable, n: the degree,a,b: the limits of translation interval
###output: translated legendre polynomial
function trans_legendre(x,n,b)
    p = legendre(x,n)
    m = 2/b
    es = sqrt((2*n+1)/2)
    q = es*sqrt(2/b)*subs(p,x => m*x-1)
    return q
end

###dirac_legendre_n:Gives the trncated expansion of the dirac delta centered at y
function dirac_legendre_n(x,b,y,n)
    polys = [trans_legendre(x,i,b) for i in 0:n]
    coefs = [polys[i+1](y) for i in 0:n]
    delta_approx = polys'coefs
    return delta_approx
end


delta_approx = dirac_legendre_n(t,3,1,10)

quadgk(delta_approx,0,3)

n = 4
b = 7
quadgk(trans_legendre(t,n+4,b)*trans_legendre(t,n,b),0,b)

xx = collect(-5:0.1:5)
yy = [delta_approx(r) for r in xx]

plot(xx,yy)

##--Grid of N points in S^1

##equidist: generates N ewuidistant points in S^1
##input: N: number of points
##outpu: Punt: the list of N equidistant points in S^1
function equidist(N)
  Punt = []
  for i in 1:N
      push!(Punt,apoint(2*(i-1)*pi/N))
  end
  return Punt
end

N = 200
Puntos = equidist(N)

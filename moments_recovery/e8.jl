##----Grid of normalized primitive roots of E8(See Section 6)

using Iterators
####
####Constructing A= simple roots of E8
  ind = collect(subsets(1:8,2))

  A = []
  for i in ind
    vec1 = zeros(8)
    vec2 = zeros(8)
    vec3 = zeros(8)
    vec4 = zeros(8)

    vec1[i] = [2; 2]
    vec2[i] = [-2; 2]
    vec3[i] = [-2; -2]
    vec4[i] = [2; -2]

    push!(A,vec1)
    push!(A,vec2)
    push!(A,vec3)
    push!(A,vec4)
  end

  ind2 = [collect(subsets(1:8,2)),collect(subsets(1:8,4)),collect(subsets(1:8,6)),collect(subsets(1:8,8))]

  indP = []
  for j in ind2
    indP = append!(indP,j)
  end

  for l in indP
    vec = ones(8)
    for i in l
      vec[i] = -1
    end
    push!(A,vec)
  end

  push!(A,ones(8))

  Puntos = [x/norm(x) for x in A]

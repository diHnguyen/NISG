using Combinatorics
function EnumX(arcs,b,numArcs,d) 
    X = findall(d.<=b)
    Len = length(X)
    X = combinations(X) |> collect
#     redunX = []
#     global 
    feasX = []
#     println(X)
    #Eliminate combinations exceeding budget
    for i in X
        cost = sum(d[x] for x in i)
#         println("X = ",i,": ",cost)
        if cost <= b
            push!(feasX, i)
        end
    end
#     println("feasible Budget = ", feasX)
#     global 
    subset = true
    
    #Eliminate non-maximally packed x-solutions
    while subset == true
        subset = false
        for i in feasX
            for j in feasX
                if i != j
                    if issubset(i, j) == true
                        subset = true
                        if length(i) < length(j)
                            filter!(x->x != i, feasX)
                        else
                            filter!(x->x != j, feasX)
                        end
                    end
                end
            end
        end    
    end
#     println("feasX = ", feasX)
    return feasX
end
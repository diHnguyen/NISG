function findQgivenP(M_set, numY, Q, path)
    global numArcs, q
    und_q = 1 .- q
    coef_q = zeros(numY)
    
    for i = 1:numY
        M = M_set[i]
        McapP = intersect(M,path)
        if isempty(McapP)
            coef_q[i] = 1
        else
            coef_q[i] = prod(und_q[McapP])
        end      
    end
    push!(Q, coef_q)
    return coef_q, Q
end
function findQgivenM(P_set, numPaths, Q, M)
    global numArcs, q
    und_q = 1 .- q
#     coef_q = 
    M_q = zeros(numPaths)
    
    for i = 1:numPaths
        coef_q = Q[i]
        path = P_set[i]
#         for path in P_set #i = 1:numY
            #M = M_set[i]
        McapP = intersect(M,path)
        if isempty(McapP)
            push!(coef_q, 1)
            push!(M_q,1)
        else
            push!(coef_q, prod(und_q[McapP]))
            push!(M_q,1)
        end      
#         end
        Q[i] = coef_q #hcat(Q, coef_q)
    end
#     return coef_q, Q
    return Q
end
# function getM(numY)
    
#     M_set = [[2, 6, 7], [6, 7, 8], [6, 7, 9], [1, 3, 5, 6, 7]] 
#     if numY > length(M_set)
#         a = M_set[numY-1]
#     else
#         a = M_set[numY]
#     end
#     return a
# end
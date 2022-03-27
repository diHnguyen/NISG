function findf_MP(arcsinPath, M)
    global q,p
    McapP = intersect(arcsinPath, M)
    PminusM = setdiff(arcsinPath, M)
#     println("McapP ", McapP)
#     println("PminusM ", PminusM)
    f_MP = 1
#     println("McapP = ", McapP)
    if isempty(McapP) == false
        f_MP = f_MP*prod((1-q[a]) for a in McapP)
    end
    if isempty(PminusM) == false
        f_MP = f_MP*prod((1-p[a]) for a in PminusM)
    end
    return f_MP
end
function updateF(Q, P_set, M_set, numY, numPaths)
    global q,p
    Q_length = length(Q)
    
    #Except for the initialization, updateQ is called when either M_set or P_set is modified, not BOTH
    if numPaths > Q_length #If P_set is modified
#         println("Adding new Paths to Q")
        for i = (Q_length+1):numPaths
            arcsinPath = P_set[i]
            f_P = zeros(numY)
            for m = 1:numY
                M = M_set[m]
                f_P[m] = findf_MP(arcsinPath,M)
            end
            push!(Q, f_P)
        end
    else #If M_set is modified, i.e., P_set is the same
#         println("Revising existing P-vector")
#         println("numPaths = ", numPaths)
        for i = 1:numPaths
#             println("Q[path] = ", Q[i])
            arcsinPath = P_set[i]
            f_P = Q[i]
            old_M_length = length(f_P)
            for m = (old_M_length+1):numY
#                 println("\tm = ", m)
                M = M_set[m]
                f_MP = findf_MP(arcsinPath,M)
                push!(Q[i], f_MP)
            end
        end
    end
    return Q
end
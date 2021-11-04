# y = [0.0, 0.0, 0, 1]
function findPathsExactPr(arcsinPath, M_set,numY, y)
#     println("y = ", y)
    
#     global feasY, numY , arcsinPath
    
#     feasY = [[2, 6, 7], [6, 7, 8], [6, 7, 9], [1, 3, 5, 6, 7]]
#     numY = length(y)
#     global arcs = [1 4; 1 5; 1 6; 2 6; 3 6; 1 2; 1 3; 4 6; 5 6]
#     global d_x = [2, 5, 2, 7, 1, 0, 0, 5, 5]
#     global q = [0.8, 0.4, 0.6, 0.9, 0.4, 0.6, 0.1, 0.4, 0.9]
#     global b_x = 5
#     global origin = 1
#     global destination = 6
    #P = [[1,8],[2,9],[3],[4,6],[5,7]]
#     p_current = 0
    #Calculating the exact probability that the evader is not detected on path P:
#     println("Exact Pr of being undetected: ")
#     p_exact = [] # zeros(length(P))
#     for path in P
        p = 0
        for m = 1:numY
            p_y = 1
    #         M_in_P = false
            for a in arcsinPath
                if a in M_set[m]
                    p_y = p_y*(1-q[a])
    #                 M_in_P = true
                end
            end
            p = p+ p_y*y[m]
#             println("")
        end
#         push!(p_exact, p)
        println("Path ", arcsinPath, ": ", p)
#         if path == 
#     end
    return p #p_current
end
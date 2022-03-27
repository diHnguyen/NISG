function getCost_SP(x_now, y_now)#, edge, q, numArcs, M_set, numY, R)
    global iter, numArcs, arcs, R, numY, M_set
    
    c_g = zeros(numArcs) #This should start at 0 not 1, the additional 1 units from H-definition only gets added ONCE per PATH not per arc
    interdicted_Arcs = findall(x_now.==1)
    for a = 1:numArcs
        if a in interdicted_Arcs
            c_g[a] = -R 
        else
            c_g[a] = log(1-p[a])
        end
    end
    
#     for a = 1:numArcs
    monitor = findall(y_now .>0)
#     for m = 1:numY
    for m in monitor
#         println("M = ", M_set[m])
        for a in M_set[m]
#             println("a = ",a)
            c_g[a] = c_g[a] + (log(1-q[a]) - log(1-p[a]))*y_now[m]
#             println(a,"c_g[a] = ",c_g[a])
        end
    end
    c_g = c_g.*(-1)
#     y_pos = findall(y_now.>0)
#     for a in y_pos
#         println(M_set[a], "\t", y_now[a])
#     end
#     for a =1:numArcs
#         if c_g[a] > 0 
#             println(a, "\tinterdict? ", x_now[a], "\tc_g ", c_g[a])
#         end
#     end
#     if iter == 1721
#         println("y_now = ", y_now)
#         println("c_g = " ,c_g[monitor])
#     end
    return c_g
end

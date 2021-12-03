function getCost_SP(x_now, y_now, edge, q, numArcs, M_set, numY, R)
    
    c_g = zeros(numArcs) #This should start at 0 not 1, the additional 1 units from H-definition only gets added ONCE per PATH not per arc
    
    interdicted_Arcs = findall(x_now.==1)
    for a in interdicted_Arcs
       c_g[a] = R 
    end
    
#     for a = 1:numArcs
    for m = 1:numY
        for a in M_set[m]
            c_g[a] = c_g[a] - log10(1-q[a])*y_now[m]
        end
    end
#     println(c_g)
    return c_g
end

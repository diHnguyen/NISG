function ConstructQD(arcs,v,C, paths)
#     global arcs
#     global v
#     global C
#     global paths

    Q = ones((length(arcs[:,1]), length(paths))).*v
    D = zeros((length(arcs[:,1]), length(paths)))
    
    for l = 1:length(paths)
#         println("PATH ",l)
        for k = 1:length(paths[l])-1
            arc_on_path_i = paths[l][k]
            arc_on_path_j = paths[l][k+1]
            arcSet = findall(arcs[:,1].==arc_on_path_i)
#             println("Arc :", arc_on_path_i, arc_on_path_j)
            for i in arcSet   
                if arcs[i,2] == arc_on_path_j
                    Q[i,l] = -C 
                    D[i,l] = 1
#                     println("Element ", i,"x", l)
                end
                
            end
        end
    end
#     for k = 1:arcNum
#        println(Q[k,:]) 
#     end
    return Q,D
end
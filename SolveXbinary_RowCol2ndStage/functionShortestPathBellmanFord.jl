function shortestPath_BellmanFord(x_now, y_now, edge)
    
    #println(f,"Current CELL's LB = ", c_L)
    #println(f,"Current CELL's UB = ", c_U)
#     println("Current interdiction x = ", x_now)
    global q, numArcs, R,M_set
    
    c_g = zeros(numArcs)
#     println
    #Get arc cost
    for a = 1:numArcs
        if x_now[a] == 1
            c_g[a] = 10^6
        else
            for m = 1:numY
                if a in M_set[m]
                    c_g[a] = c_g[a] - log(1-q[a])*y_now[m]
                end
            end
#         if x_now[a] > 0
        
        end
    end
    
#     println("c_g = ", c_g)
#     c_g = c_g + (p_y+1)*y[m]
    
    start_node = edge[:,1]
    end_node = edge[:,2]

    no_node = max(maximum(start_node), maximum(end_node) )
    no_link = length(start_node)


    function getShortestX(state, start_node, end_node, origin, destination)
        _x = zeros(Int, length(start_node))
        _path = enumerate_paths(state, destination)

        for i=1:length(_path)-1
            _start = _path[i]
            _end = _path[i+1]

            for j=1:length(start_node)
                if start_node[j]==_start && end_node[j]==_end
                _x[j] = 1
                break
                end
            end

        end
        _x
    end


    graph = SimpleDiGraph(no_node)
    distmx = Inf*ones(no_node, no_node)

    # Adding links to the graph
    for i=1:no_link
        add_edge!(graph, start_node[i], end_node[i])
        distmx[start_node[i], end_node[i]] = c_g[i]
    end

#     Run Dijkstra's Algorithm from the origin node to all nodes
#     state = dijkstra_shortest_paths(graph, origin, distmx)
    
#     Run Bellman-Ford due to negative weights
    state = bellman_ford_shortest_paths(graph, origin, distmx)
    
    label = state.dists
    pred = state.parents
    b_arc = ""
    
    for i = 1: length(state.parents)
        if state.parents[i] != 0 
            b_arc = string(b_arc, "(", state.parents[i], ",", i, ")")
        end
    end
    
    # Retrieving the shortest path
    path = enumerate_paths(state, destination)
    
    #parents = LightGraphs.DijkstraState(state, destination)

    # Retrieving the 'x' variable in a 0-1 vector
    y = getShortestX(state, start_node, end_node, origin, destination)
    #println(f,"y vector:", y)
    
    gx = sum(c_g[i]*y[i] for i = 1:no_link)    
#     SP = sum(c[i]*y[i] for i = 1:length(c))
    #If we want to return a shortest tree
#     T = Int64[]
#     for i = 1:Len
#         if pred[edge[i,2]] == edge[i,1]
#             push!(T, 1)
#         else
#             push!(T, 0)
#         end
#     end

    #println("Minimum Spanning Tree: ", T)
    #println("Shortest path y = ", y)
    #println("Indices of shortest path (edge) = ", y_index)
    #println("Nodes visited =", path)
    #println("Minimum Spanning Tree =", pred)
    #println("Node label =" ,label)
#     println("path ", findall(y.>0)," cost ", gx)
    return y, gx #, SP, T, pred, label, path
end

function gx_bound(edge, A, c, origin, x_now, T, scenNum)
    println(1)
    #println(f,"Current CELL's LB = ", c_L)
    #println(f,"Current CELL's UB = ", c_U)
    #println(f, "Current interdiction x = ", x_now)
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


    graph = Graph(no_node)
    distmx = Inf*ones(no_node, no_node)

    # Adding links to the graph
    for i=1:no_link
        add_edge!(graph, start_node[i], end_node[i])
        distmx[start_node[i], end_node[i]] = c[i]
    end

    # Run Dijkstra's Algorithm from the origin node to all nodes
    state = dijkstra_shortest_paths(graph, origin, distmx)
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
    
    y_Set = ones((scenNum, A))#Array{Int64, (scenNum, A)}
#     y_Set[1,:] = y
    for i = 2:scenNum
        y_Set[1,:] = getShortestX(state, start_node, end_node, origin, T[i])
    end
    #println(f,"y vector:", y)
    
    gx = sum(c[i]*y[i] for i = 1:no_link)    
#     SP = sum(c[i]*y[i] for i = 1:length(c))
    T = zeros(Int64,A)
    for i = 1:A
        if pred[edge[i,2]] == edge[i,1]
#             push!(T, 1)
            T[i]=1
#         else
#             push!(T, 0)
        end
    end

    y_index = findall(y .== 1)

    #println("Minimum Spanning Tree: ", T)
    #println("Shortest path y = ", y)
    #println("Indices of shortest path (edge) = ", y_index)
    #println("Nodes visited =", path)
    #println("Minimum Spanning Tree =", pred)
    #println("Node label =" ,label)

    return y, gx, T, path, y_Set
end

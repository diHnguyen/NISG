function IP_RowGen(x_sol, y_sol)
    global P_set, M_set, numPaths, numY
    global numArcs, arcs, d_x, q, b_x, origin, destination, d_y, b_y, F
    
    y_pos = findall(y_sol.>0)
    numY_pos = length(y_pos)
    TOL = 1e-5
    mod1 = Model(() -> Gurobi.Optimizer(gurobi_env))
    set_optimizer_attribute(mod1, "OutputFlag", 0)
    @variable(mod1, f[1:destination, 1:numY_pos])
    @variable(mod1, v[1:numArcs], Bin)
    A_minus_t = findall(arcs[:, 2] .!= destination)
    for a = 1:numArcs
        i = arcs[a,1]
        j = arcs[a,2]
        for _index = 1:numY_pos
            m = y_pos[_index]
            M = M_set[m]
            if a in M
                @constraint(mod1, f[i,_index] <= (1-q[a])*f[j,_index] + (1-v[a]))
            else
                @constraint(mod1, f[i,_index] <= (1-p[a])*f[j,_index] + (1-v[a]))
            end
        end
    end
    for m = 1:numY_pos
        @constraint(mod1, f[destination,m] == 1)
    end
#     println("1")
    #Flow constraint:
    outArcs = findall(arcs[:, 1].==origin)
    @constraint(mod1, sum(v[a] for a in outArcs) == 1)
    for i = (origin+1):(destination - 1)
        outArcs = findall(arcs[:, 1].==i)
        inArcs = findall(arcs[:, 2].==i)
        @constraint(mod1, sum(v[a] for a in outArcs) - sum(v[a] for a in inArcs) == 0)
    end
    #Constraint on interdicted arcs
    intArcs = findall(x_sol.>=1-TOL)
#     println("Inside Row Gen inArcs = ", intArcs)
    for a in intArcs
       @constraint(mod1, v[a] == 0) 
    end
    #Obj
    @objective(mod1, Max, sum(f[origin, m]*y_sol[y_pos[m]] for m = 1:numY_pos))
    optimize!(mod1)
#     println("2")
    if termination_status(mod1) == MOI.OPTIMAL
        f_s = JuMP.objective_value.(mod1)
    #     println(f_s)
        arcsinPath =findall(JuMP.value.(v) .>=1-TOL)
#         println("v = ", JuMP.value.(v))
#         println("preProcess arcsinPath = ", arcsinPath)
#         println("preProcess arcsinPath = ", findall(JuMP.value.(v) .> 1 - TOL))
        f_val = JuMP.value.(f)
        #If exists a separate cycle - unique won't work
        #May be able to put this if statement back if p_a is used and p_a < 1
#         _nodes = arcs[arcsinPath,1]
#         if length(_nodes) < length(unique(_nodes))
        
        c_g = ones(length(arcsinPath))
        myPath, gx = shortestPath_BellmanFord(c_g, arcs[arcsinPath,:], length(arcsinPath), origin, destination)
        arcsinPath = arcsinPath[findall(myPath.>=1-TOL)]
#         println("after ", arcsinPath)
#         end
#         if arcsinPath == [1, 5, 10, 12, 29]
#             y_pos = findall(y_sol.>0)
            
# #             println(mod1)
#             println("intArcs =",intArcs)
#             println("y_pos = ", y_sol[y_pos])
#             println("M_pos = ", M_set[y_pos])
#             println("f_s = ", f_s)
#             println("f_val = ", f_val)
#         end
    else
        f_s = -1
        arcsinPath = []
        f_val = []
    end
    
#     println(arcsinPath)
#     println(JuMP.value.(f))
    return f_s, arcsinPath, f_val
end
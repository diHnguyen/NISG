function IP_ColGen(lambda)
    global P_set, M_set, numPaths, numY
    global numArcs, arcs, d_x, q, p, b_x, origin, destination, d_y, b_y, Q
    TOL = 1e-3
#     println("lambda ", lambda)
    lambda_pos = findall(lambda.>0)
#     println("lambda_pos = ", lambda_pos)
    numPaths_pos = length(lambda_pos)
    mod2 = Model(() -> Gurobi.Optimizer(gurobi_env))
    set_optimizer_attribute(mod2, "OutputFlag", 0)
    @variable(mod2, f[1:destination, 1:numPaths_pos])
    @variable(mod2, w[1:numArcs], Bin)
    @constraint(mod2, sum(d_y[a]*w[a] for a = 1:numArcs) <= b_y)
    A_minus_t = findall(arcs[:, 2] .!= destination)
# #     println("A")
# #     println("numPaths_pos = ", numPaths_pos)
    for _index = 1:numPaths_pos
        k = lambda_pos[_index]
#         println("p = ", p)
        arcsinPath = P_set[k]
#         println("P = ", arcsinPath)
        for a in arcsinPath
            i = arcs[a,1]
            j = arcs[a,2]
            @constraint(mod2, f[i,_index] >= (1-p[a]-w[a]*(q[a]-p[a]))*f[j,_index])
        end
    end
#     println("Here1")
    for _index = 1:numPaths_pos
        @constraint(mod2, f[destination,_index] == 1)
    end
#     println("Here")
    #Obj
    @objective(mod2, Min, sum(f[origin, _index]*lambda[lambda_pos[_index]] for _index = 1:numPaths_pos))
    optimize!(mod2)
    
    if termination_status(mod2) == MOI.OPTIMAL
        f_s = JuMP.objective_value.(mod2)
    #     println(f_s)
        f_val = JuMP.value.(f)
#         w_val = Int64[]
        w_val = JuMP.value.(w)
        #Precision problem
        M = findall(w_val .>=1-TOL)
#         println("\t\tInside == : ", M )
#         println("\t\tInside > : ", findall(JuMP.value.(w) .>0) )
#         println("\t\t w[13] == 1?", (w_val[13]==1), "\t", w_val[13])
    else
        f_s = 2
        M = []
        f_val = 0
    end

#     if numPaths_pos == 3
#         println(JuMP.value.(f))
#         println(mod2)
#     end
#     println(M)
#     println(JuMP.value.(f))
    return f_s, M, f_val
end
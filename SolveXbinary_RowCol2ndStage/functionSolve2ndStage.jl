function solve2ndStage(x_fix, newElements)
    global numArcs, numPaths, numY, Q, R, P, M, P_Set, M_Set
#     global constr, conNum
    cRefNum = 200
    constr = Array{JuMP.ConstraintRef}(undef, cRefNum)
    m = Model(() -> Gurobi.Optimizer(gurobi_env))
#     @variable(m, 1>=x[1:numArcs]>=0 )
    @variable(m, y[1:numY] >= 0)
    @variable(m, u)
    @objective(m, Min, u)
    con0 = @constraint(m, sum(y[m] for m = 1:numY) == 1)
    
    x_now = zeros(numArcs)
    for a in x_fix
       x_now[a] = 1 
    end
    for i = 1:numPaths
        #V1
        constr[i] = @constraint(m, u >= -sum(sum(log(1-q[a]) for a in intersect(P_set[i],M_set[m]))*y[m] for m = 1:numY) - 10^6*(sum(x_now[a] for a in P_set[i])) )

        #V2
#         constr[i] = @constraint(m, u <= -sum(sum(log(1-q[a]) for a in intersect(P_set[i],M_set[m]))*y[m] for m = 1:numY) + 10^6*(sum(x_now[a] for a in P_set[i])) ) 

    end
    
    iter = 0
    newPathAdded = true
#     lastPath
#     P, M
    newP = true
    newM = true

    while newP == true || newM == true # newPathAdded == true
        newP = false
        newM = false
        iter = iter + 1
#         println("\nIter ", iter)
#         println("M_set = ", M_set)
#         println("P_set = ", P_set)
#         println(m)
        set_optimizer(m, ()-> Gurobi.Optimizer(gurobi_env))
        optimize!(m)
#         println("1")
    #     last_Obj = current_Obj
        current_Obj = JuMP.objective_value.(m)
        y_now = JuMP.value.(y)
#         x_now = JuMP.value.(x)
#         println("x = ", findall(x_now.>0))
#         println("y = ", findall(y_now.>0))
        lambda = zeros(numPaths)
#         println("1")
        for i = 1:numPaths
            lambda[i] = JuMP.dual(constr[i])
#             lambda[i] = JuMP.dual(myCons[i])
        end
        path, gx = shortestPath_BellmanFord(x_now, y_now, arcs)
        arcsinPath = findall(path.>0)
        
        p_exact = findPathsExactPr(arcsinPath, M_set, numY, y_now)
#         println("new path ", arcsinPath)
        if (arcsinPath in P_set) == false #&& p_exact > current_Obj
            numPaths = numPaths + 1
            push!(P_set, arcsinPath)
            println("Add path ", arcsinPath)
            newElements = true
            newP = true
#             println("path = ", arcsinPath, "; gx = ", gx)
#             println("Add Path ", arcsinPath, ": ", p)
#             println("Add new path", arcsinPath)
            coef_q, Q = findQgivenP(M_set, numY, Q, arcsinPath)
            #V1
            constr[numPaths] = @constraint(m, u >= -sum(sum(log(1-q[a]) for a in intersect(P_set[numPaths],M_set[m]))*y[m] for m = 1:numY) - 10^6*(sum(x_now[a] for a in P_set[numPaths])) )
            #V2
#             constr[numPaths] = @constraint(m, u <= -sum(sum(log(1-q[a]) for a in intersect(P_set[i],M_set[m]))*y[m] for m = 1:numY) + 10^6*(sum(x_now[a] for a in P_set[numPaths])) )
        else
#             P_term = true
            
            M = genMonitoring(lambda, P_set)
#             println("new M = ", M)
            
            #Setting names and bounds on anonymous variables must follow the right naming convention
            #https://jump.dev/JuMP.jl/stable/manual/variables/#Anonymous-JuMP-variables
            if M in M_set || isempty(M) == true
                M_term = true
            else
                println("Add M ", M)
                newM = true
                newElements = true
#                 println("Add new monitoring decision ", M)
                numY += 1
                push!(M_set, M)
                new_var = @variable(m, [numY], base_name = "y", lower_bound = 0)
#                 new_var = @variable(MP, [numY], base_name = "y", lower_bound = 0)
                push!(y, new_var[numY])
#                 M_q, 
                Q = findQgivenM(P_set, numPaths, Q, M)
#                 println("M_q = ", M_q)
                for i = 1:numPaths
                    path = P_set[i]
                    set_normalized_coefficient(constr[i], new_var[numY], -Q[i][numY])
                    set_normalized_coefficient(con0, new_var[numY], 1)
                end
            end
        end
    end
    
#     Q = []
    
    optimize!(m)

    current_Obj = JuMP.objective_value.(m)
    y_now = JuMP.value.(y)
#     x_now = JuMP.value.(x)
#     println(m)
#     println("Obj = ", current_Obj)
#     println("x_now = ", findall(x_now.>0), "; ", x_now[x_now.>0])
#     println("M_set ", M_set)
#     println("y_now = ", y_now)
#     println("P_set ", P_set)

    return M_set, P_set, newElements
end
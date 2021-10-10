function solve2ndStage(x_fix)
    global numArcs, numPaths, numY, R, P_Set, M_Set, b_y, d_y,q
    cRefNum = 200
    constr = Array{JuMP.ConstraintRef}(undef, cRefNum)
    m = Model(() -> Gurobi.Optimizer(gurobi_env))
    @variable(m, y[1:numY] >= 0)
    @variable(m, u)
    @objective(m, Min, u)
    
    con0 = @constraint(m, sum(y[m] for m = 1:numY) == 1)
    
    x_now = zeros(numArcs)
    for a in x_fix
       x_now[a] = 1 
    end
    for i = 1:numPaths
        constr[i] = @constraint(m, u >= sum(sum(log(1-q[a]) for a in intersect(P_set[i],M_set[m]))*y[m] for m = 1:numY) - R*(sum(x_now[a] for a in P_set[i])) )
    end
    
    iter = 0
    newP = true
    newM = true

    while newP == true || newM == true 
        newP = false
        newM = false
        iter = iter + 1
        println("\nIter ", iter)
        println("M_set = ", M_set)
        println("P_set = ", P_set)
#         println(m)
        set_optimizer(m, ()-> Gurobi.Optimizer(gurobi_env))
        optimize!(m)

        current_Obj = JuMP.objective_value.(m)
        y_now = JuMP.value.(y)
        println("Obj = ", current_Obj)#, "; y_now = ", y_now)
        println("y = ", findall(y_now.>0))
        lambda = zeros(numPaths)
        for i = 1:numPaths
            lambda[i] = JuMP.dual(constr[i])
        end
        pi_ = JuMP.dual(con0)
        println("pi = ", pi_)
        path, gx = shortestPath_BellmanFord(x_now, y_now, arcs, q, numArcs, M_set, numY)
        
        arcsinPath = findall(path.>0)
        println("P = ", arcsinPath, " gx = ", gx)
        if gx > current_Obj && (arcsinPath in P_set) == false 
            numPaths = numPaths + 1
            push!(P_set, arcsinPath)
#             println("Add path ", arcsinPath, " gx = ", gx, " u_hat = ", current_Obj)
            newP = true
            
            constr[numPaths] = @constraint(m, u >= sum(sum(log(1-q[a]) for a in intersect(P_set[numPaths], M_set[m]))*y[m] for m = 1:numY) - R*(sum(x_now[a] for a in P_set[numPaths])) )
        else
#             println("lambda = ", lambda)
            M, colGen_Obj = genMonitoring(lambda, pi_, P_set, M_set, numY,numPaths)
            println("M = ", M, "; ", colGen_Obj)
            #Setting names and bounds on anonymous variables must follow the right naming convention
            #https://jump.dev/JuMP.jl/stable/manual/variables/#Anonymous-JuMP-variables
            if isempty(M) == false #&& (M in M_set) == false
#                 println("Add M ", M, "; colGen_Obj = ", colGen_Obj)
                newM = true
                numY += 1
                push!(M_set, M)
                new_var = @variable(m, [numY], base_name = "y", lower_bound = 0)
                push!(y, new_var[numY])
                for i = 1:numPaths
                    path = P_set[i]
#                     println("P = ", path, "; M = ", M_set[numY])
                    McapP = intersect(P_set[i], M_set[numY])
                    if isempty(McapP) == false
                        set_normalized_coefficient(constr[i], new_var[numY], -sum(log(1-q[a]) for a in McapP) )
                    else
                        set_normalized_coefficient(constr[i], new_var[numY], 0)
                    end
                    set_normalized_coefficient(con0, new_var[numY], 1)
                end
            end
        end
    end
    
    optimize!(m)
#     println(m)
    current_Obj = JuMP.objective_value.(m)
    y_now = JuMP.value.(y)
#     x_now = JuMP.value.(x)
#     println(m)
#     println("Obj = ", current_Obj)
#     println("x_now = ", findall(x_now.>0), "; ", x_now[x_now.>0])
#     println("M_set ", M_set)
#     println("y_now = ", y_now)
#     println("P_set ", P_set)

    return M_set, P_set, current_Obj, y_now, numY, numPaths
end
function solveLP_findM(x_now)
    global numArcs, arcs, d_x, q, b_x, origin, destination, R, d_y, b_y
    global newP, newM, P_set, M_set, numPaths,numY
    #When P_set gets large enough, add ample space for constraint
    cRefNum = max(200,2*length(P_set))
    constr = Array{JuMP.ConstraintRef}(undef, cRefNum)
    lambda = []
    pi_ = 0
    y_now = []
    
    
    m = Model(() -> Gurobi.Optimizer())
    set_optimizer_attribute(m, "OutputFlag", 0)
    @variable(m, y[1:numY] >= 0)
    @variable(m, u)
    @objective(m, Min, u)
    
    con0 = @constraint(m, sum(y[m] for m = 1:numY) == 1)
#     println("q[1] = ", q[1], " --> ", log10(1-q[1]))
    for i = 1:numPaths
        
        constr[i] = @constraint(m, u >= sum(
                (1+sum(log10(1-q[a]) for a in intersect(P_set[i],M_set[m])))*y[m] for m = 1:numY) - R*(sum(x_now[a] for a in P_set[i])) )
    end
#     println(m)
    M = [-1]
    #Initialize newM and newP
    newM = true
    newP = true
#     println("\n x = ",findall(x_now.>0))
#     println("\t2nd Stage:")
    current_Obj = 0
    it = 0
    while newP == true || newM == true
        newP = false
        newM = false
        optimize!(m)
        current_Obj = JuMP.objective_value.(m)
        y_now=JuMP.value.(y)
    
#         println("\tcurrent_Obj = ",current_Obj)
#         println("\tP_set = ", P_set)
#         println("\tnumPaths = ", numPaths)
        lambda = zeros(numPaths)
        for i = 1:numPaths
            lambda[i] = JuMP.dual(constr[i])
        end
        
        pi_ = JuMP.dual(con0)
#         println("\tlambda = ", lambda)
#         println("\tpi = ", pi_)
#         println(m)
        it = it + 1
#         println(m)
#         println("\t",it, ". Obj_LP = ", current_Obj, "\t", y_now)
        c_g = getCost_SP(x_now, y_now, arcs, q, numArcs, M_set, numY, R)
        path, gx = shortestPath_BellmanFord(c_g, arcs, numArcs, origin, destination)

        arcsinPath = findall(path.>0)
#         println("\t", arcsinPath, " : " ,gx)

        if current_Obj < gx+1 && (arcsinPath in P_set) == false 
            newP = true
#             println("newP = ", newP, " add P", arcsinPath)
            numPaths = numPaths + 1
            push!(P_set, arcsinPath)
            
#             intersect(P_set[numPaths], M_set[m])
            
            constr[numPaths] = @constraint(m, u >= sum((1+sum(log10(1-q[a]) for a in intersect(P_set[numPaths], M_set[m])))*y[m] for m = 1:numY) - R*(sum(x_now[a] for a in P_set[numPaths])) )
            
        end 
            #while isempty(M) == false 
        if newP == false
#             println("\t\tnewP = ", newP, " find M")
            M = []
#             set_optimizer(m, ()-> Gurobi.Optimizer(gurobi_env))
#             optimize!(m)

#             current_Obj = JuMP.objective_value.(m)
#             y_now = JuMP.value.(y)

#             println("pi = ", pi_)
#             println("numPaths = ", numPaths)
            M, colGen_Obj = genMonitoring(lambda, pi_, P_set,numPaths)
#             println("\tM = ", M, "; ", colGen_Obj)
            #Setting names and bounds on anonymous variables must follow the right naming convention
            #https://jump.dev/JuMP.jl/stable/manual/variables/#Anonymous-JuMP-variables
            if isempty(M) == false #&& (M in M_set) == false
                newM = true
                numY += 1
#                 println("M = ", M)
#                 println("M_set = ", M_set)
                push!(M_set, M)
#                 println("M_set = ", M_set)
                new_var = @variable(m, [numY], base_name = "y", lower_bound = 0)
                push!(y, new_var[numY])
#                 println("HERE")
                for i = 1:numPaths
                    path = P_set[i]
#                         println("\tP = ", path, "; M = ", M_set[numY])
                    McapP = intersect(P_set[i], M_set[numY])
#                     println("McapP = ", McapP, log10(1-q[3]))
                    if isempty(McapP) == false
                        set_normalized_coefficient(constr[i], new_var[numY], 
                            -(1+sum(log10(1-q[a]) for a in McapP) ) )
#                         println("new coef = ",-(1+sum(log10(1-q[a]) for a in McapP))  )
                    else
                        set_normalized_coefficient(constr[i], new_var[numY], -1)
                    end
                    set_normalized_coefficient(con0, new_var[numY], 1)
                end
            end
        end
    end
#     x_now = JuMP.value.(x)
#     println(m)
#     println("Obj = ", current_Obj)
#     println("x_now = ", findall(x_now.>0), "; ", x_now[x_now.>0])
#     println("M_set ", M_set)
#     println("y_now = ", y_now)
#     println("P_set ", P_set)
#     println("\tlambda = ", lambda)
#     println("\tpi = ", pi_)
    
    
    return P_set, M_set, lambda, pi_, current_Obj, y_now
end
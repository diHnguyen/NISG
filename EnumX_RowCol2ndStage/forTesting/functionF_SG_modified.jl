function F_SG(x_now)
    global numArcs, arcs, d_x, q, b_x, origin, destination, R, d_y, b_y, F
    global newP, newM, P_set, M_set, numPaths,numY
    global t_exactRow, t_exactCol
#     global gurobi_env
    #When P_set gets large enough, add ample space for constraint
#     println("A")
    cRefNum = numPaths #length(P_set) #max(700,3*length(P_set))
    constr = Array{JuMP.ConstraintRef}(undef, cRefNum)
    lambda = []
    pi_ = 0
    y_now = []
    
#     m = direct_model(Gurobi.Optimizer(gurobi_env))
    m = Model(() -> Gurobi.Optimizer(gurobi_env))
    # set_optimizer_attribute(m,"Threads", 1)
    set_optimizer_attribute(m, "OutputFlag", 0)
    @variable(m, y[1:numY] >= 0)
    @variable(m, u >= 0)
    @objective(m, Min, u)
#     println("B")
    con0 = @constraint(m, sum(y[m] for m = 1:numY) == 1)
    for i = 1:numPaths
        constr[i] = @constraint(m, u >= sum(F[i][m]*y[m] for m = 1:numY) - (sum(x_now[a] for a in P_set[i])) )
    end

    M = [-1]
    #Initialize newM and newP
    newM = true
    newP = true
    continueApprox = true
    TOL = 1e-5
#     println("\n x = ",findall(x_now.>0))
#     println("\t2nd Stage:")
    current_Obj = 0
    it = 0
    
#     println("C")
    while newP == true || newM == true
        newP = false
        newM = false
        @timeit to "F_SG" optimize!(m)
        
        current_Obj = JuMP.objective_value.(m)
        y_now=JuMP.value.(y)
        println("\t",it)
#         println(m)
        println("\tcurrent_Obj = ",current_Obj)
        y_pos = findall(y_now.>0)
         println("\ty_pos = ", y_now[y_pos])
         println("\tM_pos = ", M_set[y_pos])
        println("\tP_set = ", P_set)
#         println("\tnumPaths = ", numPaths, " vs ", length(P_set))
        
        lambda = zeros(numPaths)
        for i = 1:numPaths
            lambda[i] = JuMP.dual(constr[i])
        end
        pi_ = JuMP.dual(con0)
        it = it + 1
#         println("\t",it, ". Obj_LP = ", current_Obj)#, "\t", y_now)
#         println("y_pos = ", y_now[y_pos])
#         println("\tM_pos = ", M_set[y_pos])
#         println("\tlambda = ", lambda)
#         println("\tpi = ", pi_)
        
#         if numPaths == 28
        lambda_pos = findall(lambda.>0)
#              println("\tlambda_pos = ", lambda[lambda_pos])
#              println("\tP_pos = ", P_set[lambda_pos])
        
#             for i in lambda_pos
#                 println("\tlambda ", lambda[i], "\t", P_set[i])
#     #         println("")
#             end
#         end
        
#         println("\tP_set = ", P_set)
#         println("PreApprox")
        if continueApprox == true
             println("\tApprox...")
            println("lambda = ", lambda)
            println("lambda_pos = ", lambda_pos)
            println("pi = ", pi_)
            c_g = getCost_SP(x_now, y_now)
#             println("\tSolving for Row")
            path, gx = shortestPath_BellmanFord(c_g, arcs, numArcs, origin, destination)
#             println("\tDone - Solving for Row")
            new_F_PM = 0
            arcsinPath = findall(path.>TOL)
            if sum(x_now[a] for a in arcsinPath) < TOL
#                 println("\tPath is not interdicted")
                for i in y_pos
                    M = M_set[i]
                    new_F_PM = new_F_PM + findf_MP(arcsinPath, M)*y_now[i]
                end
                println("\tPath ",arcsinPath, "\t", new_F_PM)
                if current_Obj + TOL < new_F_PM #sum(F[numPaths][m]*y_now[m] for m=1:numY)  
#                     println("\tPath ",arcsinPath, "\t", new_F_PM)
                    newP = true
                    numPaths = numPaths + 1
                    push!(P_set, arcsinPath)
                    F = updateF(F, P_set, M_set, numY, numPaths)
                    con = @constraint(m, u >= sum(F[numPaths][m]*y[m] for m = 1:numY) - (sum(x_now[a] for a in arcsinPath)) )
                    push!(constr, con)
                end
            end
#             println("doneApproxPath")
            if newP == false
#                 println("\tnewP = ", newP, " find M")
                M = []
                M, colGen_Obj = LP_ColGen(lambda)
                if isempty(M) == false
#                     println("\tM isn't empty")
                    new_F_PM = 0
                    for i in lambda_pos
                        P = P_set[i]
#                         println("\tF_val =",findf_MP(P, M))
                        new_F_PM = new_F_PM + findf_MP(P, M)*lambda[i]
                        
                    end
                    
                      
#                     if M == [3, 4]
#                         println("lambda_pos = ", lambda_pos)
#                         println("lambd = ", lambda[lambda_pos])
#                         println("P_pos = ", P_set[lambda_pos])
#                         println("pi_ = ", pi_)
#                     end
                    println("\tMon ",M, "\t", new_F_PM - pi_)
                    if new_F_PM - pi_< -TOL 
#                         println("\tMon ",M, "\t", new_F_PM - pi_)
#                         println("\tP_set length = ", length(P_set))
#                         println("\tP_Set =", P_set)
                        newM = true
                        numY += 1
                        push!(M_set, M)
                        new_var = @variable(m, [numY], base_name = "y", lower_bound = 0)
                        push!(y, new_var[numY])
#                         println("Pre-update F")
                        F = updateF(F, P_set, M_set, numY, numPaths)
#                         println("Post-update F")
#                         println("numPaths = ", numPaths)
#                         println("F = ", F)
#                         println("numY = ", numY)
#                         println("M_set = ", M_set)
                       
                        for i = 1:numPaths
#                             println("\tpath ", i)
#                             println(constr[i])
                            set_normalized_coefficient(constr[i], new_var[numY], -F[i][numY]) 
                            set_normalized_coefficient(con0, new_var[numY], 1)
                        end
#                         println("Post update coefs")
                    end
#                     println("\tDone setting new coef for M")
                end
#                 println("\tDone solving for M")
            end
#             println("doneApproxMon")
            if newP == false && newM == false
                continueApprox = false
            end
#             println("continueApprox = ", continueApprox)
        end
#         println("preExactPath")
        if  continueApprox == false
            println("\tExact...")
            println("lambda = ", lambda)
            println("lambda_pos = ", lambda_pos)
            println("pi = ", pi_)
#             t1 = start()
#             t2 = start()-t1
#             println("\t t = ", t2)
            @timeit to "IP_RowGen" f_s, arcsinPath, f_val = IP_RowGen(x_now, y_now)
#             println("Here")
            
#             t_exactRow = t_exactRow + start()-t1
            println("\tRow Gen: ", f_s, " \t", arcsinPath)
            if current_Obj + TOL < f_s
    #             println("\tRow Gen: ", f_s, " \t", arcsinPath)
    #             println("\t arcs ", arcs[arcsinPath,:], "\t\t??? ", (arcsinPath in P_set))
    #             if arcsinPath in P_set
    #                 println("\n")
    #                 println("\t\t", "P_set = ", P_set)
    #                 println("\t\t", "y_now = ", y_now)
    #                 println("\t\t", "M_set = ", M_set)
    #             end
    #         if current_Obj < gx+1 && (arcsinPath in P_set) == false 
                newP = true
    #             println("newP = ", newP, " add P", arcsinPath)
                numPaths = numPaths + 1
                push!(P_set, arcsinPath)

                F = updateF(F, P_set, M_set, numY, numPaths)

                con = @constraint(m, u >= sum(F[numPaths][m]*y[m] for m = 1:numY) - (sum(x_now[a] for a in arcsinPath)) )
                push!(constr, con)
    #             constr[numPaths] = @constraint(m, u >= sum((1+sum(log10(1-q[a]) for a in intersect(P_set[numPaths], M_set[m])))*y[m] for m = 1:numY) - R*(sum(x_now[a] for a in P_set[numPaths])) )
            end 
#             println("doneExactPath")
            if newP == false
#                 println("\tnewP = ", newP, " find M")
                M = []
                M, colGen_Obj = LP_ColGen(lambda)
                if isempty(M) == false
#                     println("\tM isn't empty")
                    new_F_PM = 0
                    for i in lambda_pos
                        P = P_set[i]
#                         println("\tF_val =",findf_MP(P, M))
                        new_F_PM = new_F_PM + findf_MP(P, M)*lambda[i]
                        
                    end
#                     if M == [3, 4]
#                         println("lambda_pos = ", lambda_pos)
#                         println("lambd = ", lambda[lambda_pos])
#                         println("P_pos = ", P_set[lambda_pos])
#                         println("pi_ = ", pi_)
#                     end
                    println("\tMon ",M, "\t", new_F_PM - pi_)
                    if new_F_PM - pi_< -TOL 
#                         println("\tMon ",M, "\t", new_F_PM - pi_)
#                         println("\tP_set length = ", length(P_set))
#                         println("\tP_Set =", P_set)
                        newM = true
                        numY += 1
                        push!(M_set, M)
                        new_var = @variable(m, [numY], base_name = "y", lower_bound = 0)
                        push!(y, new_var[numY])
#                         println("Pre-update F")
                        F = updateF(F, P_set, M_set, numY, numPaths)
#                         println("Post-update F")
#                         println("numPaths = ", numPaths)
#                         println("F = ", F)
#                         println("numY = ", numY)
#                         println("M_set = ", M_set)
                       
                        for i = 1:numPaths
#                             println("\tpath ", i)
#                             println(constr[i])
                            set_normalized_coefficient(constr[i], new_var[numY], -F[i][numY]) 
                            set_normalized_coefficient(con0, new_var[numY], 1)
                        end
#                         println("Post update coefs")
                    end
#                     println("\tDone setting new coef for M")
                end
#                 println("\tDone solving for M")
            end
#             println("doneExactMon")
        end
#         println("\tf_val = ", f_val)
#         c_g = getCost_SP(x_now, y_now, arcs, q, numArcs, M_set, numY, R)
#         path, gx = shortestPath_BellmanFord(c_g, arcs, numArcs, origin, destination)

#         arcsinPath = findall(path.>0)
#         println("\t", arcsinPath, " : " ,gx)
        
        
            #while isempty(M) == false 
        
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
#      println("\tt_F_SG = ", TimerOutputs.time(to["F_SG"])/10^9)
#      println("\tt_RowGen = ", TimerOutputs.time(to["IP_RowGen"])/10^9)
#      println("\tt_ColGen = ", TimerOutputs.time(to["IP_ColGen"])/10^9)
    
    return lambda, pi_, current_Obj, y_now
end
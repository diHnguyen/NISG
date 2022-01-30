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
        optimize!(m)
        
        current_Obj = JuMP.objective_value.(m)
        y_now=JuMP.value.(y)
#         println("\t",it)
#         println(m)
#         println("\tcurrent_Obj = ",current_Obj)
        y_pos = findall(y_now.>0)
#         println("\ty_pos = ", y_now[y_pos])
#         println("\tM_pos = ", M_set[y_pos])
#         println("\tP_set = ", P_set)
#         println("\tnumPaths = ", numPaths, " vs ", length(P_set))
        
        lambda = zeros(numPaths)
        for i = 1:numPaths
            lambda[i] = JuMP.dual(constr[i])
        end
        pi_ = JuMP.dual(con0)
        it = it + 1
        println("\t",it, ". Obj_LP = ", current_Obj)#, "\t", y_now)
#         println("y_pos = ", y_now[y_pos])
#         println("M_pos = ", M_set[y_pos])
#         println("\tlambda = ", lambda)
#         println("\tpi = ", pi_)
        
#         if numPaths == 28
        lambda_pos = findall(lambda.>0)
#             println("lambda_pos = ", lambda[lambda_pos])
#             println("P_pos = ", P_set[lambda_pos])
        
#             for i in lambda_pos
#                 println("\tlambda ", lambda[i], "\t", P_set[i])
#     #         println("")
#             end
#         end
        
#         println("\tP_set = ", P_set)
#         println("PreApprox")
        if continueApprox == true
            println("\tApprox...")
            c_g = getCost_SP(x_now, y_now)
            path, gx = shortestPath_BellmanFord(c_g, arcs, numArcs, origin, destination)
            new_F_PM = 0
            arcsinPath = findall(path.>TOL)
            if sum(x_now[a] for a in arcsinPath) < TOL
                for i in y_pos
                    M = M_set[i]
                    new_F_PM = new_F_PM + findf_MP(arcsinPath, M)*y_now[i]
                end
                println("\tPath ",arcsinPath, "\t", new_F_PM)
                if current_Obj + TOL < new_F_PM #sum(F[numPaths][m]*y_now[m] for m=1:numY)  
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
#             println("\t\tnewP = ", newP, " find M")
                M = []
                M, colGen_Obj = LP_ColGen(lambda)
                if isempty(M) == false
                    new_F_PM = 0
                    for i in lambda_pos
                        P = P_set[i]
                        new_F_PM = new_F_PM + findf_MP(arcsinPath, M)*lambda[i]
                    end
                    println("\tMon ",M, "\t", new_F_PM)
                    if new_F_PM - pi_< -TOL 
                        newM = true
                        numY += 1
                        push!(M_set, M)
                        new_var = @variable(m, [numY], base_name = "y", lower_bound = 0)
                        push!(y, new_var[numY])
                        F = updateF(F, P_set, M_set, numY, numPaths)
                        for i = 1:numPaths
                            set_normalized_coefficient(constr[i], new_var[numY], -F[i][numY]) 
                            set_normalized_coefficient(con0, new_var[numY], 1)
                        end
                    end
                end
            end
#             println("doneApproxMon")
            if newP == false && newM == false
                continueApprox = false
            end
        end
#         println("preExactPath")
        if  continueApprox == false
            println("\tExact...")
#             t1 = start()
#             t2 = start()-t1
#             println("\t t = ", t2)
            @timeit to "IP_RowGen" f_s, arcsinPath, f_val = IP_RowGen(x_now, y_now)
#             println("Here")
            
#             t_exactRow = t_exactRow + start()-t1
            
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
#             println("\t\tnewP = ", newP, " find M")
                M = []
#                 t1 = start()
                @timeit to "IP_ColGen" f_s, M, f_val = IP_ColGen(lambda)
#                 t_exactCol = t_exactCol + start() -t1
    #             tempq = ones(numPaths)
    #             for i = 1:numPaths
    #                 P = P_set[i]
    #                 McapP = intersect(M, P)
    #                 if isempty(McapP) ==false
    #                     tempq[i] = prod((1-q[a]) for a in McapP)
    #                 end
    #             end


    #             if sum(lambda[i]*tempq[i] for i =1:numPaths) - pi_ > f_s - pi_ + TOL
    #                 println("\tf_s ", f_s)
    #                 println("\tf_val ", f_val)
    #                 println("\tP_pos = ", P_set[lambda_pos])
    #             end
                #Setting names and bounds on anonymous variables must follow the right naming convention
                #https://jump.dev/JuMP.jl/stable/manual/variables/#Anonymous-JuMP-variables

                if f_s - pi_< -TOL #0 #isempty(M) == false &&  #&& (M in M_set) == false
                    #             println("M = ", M)
    #                 println("\tTest f_s - pi_ = ", sum(lambda[i]*tempq[i] for i =1:numPaths) - pi_)
        #             M, colGen_Obj = genMonitoring(lambda, pi_, P_set,numPaths)
        #             println()
    #                 println("\tCol Gen ",f_s - pi_, "\tM = ", M) #, "\t costs ", sum(d_y[a] for a in M))
                    newM = true
                    numY += 1
    #                 println("Add M")
    #                 println("M = ", M)
    #                 println("M_set = ", M_set)
                    push!(M_set, M)
    #                 println("M_set = ", M_set)
                    new_var = @variable(m, [numY], base_name = "y", lower_bound = 0)
                    push!(y, new_var[numY])
    #                 println("\tM_set ", M_set)
                    F = updateF(F, P_set, M_set, numY, numPaths)
    #                 println("HERE")
                    for i = 1:numPaths
    #                     path = P_set[i]
    #                         println("\tP = ", path, "; M = ", M_set[numY])
    #                     McapP = intersect(P_set[i], M_set[numY])
    #                     println("McapP = ", McapP, log10(1-q[3]))
    #                     if isempty(McapP) == false
                            set_normalized_coefficient(constr[i], new_var[numY], -F[i][numY]) 
    #                         println("new coef = ",-(1+sum(log10(1-q[a]) for a in McapP))  )
    #                     else
    #                         set_normalized_coefficient(constr[i], new_var[numY], -1)
    #                     end
                        set_normalized_coefficient(con0, new_var[numY], 1)
                    end
                end
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
    
    
    return lambda, pi_, current_Obj, y_now
end
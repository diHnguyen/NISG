#Implement a simple example of column and row generation
function LP_ColGen(lambda)
    global P_set, M_set, numPaths, numY
    global numArcs, arcs, d_x, q, b_x, origin, destination, d_y, b_y, F
#     global gurobi_env # Gurobi.Env()

    TOL = 1e-5
    y_pos = findall(y_sol.>0)
    numY_pos = length(y_pos)
#     println("P_set = ", P_set)
#     println("numArcs = ", numArcs)
#     println("numpaths = ", numPaths)
    m0 = Model(() -> Gurobi.Optimizer(gurobi_env))
    set_optimizer_attribute(m0, "OutputFlag", 0)
#     set_optimizer_attribute(m0, "LogToConsole", 0)
    @variable(m0, w[1:numArcs], Bin)
    @variable(m0, q_bar[1:numArcs])
#     println("1")
#     println("P_set[2] = ", P_set[2])
#     println(lambda[2]*sum(log10(1-q[a]) for a in P_set[2]))
    
    
    #New Objective
    @objective(m0, Min, sum(lambda[i]*sum(q_bar[a] for a in P_set[i]) for i=1:numPaths))
    #Old Objective
    #@objective(m0, Min, sum(lambda[i]*sum(log(1-q[a])*w[a] for a in P_set[i]) for i=1:numPaths))

    @constraint(m0,sum(d_y[a]*w[a] for a = 1:numArcs) <= b_y)
    
    #Constraint w new Obj:
    for a = 1:numArcs
        @constraint(m0, q_bar[a] >= log(1-p[a]) + w[a]*(log(1-q[a]) - log(1-p[a])) )
    end
    
#     println("3")
    optimize!(m0)
#     println(m0)
#     println("genMon")
    colGen_Obj = JuMP.objective_value.(m0)
    if colGen_Obj < TOL
        w_now = JuMP.value.(w)
        M = findall(w_now.>= 1-TOL)
    else
        M = []
    end
#     println("done")
    return M, colGen_Obj
end

# - pi_ + sum(lambda[i] for i = 1:numPaths)
# if termination_status(m) != MOI.OPTIMAL
#     warn("Master not optimal ($ncols patterns so far)")
# end

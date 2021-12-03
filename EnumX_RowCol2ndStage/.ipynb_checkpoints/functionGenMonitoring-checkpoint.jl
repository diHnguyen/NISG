#Implement a simple example of column and row generation
function genMonitoring(lambda, pi_, P_set, numPaths)
#     global gurobi_env # Gurobi.Env()
    global b_y, d_y, numArcs
    tol = -1e-6
#     println("P_set = ", P_set)
#     println("numArcs = ", numArcs)
#     println("numpaths = ", numPaths)
    m0 = Model(() -> Gurobi.Optimizer())
    set_optimizer_attribute(m0, "OutputFlag", 0)
    @variable(m0, w[1:numArcs], Bin)
#     println("1")
#     println("P_set[2] = ", P_set[2])
#     println(lambda[2]*sum(log10(1-q[a]) for a in P_set[2]))
    @objective(m0, Min, sum(lambda[i]*sum(log10(1-q[a])*w[a] for a in P_set[i]) for i=1:numPaths) - pi_ + sum(lambda[i] for i = 1:numPaths))
#     println("2")
    @constraint(m0,sum(d_y[a]*w[a] for a = 1:numArcs) <= b_y)
#     println("3")
    optimize!(m0)
#     println(m0)
#     println("genMon")
    colGen_Obj = JuMP.objective_value.(m0)
    if colGen_Obj < tol
        w_now = JuMP.value.(w)
        M = findall(w_now.>0)
    else
        M = []
    end
    return M, colGen_Obj
end

# if termination_status(m) != MOI.OPTIMAL
#     warn("Master not optimal ($ncols patterns so far)")
# end

#Implement a simple example of column and row generation
function genMonitoring(lambda, pi_, P_set, M_set, numY,numPaths)
    global gurobi_env # Gurobi.Env()
    global b_y, d_y, numArcs
    m0 = Model(() -> Gurobi.Optimizer(gurobi_env))
    @variable(m0, w[1:numArcs], Bin)
    @objective(m0, Min, lambda[i]*sum(sum(log(1-q[a])*w[a] for a in P_set[i]) for i=1:numPaths) - pi_ +sum(lambda[i] for i = 1:numPaths))
    @constraint(m0,sum(d_y[a]*w[a] for a = 1:numArcs) <= b_y)
    optimize!(m0)
    colGen_Obj = JuMP.objective_value.(m0)
    if colGen_Obj < 0
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

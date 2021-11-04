#Implement a simple example of column and row generation
function genMonitoring(lambda, P_set)
#     println("lambda ", lambda)
    
    global gurobi_env # Gurobi.Env()
    global b_y, d_y, numArcs # = [4,3,1,6,6]
#     println("FeasY = ", EnumX(arcs, b_y, numArcs, d_y))
    m0 = Model(() -> Gurobi.Optimizer(gurobi_env))
# set_silent(m)
# @variable(m, x[1:ncols] >= 0)
    @variable(m0, w[1:numArcs], Bin)
    @objective(m0, Min, sum(sum(log(1-q[a])*w[a] for a in P) for P in P_set))
# for a = 1:ncols
    @constraint(m0,sum(d_y[a]*w[a] for a = 1:numArcs) <= b_y)
#     @constraint(m, v >= x[a])
# end
    optimize!(m0)
    colGen_Obj = JuMP.objective_value.(m0)
    if colGen_Obj < 0
        w_now = JuMP.value.(w)
        M = findall(w_now.>0)
    else
        M = []
    end
#     println("M = ", M)
    return M
end

# if termination_status(m) != MOI.OPTIMAL
#     warn("Master not optimal ($ncols patterns so far)")
# end

# while ncols < length(b)
#     global ncols, b
    
# #     println(m)
# end
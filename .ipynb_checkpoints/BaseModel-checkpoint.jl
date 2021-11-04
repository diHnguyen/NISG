# If we want to add # in Gurobi, then we have to turn of Gurobi's own Cuts 
m0 = Model(() -> Gurobi.Optimizer(gurobi_env))
@variable(m0, 1>=x0[1:numArcs]>=0 )
@variable(m0, y0[1:numY] >= 0)
@variable(m0, u0)
@objective(m0, Min, u0)
con0 = @constraint(m0, sum(y0[m] for m = 1:numY) == 1)
@constraint(m0, sum(x0[a] for a=1:numArcs) <= 1)
constr[numPaths] = @constraint(m0, u0 >= sum(y0[m]*coef_q[m] for m = 1:numY) - R*(sum(x0[a] for a in P_set[1])) )
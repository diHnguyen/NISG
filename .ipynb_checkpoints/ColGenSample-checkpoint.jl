#Implement a simple example of column and row generation
using JuMP
using Gurobi

global gurobi_env = Gurobi.Env()
setparams!(gurobi_env, Heuristics=0.0, Cuts = 0, OutputFlag = 0)
# nwidths = length(prices)
# n = length(widths)
global b = [4,3,1,6,6]
global ncols = 2 #length(widths)


m = Model(() -> Gurobi.Optimizer(gurobi_env))
# set_silent(m)
@variable(m, x[1:ncols] >= 0)
@variable(m, v>=0)
@objective(m, Min, v)# sum(x[a] for a = 1:ncols))
for a = 1:ncols
    @constraint(m,x[a] >= b[a])
    @constraint(m, v >= x[a])
end
optimize!(m)
if termination_status(m) != MOI.OPTIMAL
    warn("Master not optimal ($ncols patterns so far)")
end

while ncols < length(b)
    global ncols, b
    ncols += 1
    #Setting names and bounds on anonymous variables must follow the right naming convention
    #https://jump.dev/JuMP.jl/stable/manual/variables/#Anonymous-JuMP-variables
    new_var = @variable(m, [ncols], base_name = "x", lower_bound = 0)
    push!(x, new_var[ncols])
    @constraint(m, x[ncols] >= b[ncols])
    @constraint(m, v >= x[ncols])
    optimize!(m)
    current_Obj = JuMP.objective_value.(m)
    println("\nIter ", ncols, " : ", current_Obj)
#     println(m)
end
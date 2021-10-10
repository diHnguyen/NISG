#Updated on 6/28
#All monitoring decisions are included in the MP
#The coefficients of each arc in the MP is exact
#The path with the least detection probability = most probability of being undetected is approximate
using JuMP
using Gurobi
using LightGraphs

global gurobi_env = Gurobi.Env()
setparams!(gurobi_env, Heuristics=0.0, Cuts = 0, OutputFlag = 0)
# outfile = "C:/Users/din/Documents/GitHub/UncertainTarget/Log_Test_Output.csv"
#     io = open(outfile, "w")
# println(io,"Instance \tPath \tExact \tApprox")
# global p_exact = []
# global p_approx = []
global ins = 1 #:10

include("./Ins"*string(ins)*".jl")
include("./functionEnumFeasX.jl")
include("./functionFindQ.jl")
include("./functionShortestPathBellmanFord.jl")
# include("./findExactUndetectionPr.jl")
include("./functionGenMonitoring.jl")
include("./functionSolve2ndStage.jl")
global numArcs = length(arcs[:,1])
global X_feas = EnumX(arcs, b_x, numArcs, d_x)
global EnumY = EnumX(arcs, b_y, numArcs, d_y)
println("EnumY = ", EnumY)
global M_set = [EnumY[1]] #EnumX(arcs, b_y, numArcs, d_y)
global numY = 1 
global R = 10^6
#Add at least one path to the initial P_set
global P_set = [[1,8]]
global numPaths = 1
# global cRefNum = 200
# global constr = Array{JuMP.ConstraintRef}(undef, cRefNum)

global x_now
global Obj_X = zeros(length(X_feas))
println("X_feas = ", X_feas)
println("Starting M_set ", M_set, ", P_set = ", P_set)
global Y_sol = []

for x_i = 1:1#length(X_feas)
    global numArcs, arcs, d_x, q, b_x, origin, destination, R, d_y, b_y #, EnumY
    x_now = X_feas[x_i]
    println("\n\nx = ", x_now)
    M_set = [EnumY[1]]
    numY = 1 #length(M_set)
    P_set = [[1,8]]
    numPaths = 1
    
    cRefNum = 200
    constr = Array{JuMP.ConstraintRef}(undef, cRefNum)
    
    M_set, P_set, current_Obj, y_now, numY, numPaths = solve2ndStage(x_now)
    println("1. Obj Val for x = ", x_now,", y = ", y_now," is ", current_Obj)
    
    
    P_set = copy(P)
    M_set = copy(EnumY)
    numPaths = length(P_set)
    numY = length(M_set)
    MP = Model(() -> Gurobi.Optimizer(gurobi_env))
    @variable(MP, x[1:numArcs], Bin)
    @variable(MP, y[1:numY] >= 0)
    @variable(MP, u)
    @objective(MP, Min, u)
    @constraint(MP, sum(y[m] for m = 1:numY) == 1)
    @constraint(MP, sum(d_x[a]*x[a] for a=1:numArcs) <= b_x)
    for a in x_now
        @constraint(MP, x[a] == 1)
    end
    for i = 1:numPaths
         @constraint(MP, u >= sum(
                sum(log(1-q[a]) for a in intersect(P_set[i],M_set[m]))*y[m] for m = 1:numY) - 10^6*(sum(x[a] for a in P_set[i])) ) 
    end
    optimize!(MP)
#     println("\nSolve once given P_set and M_set:")
#     println("P = ", P_set)
#     println("M = ", M_set)
#     println(MP)
    println("2. Obj Val for x = ", findall(JuMP.value.(x).>0),", y = ", JuMP.value.(y)," is ", JuMP.objective_value.(MP))
end
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

include("./TestInstances/Ins"*string(ins)*".jl")
include("./functionEnumFeasX.jl")
include("./functionFindQ.jl")
include("./functionShortestPathBellmanFord.jl")
include("./findExactUndetectionPr.jl")
global numArcs = length(arcs[:,1])
# global numArcs
# global p_exact, p_approx
global feasY = EnumX(arcs, b_y, numArcs, d_y)
global numY = length(feasY)
println("M_set = ", feasY)
#Add at least one path to the initial P_set
global P_set = [[1,8]]
global Q = []
global R = 1000
coef_q, Q = findQ(feasY, numY, Q, P_set[1])
println("Path's probability of not being detected for each M = ", coef_q)
# println("Q = ", Q)


m = Model(() -> Gurobi.Optimizer(gurobi_env))
# @variable(m, x[1:numArcs] >= 0)
@variable(m, y[1:numY] >= 0)
@variable(m, u >= 0)
@objective(m, Min, u)# sum(x[a] for a = 1:ncols))
@constraint(m, sum(y[m] for m = 1:numY) == 1)
# coef_q = Q[1]
@constraint(m, u >= sum(y[m]*coef_q[m] for m = 1:numY))# - R*(sum(x[a] for a in P_set[1])) )

for iter = 1:5
    global numArcs, feasY, numY, P_set, Q, R
    
    println("\nIter ", iter)
    optimize!(m)
    current_Obj = JuMP.objective_value.(m)
    y_now = JuMP.value.(y)
#     x_now = JuMP.value.(y)
    println("Obj = ", current_Obj)
    println("y_now = ", y_now)
    #Get Path Bellman
    if iter != 4
    path, gx = shortestPath_BellmanFord(zeros(numArcs), y_now, arcs)
    arcsinPath = findall(path.>0)
    push!(P, arcsinPath)
    println("path = ", arcsinPath, "; gx = ", gx)
    findPathsExactPr(y_now)
    #Calculate Q_P for each M:
    coef_q, Q = findQ(feasY, numY, Q, arcsinPath)
    println("Path's probability of not being detected for each M = ", coef_q)
    #Add Path to MP
    @constraint(m, u >= sum(y[m]*coef_q[m] for m = 1:numY))
    else
        println("Add new path")
        arcsinPath = [5,7]
        coef_q, Q = findQ(feasY, numY, Q, arcsinPath)
        @constraint(m, u >= sum(y[m]*coef_q[m] for m = 1:numY))
#         optimize!(m)
#         current_Obj = JuMP.objective_value.(m)
#         y_now = JuMP.value.(y)
#         #     x_now = JuMP.value.(y)
#         println("Obj = ", current_Obj)
#         println("y_now = ", y_now)
    end
end


# println(m)
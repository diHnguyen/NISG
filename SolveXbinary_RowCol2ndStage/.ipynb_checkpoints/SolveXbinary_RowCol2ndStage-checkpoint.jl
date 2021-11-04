#Updated on 6/28
#All monitoring decisions are included in the MP
#The coefficients of each arc in the MP is exact
#The path with the least detection probability = most probability of being undetected is approximate


#Notes: Initialization of P_set and M_Set is important to get to the optimal solution.
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
include("./findExactUndetectionPr.jl")
include("./functionGenMonitoring.jl")

include("./functionSolve2ndStage.jl")

global numArcs = length(arcs[:,1])
global X_feas = EnumX(arcs, b_x, numArcs, d_x)

global EnumY = EnumX(arcs, b_y, numArcs, d_y)
println("EnumY = ", EnumY)
global M_set = [EnumY[1]] #EnumX(arcs, b_y, numArcs, d_y)
global numY = 1 #length(M_set)
# println("M_set = ", M_set)
#Add at least one path to the initial P_set
global P_set = [[1,8]]
# println("P_set =", P_set )
global numPaths = 1
global Q = []
global R = 1 #000
coef_q, Q = findQgivenP(M_set, numY, Q, P_set[1])

# coef_q, Q = findQgivenP(M_set, numY, Q, P_set[1])

# for path in P_set
#     global Q
#     coef_q, Q = findQgivenP(M_set, numY, Q, path)
# end

global last_x
global x_now
global terminate = false
global newElements = true
while newElements == true #&& numNodes <7
    global objVal, x_now, y_now, last_x, newElements
    global numArcs, M_set, numY, P_set, numPaths, Q, R 
    global nodes, nodeStat, nodePred, nodeObj, nodeX, numNodes, openNodes
    
    println("\n\nP_set = ",P_set)
    println("M_set = ", M_set)
    println("Q = ", Q)
#     cRefNum = 200
#     constr = Array{JuMP.ConstraintRef}(undef, cRefNum)
#     global newElements
#     for I = 1:5
#     while newElements == true#Mterm == false && Pterm == false 
        MP = Model(() -> Gurobi.Optimizer(gurobi_env))
        @variable(MP, x[1:numArcs], Bin)
        @variable(MP, y[1:numY] >= 0)
        @variable(MP, u)
        @objective(MP, Min, u)
        @constraint(MP, sum(y[m] for m = 1:numY) == 1)
        @constraint(MP, sum(d_x[a]*x[a] for a=1:numArcs) <= b_x)
        for numPaths = 1:length(P_set)
            @constraint(MP, u >= sum(y[m]*coef_q[m] for m = 1:numY) - R*(sum(x[a] for a in P_set[numPaths])) )
        end
        optimize!(MP)
        objVal = JuMP.objective_value.(MP)
        x_now = JuMP.value.(x)
        y_now = JuMP.value.(y)
#         last_x = [-1]
        println("MP Obj = ", objVal, ", x_now = ", findall(x_now.>0), ", y_now = ", y_now)
        newElements = false
        M_set, P_set, newElements = solve2ndStage(x_now, newElements)
#         println("M_set ", M_set)
#         println("P_set ", P_set)
        println("newElements = ", newElements)
#         last_x = x_now
        #Update x_last only if
#         if Mterm == false && Pterm == false 
#             if x_now == x_last
#               terminate = true
#             end
#         end
#         optimize!(MP)
#         x_now = JuMP.value.(x)
#         y_now = JuMP.value.(y)
#     end
#     x_now = x_last
end
println("\n--------------------------")
println("Final obj = ", objVal)
println("Final x = ", findall(x_now.>0))
println("M_set = ", M_set)
println("Final y = ", y_now)
println("P_set = ", P_set)
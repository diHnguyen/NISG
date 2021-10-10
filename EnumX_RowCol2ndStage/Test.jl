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
# global numArcs
# global p_exact, p_approx
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
# println("Path's probability of not being detected for each M = ", coef_q)
# println("Q = ", Q)

# global ncols = 1
global cRefNum = 200
global constr = Array{JuMP.ConstraintRef}(undef, cRefNum)

# global current_Obj = 1
# global last_Obj = 2
# global last_x = ones(numArcs)
# global x_now = zeros(numArcs)
global x_now
global iter = 0
global newPathAdded = true
global lastPath
global P, M
global P_term = false
global M_term = false
# global numPaths
# for iter = 1:5
# while  last_x != x_now || newPathAdded == true
global Obj_X = zeros(length(X_feas))
println("X_feas = ", X_feas)
global Y_sol = []
# global Y
# println("Obj_X = ")

for x_i = 1:length(X_feas) #x_now in X_feas
    x_now = X_feas[x_i]
    println("\n\nx = ", x_now)
    global numArcs, M_set, numY, P_set, Q, R, numPaths, M, P, M_term, P_term, Obj_X, Y_sol
    global iter, lastPath, newPathAdded, constr, cRefNum
    
    newElements = true
    M_set, P_set, current_Obj, y_now = solve2ndStage(x_now)
    
#     current_Obj = JuMP.objective_value.(m)
#     y_now = JuMP.value.(y)
    push!(Y_sol, y_now)
    Obj_X[x_i] = current_Obj
    println("Obj Val for x = ", x_now,", y = ", y_now," is ", current_Obj)
end


#     last_Obj = current_Obj

# println(m)
println("Output: ")
println("P_set ", P_set)
println("M_set ", M_set)
for x_i =1:length(X_feas)
   println("x = ", X_feas[x_i],", y = ", Y_sol[x_i],"; obj = ", Obj_X[x_i]) 
end

println("\n\nVersus Optimal Sol : ")
# global P = [[1, 8], [3], [2, 9]]
global M = EnumY #[[2], [1, 3]]
numY = length(M)
numPaths = length(P)
m = Model(() -> Gurobi.Optimizer(gurobi_env))
@variable(m, x[1:numArcs], Bin)
@variable(m, y[1:numY] >= 0)
@variable(m, u)
@objective(m, Max, u)
con0 = @constraint(m, sum(y[m] for m = 1:numY) == 1)
@constraint(m, sum(d_x[a]*x[a] for a=1:numArcs) <= b_x)
# for a in [1,3,5,6,7]
#     @constraint(m, x[a] == 1)
# end
#     println("numY = ", numY)
#     println("numPaths = ", numPaths)
global Q = []
# global P_set = [[1, 8], [3], [2, 9]]
# global M_set = [[2], [1, 3]]

for i = 1:numPaths
#     global Q
#     coef_q, Q = findQgivenP(EnumY, numY, Q, P[i])

    constr[i] = @constraint(m, u <= -sum(sum(log(1-q[a]) for a in intersect(P[i],M[m]))*y[m] for m = 1:numY) + 10^6*(sum(x[a] for a in P[i])) ) 
#     @constraint(m, u <= -sum(sum(log(1-q[a]) for a in P[i])*y[m] for m = 1:numY)   ) 
end

optimize!(m)
println("Obj = ", JuMP.objective_value.(m))
println("x* = ", findall(JuMP.value.(x).>0))
for i = 1:numY
    if JuMP.value.(y[i]) > 0
        println("y* of ", M[i]," = ", JuMP.value.(y[i]))
    end
end
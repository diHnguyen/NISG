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
include("./functionGenMonitoring.jl")
global numArcs = length(arcs[:,1])
# global numArcs
# global p_exact, p_approx
global M_set = [[2,6,7]] #EnumX(arcs, b_y, numArcs, d_y)
global numY = 1 #length(M_set)
println("M_set = ", M_set)
#Add at least one path to the initial P_set
global P_set = [[1,8]]
println("P_set =", P_set )
global numPaths = 1
global Q = []
global R = 1 #000
coef_q, Q = findQgivenP(M_set, numY, Q, P_set[1])
# println("Path's probability of not being detected for each M = ", coef_q)
println("Q = ", Q)

# global ncols = 1
global cRefNum = 2000000
global constr = Array{JuMP.ConstraintRef}(undef, cRefNum)
m = Model(() -> Gurobi.Optimizer(gurobi_env))
@variable(m, 1>=x[1:numArcs]>=0 )
# @variable(m,x[1:numArcs], Bin)
@variable(m, y[1:numY] >= 0)
@variable(m, u) # >= 0)
@objective(m, Min, u)# sum(x[a] for a = 1:ncols))
con0 = @constraint(m, sum(y[m] for m = 1:numY) == 1)
# coef_q = Q[1]
@constraint(m, sum(x[a] for a=1:numArcs) <= 1)
constr[numPaths] = @constraint(m, u >= sum(y[m]*coef_q[m] for m = 1:numY) - R*(sum(x[a] for a in P_set[1])) )
# global current_Obj = 1
# global last_Obj = 2
global last_x = ones(numArcs)
global x_now = zeros(numArcs)
global iter = 0
global newPathAdded = true
global lastPath
global P, M
global P_term = false
global M_term = false
# for iter = 1:5
# while  last_x != x_now || newPathAdded == true
while P_term == false || M_term == false # newPathAdded == true
    global numArcs, M_set, numY, P_set, Q, R, numPaths, M, P, M_term, P_term
    global iter, lastPath, newPathAdded, last_x, x_now #terminate, current_Obj, last_Obj, 
    iter = iter + 1
    println("\nIter ", iter)
    optimize!(m)
#     last_Obj = current_Obj
    current_Obj = JuMP.objective_value.(m)
    y_now = JuMP.value.(y)
    x_now = JuMP.value.(x)
    lambda = zeros(numPaths)
    for i = 1:numPaths
        lambda[i] = JuMP.dual(constr[i])
    end
    println("Obj = ", current_Obj, "; x_now = ", findall(x_now.>0))
    println("y_now = ", y_now, "; M_set = ", M_set)
#     println("x_now = ", findall(x_now.>0))
    for i = 1:numPaths
        println("C",i," : ", constr[i])
    end
#     println("M_set = ", M_set)
#     println("P_set = ", P_set)
#     println("Q = ", Q)
    #Get Path Bellman
#     if iter != 4
    path, gx = shortestPath_BellmanFord(zeros(numArcs), y_now, arcs)
#     lastPath = arcsinPath
    arcsinPath = findall(path.>0)
#     if last_x == x_now && arcsinPath == lastPath
#         newPathAdded == false
#         lastPath = arcsinPath
#     end
#     println("P = ", arcsinPath)
#     println("P_set = ", P_set)
#     println("P in P_Set? ", arcsinPath in P_set)
    if (arcsinPath in P_set) == false 
        numPaths = numPaths + 1
        push!(P_set, arcsinPath)
        println("path = ", arcsinPath, "; gx = ", gx)
#         findPathsExactPr(y_now)
        println("Add new path", arcsinPath)
#         arcsinPath = [5,7]
        coef_q, Q = findQgivenP(M_set, numY, Q, arcsinPath)
        constr[numPaths] = @constraint(m, u >= sum(y[m]*coef_q[m] for m = 1:numY) - R*sum(x[a] for a in arcsinPath))
        println("Constr ", constr[numPaths])
    else
        P_term = true
        #newPathAdded = false
        #If no new path needs to be added
        #Then check if there is a candidate monitoring solution
        #If there is a new M, then add a new variable for M
        numY += 1
#         M = getM(numY)
        M = genMonitoring(lambda, P_set)
        println("new M = ", M)
        
        #Setting names and bounds on anonymous variables must follow the right naming convention
        #https://jump.dev/JuMP.jl/stable/manual/variables/#Anonymous-JuMP-variables
        if M in M_set
            println("M = ", M)
            println("M_set = ", M_set)
            M_term = true
        else
            push!(M_set, M)
            new_var = @variable(m, [numY], base_name = "y", lower_bound = 0)
            push!(y, new_var[numY])
            M_q, Q = findQgivenM(P_set, numPaths, Q, M)
            println("M_q = ", M_q)
            for i = 1:numPaths
                path = P_set[i]
#                 if isempty(intersect(M,path)) == false 
                    set_normalized_coefficient(
                        constr[i], new_var[numY], -Q[i][numY])
#                 else
#                     set_normalized_coefficient(
#                         constr[i], new_var[numY], -1)
#                 end
                set_normalized_coefficient(
                        con0, new_var[numY], 1)
            end
        end
#             @constraint(m, x[ncols] >= b[ncols])
#             @constraint(m, v >= x[ncols])
#             optimize!(m)
#         current_Obj = JuMP.objective_value.(m)
#         println("\nIter ", ncols, " : ", current_Obj)
    end
#     println("newPathAdded ", newPathAdded)

end


optimize!(m)
#     last_Obj = current_Obj
current_Obj = JuMP.objective_value.(m)
y_now = JuMP.value.(y)
x_now = JuMP.value.(x)
println(m)
println("P_set ", P_set)
println("M_set ", M_set)
println("y_now = ", y_now)
println("x_now = ", x_now)
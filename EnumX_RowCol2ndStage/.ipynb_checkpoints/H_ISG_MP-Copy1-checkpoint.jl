#Updated on 6/28
#All monitoring decisions are included in the MP
#The coefficients of each arc in the MP is exact
#The path with the least detection probability = most probability of being undetected is approximate
using JuMP
using Gurobi
using LightGraphs, Combinatorics
println("Done loading packages")
# const MOI = Gurobi.MOI
const gurobi_env = Gurobi.Env()
# setparams!(gurobi_env, Heuristics=0.0, Cuts = 0, OutputFlag = 0) #,LazyConstraints=1)
println("Done setting param")
# outfile = "C:/Users/din/Documents/GitHub/UncertainTarget/Log_Test_Output.csv"
#     io = open(outfile, "w")
# println(io,"Instance \tPath \tExact \tApprox")
# global p_exact = []
# global p_approx = []


global ins = 1# ARGS[1] # 0 #:10

# include("./Ins0.jl")
include("./Ins"*string(ins)*".jl")
# include("./TestInstances/N50_d25_"*string(ins)*".jl")
include("./functionEnumFeasX.jl")
include("./functionFindQ.jl")
include("./functionGetCost_SP.jl")

include("./functionShortestPathBellmanFord.jl")
include("./functionGenMonitoring.jl")
include("./functionLP_findM.jl")
global numArcs = length(arcs[:,1])
global X_feas = EnumX(arcs, b_x, numArcs, d_x)
global EnumY = EnumX(arcs, b_y, numArcs, d_y)
# println("EnumY = ", EnumY)
# someM = findall(d_y.<b_y)

global M_set =[EnumY[1]]
global numY = 1 
global R = 10^6

# path, gx = shortestPath_BellmanFord(q, arcs, numArcs, origin, destination)
# arcsinPath = findall(path.>0)
global P_set = [P[1]]#[arcsinPath]

global x_now, current_Obj, y_now
# println("X_feas = ", X_feas)
println("Starting M_set ", M_set, ", P_set = ", P_set)
global Y_sol = []

global numPaths = 1
global newP = true
global newM = true
global iter = 0
global tol = 10^(-6)
global theta_now = -1e12
global inner_Obj = 1-1e12 #set this to the smallest possible value of the obj val
MP = Model(() -> Gurobi.Optimizer(gurobi_env))
#Heuristics=0.0, Cuts = 0, OutputFlag = 0

# MP = direct_model(Gurobi.Optimizer(gurobi_env))
@variable(MP, x[1:numArcs], Bin)
@variable(MP, theta >= -1e12)
@objective(MP, Min, theta)
@constraint(MP, sum(d_x[a]*x[a] for a=1:numArcs) <= b_x)

# println(MP)

while theta_now + tol < inner_Obj  #newP == true || newM == true
    global numArcs, arcs, d_x, q, b_x, origin, destination, R, d_y, b_y
    global newP, newM, P_set, M_set, numPaths,numY, iter, theta_now, inner_Obj
    println("theta_now + tol - inner_Obj = ", theta_now + tol - inner_Obj)
    
    println(iter)
    iter += 1
    println("\nIteration number = ", iter)
#     println(MP)
    optimize!(MP)
    theta_now = JuMP.value.(theta)
    x_now=JuMP.value.(x)
#     x_now = callback_value.(Ref(cb_data), x)
#     theta_now = callback_value(cb_data, theta)

    println("\tx_now = ", findall(x_now.>0))
    println("\ttheta_now = ", theta_now)
    newP = false
    numPaths = length(P_set)
    numY = length(M_set)

    println("\tSolving SubProblem")
#     println("\tPre-solve:P_set = ", P_set, ", M_set = ", M_set)
    P_set, M_set, lambda, pi_, inner_Obj = solveLP_findM(x_now)
    numPaths = length(P_set)
    numY = length(M_set)
#     println("\tPost-solve:P_set = ", P_set, ", M_set = ", M_set)
    println("\tPost-solve: numPaths = ", numPaths, ", numY = ", numY)
    
    lambda_pos = findall(lambda.>0)
    println("\tlambda = ", lambda[lambda_pos], "; pi_ = ", pi_)
    println("\tpaths w lambda pos = ", lambda_pos)
    
#         dual_Obj = sum(lambda[i]*(-R)*sum(x_now[a] for a in P_set[i]) for i=1:numPaths)+pi_
#         println("dual ", dual_Obj)

    println("theta_now = ", theta_now, " vs subProb = ", inner_Obj)

    if theta_now + tol < inner_Obj 
        println("Adding constraint")
        c = @constraint(MP, theta - sum(lambda[i]*(-R)*sum(x[a] for a in P_set[i]) for i=1:numPaths) >= pi_)
        #MOI.set(MP, Gurobi.ConstraintAttribute("Lazy"), c, 2)
    end

#     if theta_now + tol >= inner_Obj
#        break 
#     end
end

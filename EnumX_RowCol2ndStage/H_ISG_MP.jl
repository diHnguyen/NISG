#Updated on 6/28
#All monitoring decisions are included in the MP
#The coefficients of each arc in the MP is exact
#The path with the least detection probability = most probability of being undetected is approximate
# function solveMP()
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


global ins = 0 #:10

include("./Ins"*string(ins)*".jl")
include("./functionEnumFeasX.jl")
include("./functionFindQ.jl")
include("./functionGetCost_SP.jl")

include("./functionShortestPathBellmanFord.jl")
include("./functionGenMonitoring.jl")
include("./functionLP_findM.jl")
global numArcs = length(arcs[:,1])
global X_feas = EnumX(arcs, b_x, numArcs, d_x)
global EnumY = EnumX(arcs, b_y, numArcs, d_y)
println("EnumY = ", EnumY)
global M_set = [EnumY[1]] #EnumX(arcs, b_y, numArcs, d_y)
global numY = 1 
global R = 10^6
global P_set = [P[1]]

global x_now, current_Obj, y_now
println("X_feas = ", X_feas)
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

#MASTER PROBLEM
# MP = direct_model(Gurobi.Optimizer(gurobi_env))
@variable(MP, x[1:numArcs], Bin)
@variable(MP, theta >= -1e12)
@objective(MP, Min, theta)
@constraint(MP, sum(d_x[a]*x[a] for a=1:numArcs) <= b_x)

println(MP)
# while outer_Obj < inner_Obj - tol #newP == true || newM == true
#     global numArcs, arcs, d_x, q, b_x, origin, destination, R, d_y, b_y
#     global newP, newM, P_set, M_set, numPaths,numY, iter, outer_Obj, inner_Obj
println(iter)
#     cb_calls = Int32[]
#     cb_calls = Cint[]
function my_callback_function(cb_data) #, cb_where::Cint)
#     function my_callback_function(cb_data::Gurobi.CallbackData, cb_where::Int32)
    global numArcs, arcs, d_x, q, b_x, origin, destination, R, d_y, b_y
    global P_set, M_set, numPaths,numY, iter, theta_now, inner_Obj
#         global iter, P_set, M_set#, lambda, pi_, inner_Obj
#         global iter_num
    iter += 1
    println("\nIteration number = ", iter)
    println(MP)
    x_now = callback_value.(Ref(cb_data), x)
    theta_now = callback_value.(Ref(cb_data), theta)

    println("\tx_now = ", x_now)
    println("\ttheta_now = ", theta_now)
    newP = false
    numPaths = length(P_set)
    numY = length(M_set)

    println("\tSolving SubProblem")
#     println("\tPre-solve:P_set = ", P_set, ", M_set = ", M_set)
    P_set, M_set, lambda, pi_, inner_Obj = solveLP_findM(x_now)
    println("\tPost-solve:P_set = ", P_set, ", M_set = ", M_set)
    println("\tlambda = ", lambda, "; pi_ = ", pi_)

    numPaths = length(P_set)
    numY = length(M_set)
#         dual_Obj = sum(lambda[i]*(-R)*sum(x_now[a] for a in P_set[i]) for i=1:numPaths)+pi_
#         println("dual ", dual_Obj)
    TOL = 1e-6
    println("theta_now = ", theta_now, " vs subProb = ", inner_Obj)
    if theta_now + TOL >= inner_Obj  #fs_x_current ≈ fm_current # we are done
        @info("No additional constraint from the subproblem")
    end
    
    if theta_now + TOL < inner_Obj 
        println("Adding constraint")
        con = @build_constraint(theta - sum(lambda[i]*(-R)*sum(x[a] for a in P_set[i]) for i=1:numPaths) >= pi_)
        println(" ", con)
#         display(con)
        MOI.submit(MP, MOI.LazyConstraint(cb_data),con)
    end
    if iter > 3
        # You can terminate the callback as follows:
        GRBterminate(backend(MP))
    end
end

MOI.set(MP, MOI.RawParameter("LazyConstraints"), 1)

MOI.set(
    MP,
    MOI.LazyConstraintCallback(),
    my_callback_function,
)

optimize!(MP)


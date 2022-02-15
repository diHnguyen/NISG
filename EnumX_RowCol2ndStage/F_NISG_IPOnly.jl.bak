println("\n==============================")
println("Integer Program: F_NISG_IPOnly")
println("==============================")
using JuMP, Gurobi, Test, Combinatorics, LightGraphs, TimerOutputs, Dates


myRun = Dates.format(now(), "HH:MM:SS")

const gurobi_env = Gurobi.Env()
# const to = TimerOutput()
ins = ARGS[1] #:10
nodeSet = "N20"
density = "d20"
dataSet = nodeSet*string("_")*density*string("_")#"N10_d10_"
println(dataSet*string(ins)*"\tRunning..."*myRun)
#     include("./TestInstances/"*fileName*".jl")
# global q = []
include("./TestInstances/p_Instances/"*nodeSet*string("/")*dataSet*string(ins)*".jl")
global numArcs = length(arcs[:,1])
# include("./functionEnumFeasX.jl")
# include("./functionFindQ.jl")
include("./functionGetCost_SP.jl")
include("./functionShortestPathBellmanFord.jl")
# include("./functionGenMonitoring.jl")

#This is the better function to find F (or Q)
include("./functionFindFValues.jl")

#Exact Row Gen:
include("./functionIP_RowGen.jl")

#Exact Col Gen:
include("./functionIP_ColGen.jl")

#Approx Col Gen
include("./functionLP_ColGen.jl")
# println("1")
include("./functionF_SG_IPOnly.jl")

#Q has numPaths elements, each element is of length = numY
global F = []
# f_MP = findf_MP([1,3], [1,3])
# println(f_MP)
# break
M_start = findall(d_x.<=b_x)
global M_set = [[M_start[1]]]#[EnumY[1]] #EnumX(arcs, b_y, numArcs, d_y)
global numY = 1 
# println("1")
path, gx = shortestPath_BellmanFord(q, arcs, numArcs, origin, destination)
arcsinPath = findall(path.>0)
global P_set = [arcsinPath]
global y_sol = []
global numPaths = 1
global newP = true
global newM = true
global iter = 0
# println("1")
global R = sum(-log(1-p[a]) for a = 1:numArcs)
println(R)
println(P_set)
println(M_set)
F = updateF(F, P_set, M_set, numY, numPaths)
println(d_x[[2,5]])
# conRefNum = 20# max(200,2*length(P_set))
# con = Array{JuMP.ConstraintRef}(undef, conRefNum)
to = TimerOutput()
global k = 0
global rootOptGap = 0
model = direct_model(Gurobi.Optimizer())
#     set_optimizer_attribute(model, "Presolve", 0)
#     set_optimizer_attribute(model, "Precrush", 1)
#     set_optimizer_attribute(model, "DualReductions", 0)
#     set_optimizer_attribute(model, "PreQLinearize", 0)
#     set_optimizer_attribute(model, "Heuristics", 0)
#     set_optimizer_attribute(model, "MIPGap", 0)
set_optimizer_attribute(model, "LazyConstraints", 1)
set_optimizer_attribute(model, "OutputFlag", 0)
# set_optimizer_attribute(model, "TimeLimit", 3600.0)
global loc = "./Documents/GitHub/UncertainTarget/EnumX_RowCol2ndStage" #"./Documents/GitHub/UncertainTarget/EnumX_RowCol2ndStage/TestInstances/"
# println("LogFile", loc*dataSet*string(ins)*".log")
# set_optimizer_attribute(model, "LogFile", loc*dataSet*string(ins)*".log")
# println("1")
# Presolve = 0,
#         PreCrush = 1,
#         Heuristics = 0,

@variable(model, x[1:numArcs], Bin)
@variable(model, theta >= 0)#1e6)
@objective(model, Min, theta)
@constraint(model, sum(d_x[a]*x[a] for a=1:numArcs) <= b_x)
global last_theta = -1
global start = time()
cb_calls = Cint[]

global t_exactRow = 0.0
global t_exactCol = 0.0
global MIPGAP = 0
global lambda = []
global lambda_pos = []
# println("2")
function my_callback_function(cb_data, cb_where::Cint)
    global numArcs, arcs, d_x, q, b_x, origin, destination, R, d_y, b_y, y_sol
    global newP, newM, P_set, M_set, numPaths,numY, iter, theta_now, inner_Obj
    global k, last_theta, t_exactRow, t_exactCol,lambda_pos,lambda, MIPGAP, nodeCount
    # You can reference variables outside the function as normal
    push!(cb_calls, cb_where)
    iter += 1
#      println("iter ", iter)
#     println("cb_where = ", cb_where)
#     println("GRB_CB_MIPSOL = ", GRB_CB_MIPSOL)
#     println("GRB_CB_MIPNODE = ", GRB_CB_MIPNODE)
    # You can query a callback attribute using GRBcbget
#     MIPGAP = 0
#     if cb_where == GRB_CB_MIPSOL #|| cb_where == GRB_CB_MIPNODE
        
#     end
    if cb_where == GRB_CB_MIPNODE
        resultN = Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_NODCNT, resultN)
        nodeCount = resultN[]
#         println("\n\tMIPNODE_NODCNT ", nodeCount)
        if nodeCount == 0
            resultO = Ref{Cdouble}()
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_OBJBST, resultO)
            obj = resultO[]
#             println("\tMIPNODE_OBJBST ", obj)
            resultO = Ref{Cdouble}()
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_OBJBND, resultO)
            objBound = resultO[]
#             println("\tMIPNODE_OBJBND ", objBound)
            MIPGAP = (obj - objBound)/obj*100
#             println("\tMIPGAP = ", MIPGAP)
#             if resultP[] != GRB_OPTIMAL
#                 return  # Solution is something other than optimal.
#             end
        end
    end
    if cb_where == GRB_CB_MIPSOL
#         resultP = Ref{Cint}()
#         GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP) 
        
        resultN = Ref{Cdouble}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_NODCNT, resultN)
        nodeCount = resultN[]
#         println("\n\tMIPSOL_NODCNT ", nodeCount)
        if nodeCount == 0
            resultO = Ref{Cdouble}()
            GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_OBJBST, resultO)
            obj = resultO[]
#             println("\tMIPSOL_OBJBST ", obj)
            resultO = Ref{Cdouble}()
            GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_OBJBND, resultO)
            objBound = resultO[]
#             println("\tMIPNODE_OBJBND ", objBound)
            MIPGAP = (obj - objBound)/obj*100
#             println("\tMIPGAP = ", MIPGAP)
#             if resultP[] != GRB_OPTIMAL
#                 return  # Solution is something other than optimal.
#             end
        end
        
        # Before querying `callback_value`, you must call:
        Gurobi.load_callback_variable_primal(cb_data, cb_where)
        x_val = callback_value.(Ref(cb_data), x)
        theta_val = callback_value(cb_data, theta)
#          println("x_val = ", findall(x_val.>0), "\t theta_val", theta_val)
#          println("M = ", M_set)
#          println("P = ", P_set)
        #Solve sub problem
#         println("3")
        
#         P_set, M_set, lambda, pi_, inner_Obj, y_sol, findF = H_SG(x_val)
#         if findF == true
#         println("preSG")
        lambda, pi_, inner_Obj, y_sol = F_SG(x_val)
#         println("postSG")
#         end
        lambda_pos = findall(lambda.>0)
        numPaths = length(P_set)
        numY = length(M_set)
#         println("\ninner_Obj = ", inner_Obj)
#         println("lambda ", lambda)
#         println("pi_ ", pi_)
#     println("4")
#         
#         println(Q)
        TOL = 10^(-5)
        if (inner_Obj - theta_val) > TOL  
            k+=1 
            con_ = @build_constraint(theta >= -sum(lambda[i]*sum(x[a] for a in P_set[i]) for i=1:numPaths) + pi_)
#             println(con_)
            w_con = MOI.submit(model, MOI.LazyConstraint(cb_data), con_)
#                 println("Typeof = ", typeof(w_con))
#             MOI.set(model, Gurobi.Lazy(), con[k], 2)
#             MOI.set(model, Gurobi.ConstraintAttribute("lazy"), con, 2)
            end
    end
    return
end
#     You _must_ set this parameter if using lazy constraints.
MOI.set(model, MOI.RawParameter("LazyConstraints"), 1)
MOI.set(model, Gurobi.CallbackFunction(), my_callback_function)
optimize!(model)
    
#     println(model)
totalTime = time() - start
x_sol = JuMP.value.(x)
theta_sol = JuMP.value.(theta)
y_pos = findall(y_sol .> 0)
nodeCount = MOI.get(model, MOI.NodeCount())


timesFile = open(loc*"/TestInstances/p_Instances_Output/IP_Only/O_"*dataSet*string(ins)*"_IPOnly.txt", "a")
println(timesFile, "dataSet = ", dataSet)
println(timesFile, "ins = ", ins)
println(timesFile, "theta = ", theta_sol)
println(timesFile, "totalTime = ", totalTime)
println(timesFile, "t_RowGen = ", TimerOutputs.time(to["IP_RowGen"])/10^9)
println(timesFile, "n_RowGen = ", TimerOutputs.ncalls(to["IP_RowGen"]))
println(timesFile, "t_ColGen = ", TimerOutputs.time(to["IP_ColGen"])/10^9)
println(timesFile, "n_ColGen = ", TimerOutputs.ncalls(to["IP_ColGen"]))
println(timesFile, "x = ", x_sol)
println(timesFile, "x_i = ", findall(x_sol.>0))
println(timesFile, "y = ", y_sol[y_pos])
println(timesFile, "lambda = ", lambda[lambda_pos])
println(timesFile, "M_pos = ", M_set[y_pos])
println(timesFile, "P_pos = ", P_set[lambda_pos])
println(timesFile, "M_set = ", M_set)
println(timesFile, "P_set = ", P_set)
println(timesFile,"Node count = ", nodeCount)
println(timesFile,"MIPGAP_root = ", MIPGAP)
close(timesFile)

println("dataSet ", dataSet, ins)
println("totalTime = ", totalTime)
println("\tt_RowGen = ", TimerOutputs.time(to["IP_RowGen"])/10^9)
println("\tn_RowGen = ", TimerOutputs.ncalls(to["IP_RowGen"]))
println("\tt_ColGen = ", TimerOutputs.time(to["IP_ColGen"])/10^9)
println("\tn_ColGen = ", TimerOutputs.ncalls(to["IP_ColGen"]))
println("x = ", x_sol)
println("x_i = ", findall(x_sol.>0))
println("theta = ", theta_sol)
println("M_set = ", length(M_set))
println("y = ", y_sol[y_pos])
println("lambda = ", lambda[lambda_pos])
println("M_pos = ", M_set[y_pos])
println("P_pos = ", P_set[lambda_pos])
println("M_set = ", M_set)
println("P_set = ", P_set)
println("Node count = ", nodeCount)
println("MIPGAP_root = ", MIPGAP)
using JuMP, Gurobi, Test, Combinatorics, LightGraphs, TimerOutputs
const gurobi_env = Gurobi.Env()
ins = ARGS[1] #:10
dataSet = "N50_d10_"
println("Running...")
#     include("./TestInstances/"*fileName*".jl")
# global q = []
include("./TestInstances/Prelim_Ins/"*dataSet*string(ins)*".jl")
global numArcs = length(arcs[:,1])
include("./functionEnumFeasX.jl")
include("./functionFindQ.jl")
include("./functionGetCost_SP.jl")
include("./functionShortestPathBellmanFord.jl")
include("./functionGenMonitoring.jl")

#This is the better function to find Q
include("./functionUpdateQ.jl")

#Exact Row Gen:
include("./functionIP_RowGen.jl")

#Exact Col Gen:
include("./functionIP_ColGen.jl")

include("./functionF_SG.jl")

#Q has numPaths elements, each element is of length = numY
global Q = []
# f_MP = findf_MP([1,3], [1,3])
# println(f_MP)
# break
M_start = findall(d_x.<=b_x)
global M_set = [[M_start[1]]]#[EnumY[1]] #EnumX(arcs, b_y, numArcs, d_y)
global numY = 1 
path, gx = shortestPath_BellmanFord(q, arcs, numArcs, origin, destination)
arcsinPath = findall(path.>0)
global P_set = [arcsinPath]
global y_sol = []
global numPaths = 1
global newP = true
global newM = true
global iter = 0
println(P_set)
println(M_set)
Q = updateQ(Q, P_set, M_set, numY, numPaths)
# conRefNum = 20# max(200,2*length(P_set))
# con = Array{JuMP.ConstraintRef}(undef, conRefNum)
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
# set_optimizer_attribute(model, "OutputFlag", 0)

global loc = "./Documents/GitHub/UncertainTarget/EnumX_RowCol2ndStage/TestInstances/"
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
global MIPGAP = 0
# println("2")
function my_callback_function(cb_data, cb_where::Cint)
    global numArcs, arcs, d_x, q, b_x, origin, destination, R, d_y, b_y, y_sol
    global newP, newM, P_set, M_set, numPaths,numY, iter, theta_now, inner_Obj
    global k, last_theta
    # You can reference variables outside the function as normal
    push!(cb_calls, cb_where)
    iter += 1
#     println("iter ", iter)
    
    # You can query a callback attribute using GRBcbget
    if cb_where == GRB_CB_MIPSOL #|| cb_where == GRB_CB_MIPNODE
        resultP = Ref{Cint}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP) 
        MIPGAP = 0
        if cb_where == GRB_CB_MIPNODE
            resultP = Ref{Cdouble}()
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_NODCNT, resultP)
            node = resultP[]
#             println("\n\tMIPNODE_NODCNT ", node)
            if node == 0
                resultP = Ref{Cdouble}()
                GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_OBJBST, resultP)
                obj = resultP[]
#                 println("\tMIPNODE_OBJBST ", obj)
                resultP = Ref{Cdouble}()
                GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_OBJBND, resultP)
                objBound = resultP[]
#                 println("\tMIPNODE_OBJBND ", objBound)
                MIPGAP = (obj - objBound)/obj*100
#                 println("\tMIPGAP = ", MIPGAP)
                if resultP[] != GRB_OPTIMAL
                    return  # Solution is something other than optimal.
                end
            end
        end
        # Before querying `callback_value`, you must call:
        Gurobi.load_callback_variable_primal(cb_data, cb_where)
        x_val = callback_value.(Ref(cb_data), x)
        theta_val = callback_value(cb_data, theta)
#         println("x_val = ", findall(x_val.>0), "\t theta_val", theta_val)
        #Solve sub problem
#         println("3")
        P_set, M_set, lambda, pi_, inner_Obj, y_sol = F_SG(x_val)
        numPaths = length(P_set)
        numY = length(M_set)
#         println("lambda ", lambda)
#     println("4")
#         println("\ninner_Obj = ", inner_Obj)
#         println(Q)
        TOL = 10^(-3)
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

timesFile = open(loc*"O_"*dataSet*string(ins)*".txt", "a")
# println(timesFile, dataSet*prob, "; Ins ", ins, "; Time ", totalTime, "; theta ", theta_sol,"; x_sol ", findall(x_sol.>0),"; M_select ", M_set[y_pos], "; y_sol ", y_sol[y_pos], "; Iter ", iter, "; M_len ", length(M_set), "; P_len ", length(P_set), "; NodeCount ", nodeCount)
println(timesFile, "dataSet = ", dataSet)
println(timesFile, "ins = ", ins)
println(timesFile, "theta = ", theta_sol)
println(timesFile, "totalTime = ", totalTime)
println(timesFile, "x = ", x_sol)
println(timesFile, "x_i = ", findall(x_sol.>0))
println(timesFile, "y = ", y_sol[y_pos])
println(timesFile, "M_pos = ", M_set[y_pos])
println(timesFile, "M_set = ", M_set)
println(timesFile, "P_set = ", P_set)
println(timesFile,"Node count = ", nodeCount)
println(timesFile,"MIPGAP_root = ", MIPGAP)
close(timesFile)

println("dataSet = ", dataSet)
println("totalTime = ", totalTime)
println("x = ", x_sol)
println("x_i = ", findall(x_sol.>0))
println("theta = ", theta_sol)
println("M_set = ", length(M_set))
println("y = ", y_sol[y_pos])
println("M_pos = ", M_set[y_pos])
println("M_set = ", M_set)
println("P_set = ", P_set)
println("Node count = ", nodeCount)
println("MIPGAP_root = ", MIPGAP)
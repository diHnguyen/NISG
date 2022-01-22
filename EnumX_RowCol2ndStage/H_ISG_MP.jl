#Relevant as of Jan 11th.
using JuMP, Gurobi, Test, Combinatorics, LightGraphs, TimerOutputs

# function solveMP()
const gurobi_env = Gurobi.Env()
# setparams!(gurobi_env, OutputFlag = 0)

ins = ARGS[1] #:10
dataSet = "N50_d40_"
println("Running...")
global prob = "L" #"L"
#     include("./TestInstances/"*fileName*".jl")
global q
include("./TestInstances/"*dataSet*string(ins)*".jl")
global numArcs = length(arcs[:,1])
global lb = 1
global ub = 999
# println("1")
function setQ(prob)
    global lb, ub, numArcs, q
    if prob == "L"
        lb = 1
        ub = 333
    elseif prob == "M"
        lb = 334
        ub = 666
    elseif prob == "H"
        lb = 667
        ub = 999
    end
    new_q = rand(lb:ub, numArcs)./1000
       
    return new_q
end
# println("2")
if prob != ""
    q = setQ(prob)
end
println("q = ", q)
# q = [0.28, 0.23, 0.033, 0.257, 0.208, 0.037, 0.009, 0.269, 0.302, 0.167, 0.098, 0.288, 0.187]
# println("3")
#     global numArcs, arcs, d_x, q, b_x, origin, destination, R, d_y, b_y
#     global newP, newM, P_set, M_set, numPaths,numY
    
    include("./functionEnumFeasX.jl")
    include("./functionFindQ.jl")
    include("./functionGetCost_SP.jl")

    include("./functionShortestPathBellmanFord.jl")
    include("./functionGenMonitoring.jl")
    include("./functionLP_findM.jl")
    

        
#     global X_feas = EnumX(arcs, b_x, numArcs, d_x)
#     global EnumY = EnumX(arcs, b_y, numArcs, d_y)
    M_start = findall(d_x.<=b_x)
    global M_set = [[M_start[1]]]#[EnumY[1]] #EnumX(arcs, b_y, numArcs, d_y)
    global numY = 1 
    global R = -sum(log10(1-q_) for q_ in q)  #10^6#0^6
    println("R = ", R)
    path, gx = shortestPath_BellmanFord(q, arcs, numArcs, origin, destination)
    arcsinPath = findall(path.>0)
    global P_set = [arcsinPath]
#     global P_set = [P[1]]

#     global x_now, current_Obj, y_now

#     println("X_feas = ", X_feas)
#     println("Starting M_set ", M_set, ", P_set = ", P_set)

global y_sol = []
#     global 
global numPaths = 1
global newP = true
global newM = true
global iter = 0
#     global theta_now = -1e12
#     global inner_Obj = 1-1e12
conRefNum = 20# max(200,2*length(P_set))
con = Array{JuMP.ConstraintRef}(undef, conRefNum)
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
# println("Here")
set_optimizer_attribute(model, "LogFile", loc*"L_"*dataSet*prob*"_"*string(ins)*".log")
# Presolve = 0,
#         PreCrush = 1,
#         Heuristics = 0,

@variable(model, x[1:numArcs], Bin)
@variable(model, theta >= 1-R)#1e6)
@objective(model, Min, theta)
@constraint(model, sum(d_x[a]*x[a] for a=1:numArcs) <= b_x)
    
#     @variable(model, 0 <= x <= 2.5, Int)
#     @variable(model, 0 <= y <= 2.5, Int)
#     @objective(model, Max, y)
#     println("\nM_set = ", M_set," \t P_set = ", P_set)
global last_theta = -1e6 - 1
global start = time()
cb_calls = Cint[]
function my_callback_function(cb_data, cb_where::Cint)
    global numArcs, arcs, d_x, q, b_x, origin, destination, R, d_y, b_y, y_sol
    global newP, newM, P_set, M_set, numPaths,numY, iter, theta_now, inner_Obj
    global k, last_theta
        # You can reference variables outside the function as normal
        push!(cb_calls, cb_where)
        iter += 1
    
#         println("\niter = ", iter)
#         println("\nM_set = ", M_set," \t P_set = ", P_set)
        # You can select where the callback is run
#         if cb_where != GRB_CB_MIPSOL && cb_where != GRB_CB_MIPNODE
#             return
#         end
        
        # You can query a callback attribute using GRBcbget
        if cb_where == GRB_CB_MIPSOL #|| cb_where == GRB_CB_MIPNODE
            resultP = Ref{Cint}()
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
#         println(, " ", MOI.get(model, MOI.NodeCount())
#         println("cb_data ", cb_data)
#         println("cb_where ", cb_where)
#         println("GRB_CB_MIPNODE_STATUS ", GRB_CB_MIPNODE_STATUS)   
        # Before querying `callback_value`, you must call:
        Gurobi.load_callback_variable_primal(cb_data, cb_where)
#         if resultP[] != GRB_OPTIMAL
#             return  # Solution is something other than optimal.
#         end
#         if resultP == GRB_OPTIMAL
#             println("\niter = ", iter)
            x_val = callback_value.(Ref(cb_data), x)
            theta_val = callback_value(cb_data, theta)
#             println("x_val = ", findall(x_val.>0))#, "\t theta_val", theta_val)
            P_set, M_set, lambda, pi_, inner_Obj, y_sol = solveLP_findM(x_val)
#             println("M = ", length(M_set), "\t P = ", length(P_set)) 
#             println("P = ", P_set)
#             println("theta_val = ", theta_val, "\tinner_Obj = ", inner_Obj)
        
#             println("lambda = ", lambda, "\t pi_ = ", pi_)
            numPaths = length(P_set)
            numY = length(M_set)

            TOL = 10^(-3)
#             println("(inner_Obj - theta_val) = ", (inner_Obj - theta_val))
            
#             println("TOL = ", TOL)
            if (inner_Obj - theta_val) > TOL  
                k+=1 
                con_ = @build_constraint(theta >= sum(lambda[i]*(-R)*sum(x[a] for a in P_set[i]) for i=1:numPaths) + pi_)
#                 println("Added con")
#                 println("Type = ", typeof(con_))
#             cut_vio = (sum(lambda[i]*(-R)*sum(x_val[a] for a in P_set[i]) for i=1:numPaths) + pi_) - theta_val
#                 println("Add con : cut violation = ", cut_vio)
#                 println(theta, ">=", sum(lambda[i]*(-R)*sum(x[a] for a in P_set[i]) for i=1:numPaths) + pi_)
                
                w_con = MOI.submit(model, MOI.LazyConstraint(cb_data), con_)
#                 println("Typeof = ", typeof(w_con))
#             MOI.set(model, Gurobi.Lazy(), con[k], 2)
#             MOI.set(model, Gurobi.ConstraintAttribute("lazy"), con, 2)
            end
#         if last_theta > theta_val
#             GRBterminate(backend(model))
#         else
#             last_theta = theta_val
#         end
#         end 
    end
        
       
        
#         if iter > 12
#             # You can terminate the callback as follows:
#             GRBterminate(backend(model))
#         end
    return
end
#     You _must_ set this parameter if using lazy constraints.
    MOI.set(model, MOI.RawParameter("LazyConstraints"), 1)
    MOI.set(model, Gurobi.CallbackFunction(), my_callback_function)
#     println(model)
    optimize!(model)
    
#     println(model)
    totalTime = time() - start
    x_sol = JuMP.value.(x)
    theta_sol = JuMP.value.(theta)
    println("dataSet ", dataSet*prob)
    println("x = ", x_sol)
    println("x_i = ", findall(x_sol.>0))
    println("theta = ", theta_sol)
    println("M_set = ", length(M_set))
    y_pos = findall(y_sol .> 0)
    println("y = ", y_sol[y_pos])
    println("M_pos = ", M_set[y_pos])
    nodeCount = MOI.get(model, MOI.NodeCount())
    println("Node count = ", nodeCount)

timesFile = open(loc*"O_"*dataSet*prob*".txt", "a")
println(timesFile, dataSet*prob, "; Ins ", ins, "; Time ", totalTime, "; theta ", theta_sol,"; x_sol ", findall(x_sol.>0),"; M_select ", M_set[y_pos], "; y_sol ", y_sol[y_pos], "; Iter ", iter, "; M_len ", length(M_set), "; P_len ", length(P_set), "; NodeCount ", nodeCount)
close(timesFile)

if prob != ""
    qFile = open(loc*"q_"*dataSet*prob*".txt", "a")
    println(qFile, dataSet*prob, "; Ins ", ins, "; q ", q)
    close(qFile)
end
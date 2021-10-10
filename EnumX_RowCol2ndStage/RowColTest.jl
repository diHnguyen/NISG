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
# include("./functionFindQ.jl")
include("./functionShortestPathBellmanFord.jl")
# include("./findExactUndetectionPr.jl")
include("./functionGenMonitoring.jl")
# include("./functionSolve2ndStage.jl")
function solve2ndStage(x_fix, arcs, numArcs, numPaths, numY, P_set, M_set)
#     global constr, conNum
    cRefNum = 200
    constr = Array{JuMP.ConstraintRef}(undef, cRefNum)
    m = Model(() -> Gurobi.Optimizer(gurobi_env))
#     @variable(m, 1>=x[1:numArcs]>=0 )
    @variable(m, y[1:numY] >= 0)
    @variable(m, u)
    #V1
    @objective(m, Min, u)
    
    #V2
#     @objective(m, Max, u)
    con0 = @constraint(m, sum(y[m] for m = 1:numY) == 1)
    
    x_now = zeros(numArcs)
    for a in x_fix
       x_now[a] = 1 
    end
    for i = 1:numPaths
        #V1
        constr[i] = @constraint(m, u >= -sum(sum(log(1-q[a]) for a in intersect(P_set[i],M_set[m]))*y[m] for m = 1:numY) - 10^6*(sum(x_now[a] for a in P_set[i])) )

        #V2
#         constr[i] = @constraint(m, u <= -sum(sum(log(1-q[a]) for a in intersect(P_set[i],M_set[m]))*y[m] for m = 1:numY) + 10^6*(sum(x_now[a] for a in P_set[i])) ) 

    end
    
    iter = 0
    newPathAdded = true
#     lastPath
#     P, M
    newP = true
    newM = true

    while newP == true || newM == true 
        newP = false
        newM = false
        iter = iter + 1
#         println("\nIter ", iter)
        println("\nM_set = ", M_set)
        println("P_set = ", P_set)
#         println(m)
        set_optimizer(m, ()-> Gurobi.Optimizer(gurobi_env))
        optimize!(m)

        current_Obj = JuMP.objective_value.(m)
        y_now = JuMP.value.(y)
        println("y = ", y_now) #findall(y_now.>0))
        lambda = zeros(numPaths)
#         println("1")
        for i = 1:numPaths
            lambda[i] = JuMP.dual(constr[i])
        end
        path, gx = shortestPath_BellmanFord(x_now, y_now, arcs, q, numArcs, M_set, numY)
        arcsinPath = findall(path.>0)
        println("arcsinPath ", arcsinPath," cost ", gx)
        if (arcsinPath in P_set) == false 
            numPaths = numPaths + 1
            push!(P_set, arcsinPath)
            println("Add path ", arcsinPath)
            newP = true
#             coef_q, Q = findQgivenP(M_set, numY, Q, arcsinPath)
            #V1
            constr[numPaths] = @constraint(m, u >= -sum(sum(log(1-q[a]) for a in intersect(P_set[numPaths], M_set[m]))*y[m] for m = 1:numY) - 10^6*(sum(x_now[a] for a in P_set[numPaths])) )
            #V2
#             constr[numPaths] = @constraint(m, u <= -sum(sum(log(1-q[a]) for a in intersect(P_set[numPaths],M_set[m]))*y[m] for m = 1:numY) + 10^6*(sum(x_now[a] for a in P_set[numPaths])) )
        else
            M = genMonitoring(lambda, P_set)
            
            #Setting names and bounds on anonymous variables must follow the right naming convention
            #https://jump.dev/JuMP.jl/stable/manual/variables/#Anonymous-JuMP-variables
            if M in M_set || isempty(M) == true
                M_term = true
            else
                println("Add M ", M)
                newM = true
                numY += 1
                push!(M_set, M)
                new_var = @variable(m, [numY], base_name = "y", lower_bound = 0)
                push!(y, new_var[numY])
#                 M_q, 
#                 Q = findQgivenM(P_set, numPaths, Q, M)
#                 println("M_q = ", M_q)
                for i = 1:numPaths
                    path = P_set[i]
#                     println("P = ", path, "; M = ", M_set[numY])
                    McapP = intersect(P_set[i], M_set[numY])
                    if isempty(McapP) == false
                        set_normalized_coefficient(constr[i], new_var[numY], sum(log(1-q[a]) for a in McapP) )
                    else
                        set_normalized_coefficient(constr[i], new_var[numY], 0)
                    end
                    set_normalized_coefficient(con0, new_var[numY], 1)
                end
            end
        end
    end
    
    optimize!(m)
#     println(m)
    current_Obj = JuMP.objective_value.(m)
    y_now = JuMP.value.(y)
#     x_now = JuMP.value.(x)
#     println(m)
#     println("Obj = ", current_Obj)
#     println("x_now = ", findall(x_now.>0), "; ", x_now[x_now.>0])
    println("M_set ", M_set)
#     println("y_now = ", y_now)
    println("P_set ", P_set)

    return M_set, P_set, current_Obj, y_now, numY, numPaths
end
global numArcs = length(arcs[:,1])
# global numArcs
# global p_exact, p_approx
global X_feas = EnumX(arcs, b_x, numArcs, d_x)

global EnumY = EnumX(arcs, b_y, numArcs, d_y)
println("EnumY = ", EnumY)
for x_i = 4:4#length(X_feas)
    global arcs, d_x, q, b_x, origin, destination, P, d_y, b_y, EnumY
    x_now = X_feas[x_i]
    println("\n\nx = ", x_now)
    M_set = [EnumY[1]]
    numY = 1 #length(M_set)
    P_set = [[1,8]]
    numPaths = 1
    
    cRefNum = 200
    constr = Array{JuMP.ConstraintRef}(undef, cRefNum)
    global newPathAdded = true
    global lastPath
    
#     P_term = false
#     M_term = false
    M_set, P_set, current_Obj, y_now, numY, numPaths = solve2ndStage(x_now, arcs, numArcs, numPaths, numY, P_set, M_set)
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
#     for i = 1:numPaths
#         log_q = []
#         for m = 1:numY
#             mySet = intersect(P_set[i],M_set[m])
#             println("Path ", P_set[i], ", Monitor", M_set[m])
#             println(mySet, " : ", q[mySet])
#             myq = 0
#             if isempty(mySet) == false
#                 myq = sum(log(1-q[a]) for a in mySet)
#                 println("myq = ", myq)
#             end
#             push!(log_q, myq)
#         end
#         println("Path ", i, ":",log_q)
    #     end 
    for i = 1:numPaths
         @constraint(MP, u >= -sum(
                sum(log(1-q[a]) for a in intersect(P_set[i],M_set[m]))*y[m] for m = 1:numY) - 10^6*(sum(x[a] for a in P_set[i])) ) 
    end
    optimize!(MP)
    println("\nSolve once given P_set and M_set:")
    println("P = ", P_set)
    println("M = ", M_set)
#     println(MP)
    println("2. Obj Val for x = ", findall(JuMP.value.(x).>0),", y = ", JuMP.value.(y)," is ", JuMP.objective_value.(MP))
end

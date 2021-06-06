using JuMP
using Gurobi
using Combinatorics



global gurobi_env = Gurobi.Env()
setparams!(gurobi_env, Heuristics=0.0, Cuts = 0, OutputFlag = 0)
global ins = 3
include("./Ins"*string(ins)*".jl")
println("Here")
println("Arcs = ", arcs)
global numArcs = length(arcs[:,1])
global Len = numArcs
include("C:/Users/din/Documents/GitHub/IntCVaR/functionGbound.jl")
include("C:/Users/din/Documents/GitHub/UncertainTarget/functionEnumFeasX.jl")

global X 

function solve2ndStage(x_ind, arcs, numArcs, q, d, origin, destination)
    global gurobi_env
#     x = zeros(numArcs)
#     for a in x_ind
#        x[a] = 1 
#     end
#     println("x = ", x[x.>0])
    m = Model(() -> Gurobi.Optimizer(gurobi_env))

    @variable(m, z[1:numArcs]>=0)
    @variable(m, v >= 0)
    @variable(m, s >= 0)
    
    #Flow Constraints:
#     println("Flow")
    for i = 1:(destination-1)
        inflow = findall(arcs[:,2].==i)
        outflow = findall(arcs[:,1].==i)
#         println("Fine")
        if i == origin
            con = @constraint(m, sum(z[a] for a in outflow) + s == 1)
#             println(con)
        else
            con = @constraint(m, sum(z[a] for a in outflow) - sum(z[a] for a in inflow) == 0)
#             println(con)
        end
    end
#     println("a")
    #Arc Capacity
    for a = 1:numArcs
        if q[a] > 0
            con = @constraint(m, z[a] <= v/q[a])
#             println(con)
#         else
#             con = @constraint(m, z[a] <= v*(10^6))
#             println(con)
        end
        if a in x_ind
            con = @constraint(m, z[a] <= 0)#(10^6)*(1-x[a])) 
        
#             println(con)
        end
    end
#     @constraint(m, sum(d[i]*x[i] for i =1:numArcs)<= b)
#     println(findall(P.>0))
#     @constraint(m, v <= sum(x[a] + q[a]*y[a] for a in findall(P.>0)))
#     @constraint(m, sum(y[a] for a =1:numArcs) == 1)
    @objective(m, Min, v + 10^6*s)
#     println(m)
    optimize!(m)
    Obj_Val = JuMP.objective_value.(m)
    current_z = JuMP.value.(z)
    println("z = ", findall(current_z.>0), " ; ", current_z[current_z.>0])
    f = 1/(sum(q[a] for a = 1:numArcs))
    y = f./q
#     println("Optimal y: ",y)
    current_v = JuMP.value.(v)
    println("Solved 2ndStage ", Obj_Val, " ; s = ", JuMP.value.(s))
    return Obj_Val, current_z, current_v
end

function retrieveMinCut(z, v, numArcs)
#     println("In retrieveMinCut")
    global gurobi_env, destination, origin
#     println("2")
    m = Model(() -> Gurobi.Optimizer(gurobi_env))
    
    P = findall(z.>0)
#     println("P = ", P)
#     sizeP = length(P)
    cap = Array{JuMP.ConstraintRef}(undef, numArcs)
    flow = Array{JuMP.ConstraintRef}(undef, destination)
    @variable(m, y[1:numArcs]>=0)
    
    #Flow Constraints:
    for i = 1:(destination-1)
        inflow = findall(arcs[:,2].==i)
        outflow = findall(arcs[:,1].==i)
        if i != origin && i != destination
            flow[i] = @constraint(m, sum(y[a] for a in outflow if a in P) - sum(y[a] for a in inflow if a in P) == 0)
#             println(con)
        end
    end
    #Capacity Constraints
    for a = 1:numArcs
        if a in P && q[a] > 0
            cap[a] = @constraint(m, y[a] <= v/q[a])
        end
    end
#     println("1")
    outflow = findall(arcs[:,1].==1)
    @objective(m, Max, sum(y[a] for a in outflow if a in P))
    optimize!(m)
    
    cap_dual = zeros(numArcs)
    flow_dual = zeros(destination)
    if termination_status(m) == MOI.OPTIMAL
        Obj_Val = JuMP.objective_value.(m)
        for a in P
            if q[a] > 0
                cap_dual[a] = JuMP.dual(cap[a])
            end
        end
    end
    
#     println(Obj_Val)
#     println("flow_dual = ", cap_dual)
    return findall(cap_dual.<0)
end


X = EnumX(arcs,b,numArcs,d)
println("Set of maximally packed x-solution = ", X)
global x_sol = zeros(numArcs)
global current_z, current_v


outfile = "C:/Users/din/Documents/GitHub/UncertainTarget/TestInstances/Output_Ins"*string(ins)*".jl"
io = open(outfile, "a") #do io
println(io, "#Verify")
println(io, "Set of maximally packed x-solution = ", X)
if length(X) > 0
    for i = 1:length(X)
        x_sol = X[i]
        println("\nx = ", x_sol)
        Obj, current_z, current_v =  solve2ndStage(x_sol, arcs, numArcs, q, d, origin, destination)
        
        if Obj <= 1
            A_C = retrieveMinCut(current_z, current_v, numArcs)
            if length(A_C) > 0
                f_inv = 1/sum(1/q[a] for a in A_C)
                println(io, "x = ", x_sol, " ; v** = ",f_inv ," ; minCut = ", A_C, " ; ", [1/q[a]*f_inv for a in A_C] )
                println("W&W says... v** = f_inv = ", f_inv)
                println("minCut = ", A_C)
                println("y Monitoring = ", [1/q[a]*f_inv for a in A_C])
            else
                println(io, "x = ", x_sol, " ; The evader always reaches the target destination.")
                println("There's always a feasible path (e.g., q[a]=0 for all a in P)")
            end
        else
            println("No monitoring is required.")
            println(io,"x = ", x_sol, " ; No monitoring is required.")
        end     
    end
else
    x_sol = zeros(numArcs)
    println("\nx = ", x_sol)
    Obj, current_z, current_v =  solve2ndStage(x_sol, arcs, numArcs, q, d, origin, destination)
    A_C = retrieveMinCut(current_z, current_v, numArcs)
    f_inv = 1/sum(1/q[a] for a in A_C)
    println(io, "x = ", x_sol, " ; v** = ",f_inv ," ; minCut = ", A_C, " ; ", [1/q[a]*f_inv for a in A_C] )
    
    println("W&W says... v** = f_inv = ", f_inv)
    println("minCut = ", A_C)
    println("y Monitoring = ", [1/q[a]*f_inv for a in A_C])
end

close(io)
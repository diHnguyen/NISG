using Gurobi
using JuMP
using LightGraphs

global arcs
global d
global q

# global b = 5

include("C:/Users/din/Documents/GitHub/IntCVaR/functionGbound.jl")
include("C:/Users/din/Documents/GitHub/UncertainTarget/functionEnumFeasX.jl")

function longestPathLP(arcs, numArcs, costs)
    global gurobi_env, origin, destination
    
    #Taking the negative of log:
    c = zeros(numArcs)
    for i = 1:numArcs
        if costs[i] != 0
            c = -log(costs[i])
        else
            c = 10^6
        end
    end
    m = Model(() -> Gurobi.Optimizer(gurobi_env))
    
    @variable(m, z[1:numArcs], Bin)
    
    #FlowConstraints
    for i = 1:destination
        outflow = findall(arcs[:,1].==i)
        if i == origin
            @constraint(m, sum(z[i] for i in outflow) == 1)
        elseif i != destination
            inflow = findall(arcs[:,2].==i)
            @constraint(m, sum(z[i] for i in outflow)-sum(z[i] for i in inflow) == 0)
        end
    end
    
    @objective(m, Max, costs[i]*z[i] for i = 1:numArcs)
    optimize!(m)
    z_sol = JuMP.value.(z)
    P = findall(z_sol.>0)
    return P 
end

global origin = 1
global destination = 6
ins = 1

# for ins = 1:10
    global gurobi_env = Gurobi.Env()
#     global arcs
    setparams!(gurobi_env, Heuristics=0.0, Cuts = 0, OutputFlag = 0)
    println("\n### Instance ", ins)
    include("./TestInstances/Ins"*string(ins)*".jl")
    
    global numArcs = length(arcs[:,1])
    global Len = numArcs
    
    global P, gx
    global feasX, feasY, numX, numY, Q, P_set, numPaths
    feasX = EnumX(arcs,b,numArcs,d)
    Q = []
    P_set = [[1,8],[2,9],[3],[4,6],[5,7]]
#     feasX = []
    feasY = EnumX(arcs,b,numArcs,q.*10) #Using 10q as monitoring cost since all detection prob = 1
    numX = length(feasX)
    numY = length(feasY)
    numPaths = length(P_set)

    println("feasX = ", feasX)
    println("feasY = ", feasY)
    println("Path P = ", P_set)

    
    println("1")
    for x_i = 1:numX
        global Q, feasX, feasY, numY, P_set, numPaths
        q_I = zeros(numY,numPaths)
        for j = 1:numPaths # P_set
                P = P_set[j]
            for a in feasX[x_i]
                if a in P
                    q_I[:,j] = ones(numY)
                end
            end
            for i = 1:numY # in feasY
                M = feasY[i]
#                 println("Monitor ", M, " ; Path ", P)
                M_cap_P = findall(x->x in M, P)
                if length(M_cap_P) > 0
                    q_I[i,j] = 1.0
                end
            end
        end
        push!(Q, q_I)
    end
    for i = 1:numX
        println("x = ", feasX[i])
        println("Q = ", Q[i])
    end
#     break
    m = Model(() -> Gurobi.Optimizer(gurobi_env))
    
    @variable(m, x[1:numArcs], Bin)
    @variable(m, y[1:numX, 1:numY]>=0)
    @variable(m, v <= 10^6)
    
    @constraint(m, sum(d[i]*x[i] for i =1:numArcs)<= b)
    for p = 1:numPaths
        @constraint(m, v <= sum(Q[i][m,p]*y[i,m] for i = 1:numX for m = 1:numY))
    end
    for a = 1:numArcs
        global feasX, numX, feasY, numY
        for i = 1:numX
            if a in feasX[i]
                @constraint(m, sum(y[i,m] for m=1:numY) <= x[a])
            else
                @constraint(m, sum(y[i,m] for m=1:numY) <= 1-x[a])
            end
        end
    end
    @constraint(m, sum(y[i,m] for i = 1:numX for m = 1:numY) == 1)
    @objective(m, Max, v)
#     println(m)
    optimize!(m)
    current_Obj = JuMP.objective_value.(m)
    current_x = JuMP.value.(x)
#     X = findall(current_x.>0)
#     println(X)
#     M = findall(feasX .== X)[1]
#     as_ints(a::AbstractArray{CartesianIndex{L,2}}) where L = reshape(reinterpret(Int, a), (L, size(a)...))
    global myY = JuMP.value.(y)
    println("myY = ", myY)
    for i = 1:numX
        global myY
#         current_y = myY[i,:]
        println("Interdiction # = ", i, ": ", myY[i,:])
    end
#     current_y = as_ints(current_y)
#     global last_Obj = -1.0
#     global current_x
#     global current_y
#     global iter = 0
#     gx = -1
#     while iter < 5
# #     while current_Obj - gx > 0.001 #&& iter < 5
#         global arcs, q, d
#         global iter,current_Obj, current_x, current_y,gx
#         global numY, numPaths, numX, Q
#         iter +=1
#         optimize!(m)
#         current_Obj = JuMP.objective_value.(m)
#         current_x = JuMP.value.(x)
#         current_y = JuMP.value.(y)
#         println("\nIter = ", iter)
#         println("Obj = ", current_Obj)
#         println("x = ", current_x)
#         println("y = ", current_y)
#         interdict = findall(current_x.>0)
#         if (interdict in feasX) == false
#             push!(feasX, interdict)
#         end
        
#         i = findall(feasX, interdict)
#         println("Interdiction Sol #", i)
# #         dx = current_x + current_y
#         #Given current set of Monitoring decisions, solve for SP:
# #         P, gx = gx_bound(d, dx, arcs)
#         undetectionProb = ones(numArcs)
#         X = feasX[i]
#         for a in X
#             undetectionProb[a] -= 1
#         end
#         Y = findall(current_y[i,:].>0)
#         println("Y = ", Y)
#         for M in Y
#             for a in feasY[M]
#                 undetectionProb[i] -= current_y[i,M]
#             end
#     end
#             if undetectionProb[i] < 0
#                 undetectionProb = 0
#             end
#         end
#         println("undetectionP = ", undetectionProb)
#     break
#         P = longestPathLP(arcs, numArcs, costs)
# #         println("gx = ", gx)

#         @constraint(m, v <= sum(x[a] + q[a]*y[a] for a in findall(P.>0)))
# #         @printf("current_Obj - gx = %1.4f", current_Obj - gx)
#     end
#     println("\n=====OPTIMAL=====")
#     outfile = "C:/Users/din/Documents/GitHub/UncertainTarget/TestInstances/Output_Ins"*string(ins)*".jl"
#     io = open(outfile, "w") #do io
#     println(io, "Obj = ", current_Obj)
#     println(io, "x = ", findall(current_x.>0))
#     println(io, "y = ", findall(current_y.>10^(-6))," ; ", current_y[current_y.>10^(-6)])
#     close(io)
    println("Obj = ", current_Obj)
    println("x = ", findall(current_x.>0))
#     println("y = ", findall(current_y.>10^(-6))," ; ", current_y[current_y.>10^(-6)])

# end


using JuMP
using Gurobi
using Combinatorics

global arcs
global d
global q

global b = 5


global gurobi_env = Gurobi.Env()
setparams!(gurobi_env, Heuristics=0.0, Cuts = 0, OutputFlag = 0)

global origin = 1
global destination = 6
include("./Ins3.jl")
global numArcs = length(arcs[:,1])
global Len = numArcs
include("C:/Users/din/Documents/GitHub/IntCVaR/functionGbound.jl")

global X 

function EnumX(arcs,b,numArcs,d)
    X = findall(d.<=b)
    Len = length(X)
    X = combinations(X) |> collect
#     redunX = []
#     global 
    feasX = []
    println(X)
    #Eliminate combinations exceeding budget
    for i in X
        cost = sum(d[x] for x in i)
#         println("X = ",i,": ",cost)
        if cost <= b
            push!(feasX, i)
        end
    end
    println("feasible Budget = ", feasX)
#     global 
    subset = true
    
    #Eliminate non-maximally packed x-solutions
    while subset == true
        subset = false
        for i in feasX
            for j in feasX
                if i != j
                    if issubset(i, j) == true
                        subset = true
                        if length(i) < length(j)
                            filter!(x->x != i, feasX)
                        else
                            filter!(x->x != j, feasX)
                        end
                    end
                end
            end
        end    
    end
    println("feasX = ", feasX)
    return feasX
end

function solve2ndStage(x_ind, arcs, numArcs, q, d, origin, destination)
    global gurobi_env
    x = zeros(numArcs)
    for a in x_ind
       x[a] = 1 
    end
#     println("x = ", x)
    m = Model(() -> Gurobi.Optimizer(gurobi_env))

    @variable(m, y[1:numArcs]>=0)
    @variable(m, v >= 0)
    
    #Flow Constraints:
#     println("Flow")
    for i = 1:(destination-1)
        inflow = findall(arcs[:,2].==i)
        outflow = findall(arcs[:,1].==i)
#         println("Fine")
        if i == origin
            con = @constraint(m, sum(y[a] for a in outflow) == 1)
            println(con)
        else
            con = @constraint(m, sum(y[a] for a in outflow) - sum(y[a] for a in inflow) == 0)
            println(con)
        end
    end
#     println("a")
    #Arc Capacity
    for a = 1:numArcs
        if q[a] > 0
            con = @constraint(m, y[a] <= v*1/q[a])
            println(con)
        else
            con = @constraint(m, y[a] <= v*(10^6))
            println(con)
        end
        con = @constraint(m, y[a] <= (10^6)*(1-x[a]))  
        println(con)
    end
#     @constraint(m, sum(d[i]*x[i] for i =1:numArcs)<= b)
#     println(findall(P.>0))
#     @constraint(m, v <= sum(x[a] + q[a]*y[a] for a in findall(P.>0)))
    @constraint(m, sum(y[a] for a =1:numArcs) == 1)
    @objective(m, Min, v)
    
    optimize!(m)
    Obj_Val = JuMP.objective_value.(m)
    return Obj_Val
end

X = EnumX(arcs,b,numArcs,d)
println("X = ", X)
Obj =  solve2ndStage(X[1], arcs, numArcs, q, d, origin, destination)
println(Obj)
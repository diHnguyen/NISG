using Gurobi
using JuMP
using LightGraphs
using Printf

# global arcs
# global d
# global q

# global b = 5

include("C:/Users/din/Documents/GitHub/IntCVaR/functionGbound.jl")
# global origin = 1
# global destination = 6
for ins = 1:10
    global gurobi_env = Gurobi.Env()
#     global arcs
    setparams!(gurobi_env, Heuristics=0.0, Cuts = 0, OutputFlag = 0)
    println("\n### Instance ", ins)
    include("./TestInstances/Ins"*string(ins)*".jl")
    global numArcs = length(arcs[:,1])
    global Len = numArcs
    
    global P, gx
    P, gx = gx_bound(d, d, arcs)
#     println(P)
#     println(arcs[P.>0,:])
    # for a in path
    #    print(arcs[a,:], " ") 
    # end
    m = Model(() -> Gurobi.Optimizer(gurobi_env))

    @variable(m, x[1:numArcs], Bin)
    @variable(m, y[1:numArcs]>=0)
    @variable(m, v <= 10^6)

    @constraint(m, sum(d[i]*x[i] for i =1:numArcs)<= b)
#     println(findall(P.>0))
    @constraint(m, v <= sum(x[a] + q[a]*y[a] for a in findall(P.>0)))
    @constraint(m, sum(y[a] for a =1:numArcs) == 1)
    @objective(m, Max, v)
    global last_Obj = -1.0
    global current_Obj = 0.0
    # println(typeof(last_Obj))
    global current_x
    global current_y
    global iter = 0
    gx = -1
    while current_Obj - gx > 0.001 #&& iter < 5
        global arcs, q, d
        global iter,current_Obj, current_x, current_y,gx
        iter +=1
#         println("\nIter ", iter)
        optimize!(m)
        current_Obj = JuMP.objective_value.(m)
#         println("current_Obj = ", current_Obj)
        current_x = JuMP.value.(x)
        current_y = JuMP.value.(y)
#         println("x = ", findall(current_x.>0))
#         println("y = ", findall(current_y.>0))

        dx = current_x + q.*current_y
        P, gx = gx_bound(d, dx, arcs)
#         println("gx = ", gx)

        @constraint(m, v <= sum(x[a] + q[a]*y[a] for a in findall(P.>0)))
#         @printf("current_Obj - gx = %1.4f", current_Obj - gx)
    end
#     println("\n=====OPTIMAL=====")
    outfile = "C:/Users/din/Documents/GitHub/UncertainTarget/TestInstances/Output_Ins"*string(ins)*".jl"
    io = open(outfile, "w") #do io
    println(io, "Obj = ", current_Obj)
    println(io, "x = ", findall(current_x.>0))
    println(io, "y = ", findall(current_y.>10^(-6))," ; ", current_y[current_y.>10^(-6)])
    close(io)
    println("Obj = ", current_Obj)
    println("x = ", findall(current_x.>0))
    println("y = ", findall(current_y.>10^(-6))," ; ", current_y[current_y.>10^(-6)])

end


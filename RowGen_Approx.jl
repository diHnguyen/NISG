using Gurobi
using JuMP
using LightGraphs

global gurobi_env = Gurobi.Env()
outfile = "C:/Users/din/Documents/GitHub/UncertainTarget/Log_Test_Output.csv"
    io = open(outfile, "w")
println(io,"Instance \tPath \tExact \tApprox")
global p_exact = []
global p_approx = []
for ins = 1:10
#     ins = 3
    include("./TestInstances/Ins"*string(ins)*".jl")
    include("./functionEnumFeasX.jl")
    global numArcs = length(arcs[:,1])
    global Len = numArcs
    global p_exact, p_approx

    global feasY = EnumX(arcs, b, numArcs, d)
    println(feasY)
    numY = length(feasY)
    for count = 1:5
        global q, P, feasY
        println("\nIter ", count)
        y = rand(numY)
        temp = sum(y[i] for i = 1:numY)
        y = y./temp
        # global y = [0.198021, 0.220998, 0.247267, 0.333714]
        
        println("y = ", y)

    #     P = [[1,7],[1,5,10], [1,6,9], [1,4,10], [3,9]]

        #Calculating the exact probability that the evader is not detected on path P:
        println("Exact Pr of being undetected: ")
        p_exact = [] # zeros(length(P))
        for path in P
            p = 0
            for m = 1:numY
                p_y = 1
        #         M_in_P = false
                for a in path
                    if a in feasY[m]
                        p_y = p_y*(1-q[a])
        #                 M_in_P = true
                    end
                end
                p = p+ p_y*y[m]
                
            end
            push!(p_exact, p)
            println("Path ", path, ": ", p)
        end

        #Approx
        println("Approx. Pr of being undetected")
        p_approx =[]
        for path in P
            p = 0
            for m = 1:numY
                p_y = 0
        #         M_in_P = false
                for a in path
                    if a in feasY[m]
                        p_y = p_y + log(1-q[a])
        #                 M_in_P = true
                    end
                end
                p = p + (p_y+1)*y[m]
                
            end
            push!(p_approx, p)
            println("Path ", path, ": ", p)
            
        end
        println("p_exact ", p_exact)
        println("p_approx ", p_approx)
        for j = 1:length(P)
            println(io,ins,";", count , ";", string(P[j]), ";", p_exact[j],";", p_approx[j])
        end
    end
     #do io
#     println(io, "Obj = ", current_Obj)
#     println(io, "x = ", findall(current_x.>0))
#     println(io, "y = ", findall(current_y.>10^(-6))," ; ", current_y[current_y.>10^(-6)])
    
#     println("Obj = ", current_Obj)
#     println("x = ", findall(current_x.>0))
#     println("y = ", findall(current_y.>10^(-6))," ; ", current_y[current_y.>10^(-6)])

end
close(io)
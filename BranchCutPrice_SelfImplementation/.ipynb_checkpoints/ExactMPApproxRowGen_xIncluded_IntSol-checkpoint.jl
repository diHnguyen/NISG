#Updated on 6/28
#All monitoring decisions are included in the MP
#The coefficients of each arc in the MP is exact
#The path with the least detection probability = most probability of being undetected is approximate


#Notes: Initialization of P_set and M_Set is important to get to the optimal solution.
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

include("./TestInstances/Ins"*string(ins)*".jl")
include("./functionEnumFeasX.jl")
include("./functionFindQ.jl")
include("./functionShortestPathBellmanFord.jl")
include("./findExactUndetectionPr.jl")
include("./functionGenMonitoring.jl")

include("./functionSolveNode.jl")

global objVal_inc = 2
global x_inc = []
global y_inc = []

global numArcs = length(arcs[:,1])
global M_set = [[1,7]] 


# M_set = EnumX(arcs, b_y, numArcs, d_y)


global numY = length(M_set) #
# println("M_set = ", M_set)
global P_set = [[1,8]]
# global P_set = P
# println("P_set =", P_set )
global numPaths = length(P_set)#1
global Q = [] #length(Q) = numPaths
global R = 1 
# coef_q, Q = findQgivenP(M_set, numY, Q, P_set[1])
for path in P_set
    global Q
    coef_q, Q = findQgivenP(M_set, numY, Q, path)
end
# println("Q = ", Q)

global allNodes = [1]
global nodeStat = [1]
global nodePred = [0]
global nodeObj = [0.0]
global nodeX = [[]]
global nodeA = [[]]

global numNodes = 1
# for myCount = 1:4
global openNodes = findall(nodeStat.==1)
while isempty(openNodes) == false #&& numNodes <7
    global objVal_inc, x_inc, y_inc
    global numArcs, M_set, numY, P_set, numPaths, Q, R 
    global nodes, nodeStat, nodePred, nodeObj, nodeX, numNodes, openNodes
    
    #Select highest node num in the set of all open nodes
    node = findall(nodeStat.==1)[end]
    println("\nSOLVING NODE ", node)
    println("P_set = ",P_set)
    println("M_set = ", M_set)
    println("Q = ", Q)
    
    #Arrive at a node (var node) 
    obj, x, y = solveNode(nodeA[node], nodeX[node])
    #Set of arcs with non-binary x-variables
    A_split = findall((x.>0) .& (x.<1))
    
    #If incumbent solution > objective val at node
    #This condition needs to be revised
    if objVal_inc > obj
#         println("nodeA = ", nodeA)
#         println("nodeX = ", nodeX)
        #If we still have non-binary x variables
        if length(A_split) > 0 
            println("Not integer")
            if length(allNodes) < node
            else
                nodeObj[node] = obj
    #             nodeX[node] = x
            end
    #         println("nodeA = ", nodeA)
            
            #Get the set of arcs that have been branched on up to the current node
            #AKA the set of x-constraints of either x=1 or x=0
            newA = copy(nodeA[node])
    #         println("newA = ", newA)
    #         x_split = argmax(x)
            arc_split = 0
            x_split = 0
            
            #Select arc arc_split that has largest non-binary x-variable 
            for a in A_split
                if x[a] > x_split
                    arc_split = a
                    x_split = x[a]
                end
            end
            
            #Update the set of arcs branched on for the child-nodes
            push!(newA, arc_split)
    #         println("newA = ", newA)
    #         println("nodeA = ", nodeA)
            
            #Update the corresponding x-variables for each child-node
            currentX0 = copy(nodeX[node])
            currentX1 = copy(nodeX[node])
    #         nodeStat[node] = 0
            push!(currentX0, 0)
            push!(currentX1, 1)
            push!(nodeX, currentX0)
            push!(nodeX, currentX1)
            
            #For each of the new node
            for n = 1:2
                numNodes = numNodes + 1 #Increase the total number of nodes in the tree
                push!(allNodes, numNodes) #Add node index to list allNodes
                push!(nodeStat,1) #Set child-node status to open
                push!(nodePred,node) #Set child-node's parent to current node
                push!(nodeObj,0) #Set node objective function value to 0 (doesn't matter)
                push!(nodeA, newA) #Add the corresponding set of have-been-branched-on arcs for the new child-node to list
            end
        else
            println("Integer - Update incumbent sol?")
            
            if objVal_inc > obj
                objVal_inc = obj
                x_inc = x
                y_inc = y
                println("New inc obj = ", objVal_inc)
                println("New inc x = ", x_inc)
                println("New inc y = ", y_inc)
            end

        end
    else
        println("Node not better than incumbent. Don't split node.")
    end
    
    #Set current (parent) node status to closed
    nodeStat[node] = 0
    
    #Find next node to get to
    #Can remove this actually (since we go back to beginning of lop)
    openNodes = findall(nodeStat.==1)
#     println("Open Nodes ", openNodes)
#     println(nodeA[openNodes])
#     println(nodeX[openNodes])
    
end
println("\n--------------------------")
println("Final obj = ", objVal_inc)
println("Final x = ", x_inc)
println("Final y = ", y_inc)
println("nodePred = ", nodePred)
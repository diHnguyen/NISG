using GraphPlot
using JuMP
using LightGraphs, SimpleWeightedGraphs
using Combinatorics

ins = 10

include("./TestInstances/Ins"*string(ins)*".jl")
# global interdict = [7, 8]
global interdict = []
global N = maximum(arcs[:,2])
global numArcs = length(arcs[:,2])

include("C:/Users/din/Documents/GitHub/UncertainTarget/functionEnumFeasX.jl")
global X = EnumX(arcs,b,numArcs,d)
println(X)

# global newArcs = arcs[setdiff(1:end, interdict),:]
# println(newArcs)

# nodelabel = collect(1:N)
# global interdict = [3]
outfile = "C:/Users/din/Documents/GitHub/UncertainTarget/TestInstances/Output_Ins"*string(ins)*".jl"
io = open(outfile, "a") #do io
println(io, "#Enumerate Cuts")
# println(io, "Set of maximally packed x-solution = ", X)
for interdict in X
#     global g
    global N
    g = SimpleWeightedDiGraph(N)
    for i = 1:numArcs #length(newArcs[:,1])
        if (i in interdict) == false
            add_edge!(g,arcs[i,1],arcs[i,2], 1) #d[i]) 
        end
    end
    A_C = mincut(g,weights(g))
    println(io, "x = ", interdict, " : Size of min Cut on Residual Network = ", A_C[2])
end
close(io)


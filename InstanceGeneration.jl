using LightGraphs
using GraphPlot
# using Compose
using FileIO
using CSV
# g = SimpleDiGraph(2);
# add_edge!(g, 1, 2);
# add_edge!(g, 2, 1);
# # gplot(g)
# draw(PNG("mygraph.png", 8cm, 8cm), gplot(g))


#Specify number of arcs and nodes

global N = 6

#Probability an arc (i,j) is present/Density
global p = 0.4

function makeArcs(N,p)
    global arcs
    for i = 1:N-1
        for j = i:N
            a = rand()
            if a <= p && i != j
                arcs = vcat(arcs, [i j])
            end
        end
    end
    return arcs
end
function testGraph()
    #Source 1:
    global N, arcs
    use = true
    outgoing = arcs[:,1]
    incoming = arcs[:,2]
    toSource = []
    toSink = []
    if length(outgoing[outgoing.==1]) >= 1
#         println("Node 1", 1)
        if length(incoming[incoming.==N]) >= 1
            for i = 2:N-1
                if length(outgoing[outgoing.==i]) == 0
                    use = false
                    push!(toSink, i)
                end
                if length(incoming[incoming.==i]) == 0
                    use = false
                    push!(toSource, i)
                end
            end
        end
    end
    return use, toSource, toSink
end
function addAdditionalArcs(toSource,toSink)
    global N, arcs
    for j in toSource
        arcs = vcat(arcs, [1 j])
    end
    for i in toSink
       arcs = vcat(arcs, [i N]) 
    end
    return arcs
end

for ins = 1:10
    global arcs = Array{Int64}(undef, (0,2))
    arcs = makeArcs(N,p)
    global numArcs
    numArcs = length(arcs[:,1])
    println(arcs)
    println("numArcs = ", numArcs)

    # arcs = addAdditionalArcs(toSource, toSink)
    # println(arcs)



    use, toSource, toSink = testGraph()
    println("use = ", false)
    println("toSource = ", toSource)
    println("toSink = ", toSink)

    arcs = addAdditionalArcs(toSource,toSink)

    println(arcs)
    d = rand(0:10, length(arcs[:,1]))
    q = rand(0:10, length(arcs[:,1]))./10

    outfile = "C:/Users/din/Documents/GitHub/UncertainTarget/TestInstances/Ins"*string(ins)*".jl"
    io = open(outfile, "w") #do io
    #     write(io, "arcs = ", arcs )
    #     write(io, "d = ", d )
    #     write(io, "q = ", q )
    # end
    println(io, "arcs = ", arcs )
    println(io, "d = ", d )
    println(io, "q = ", q )
    close(io)
end

# x = length(arcs[:,1])
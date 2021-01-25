arcs = [1 2; 1 5; 2 3; 2 4; 3 6; 4 7; 5 4; 5 8; 6 10; 7 10; 8 7; 8 9; 9 10]
Len = length(arcs[:,1])
N1 = [2,5]

#c = ones(length(arcs))
c = [0.0, 0.0, 17.0, 7.0, 1.0, 10.0, 12.0, 12.0, 0.0, 0.0, 8.0, 2.0, 0.0]
b = 2 
T = [6, 7, 9]
v = [10,20,30]
d = ones(Len).*100
scenNum = length(T)
p_init = ones(scenNum)./scenNum
origin= 1
destination = 10

#Get stage 2 starting nodes:
N2 = Int64[]
for n in N1
    outArcs = findall(arcs[:,1] .== n)
    for i in outArcs
        if (arcs[i, 2] in N2) == false
            push!(N2, arcs[i, 2])
        end
    end
end
sort!(N2)

#Get A_star
A_star = findall(arcs[:,1] .== origin)
A_star = vcat(A_star, findall(arcs[:,2] .== destination))

arcs = [1 2; 1 5; 2 3; 2 4; 3 6; 4 7; 5 4; 5 9; 6 7; 6 10; 7 8; 7 10; 8 10; 9 7; 9 8; 9 10]
A = length(arcs[:,1])
c = ones(length(arcs))
b = 2 
T = 4
p = ones(T)./T
scen = [6, 7, 8, 9]
source = 1
sink = 10
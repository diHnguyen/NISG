# A Two-Stage Network Interdiction-Monitoring Game - #NISG
# Di H. Nguyen, Yongjia Song, and J. Cole Smith

Link to paper: TBD

Abstract: We study a network interdiction problem involving two agents: a defender and an evader. 
The evader seeks to traverse a path from a source node to a sink node in a directed network without being detected. 
The game takes place in two stages. In the first stage, the defender removes a set of arcs in the network. 
In the second stage, the defender and evader play a simultaneous game. 
The defender monitors a set of arcs, thus increasing the probability that the evader will be detected on that arc (if the evader uses the arc). 
The evader selects a source-sink path. 
Because the second stage is played simultaneously, both agents use mixed-strategy solutions. 
We approach the solution of the second-stage problem by proposing a constraint-and-column generation algorithm. 
We show that both the constraint-generation and column-generation problems are NP-hard. 
Accordingly, we prescribe approximate versions of these problems that can be solved as linear programs. 
Our algorithm relies on solving these linear programming problems until an exact solution of the constraint-generation and column-generation problems is required to prove optimality. 
Then, to link the first- and second-stage problems, we model the original problem using an epigraph reformulation, which we solve using a Benders-decomposition based approach. 
The efficacy of our approach is demonstrated on a set of randomly generated test instances.


In this repository you will find the following: Source Code + Test Instances 

(Some of which share the same "base" instance with some parameter(s) modified. Please refer to the Computational Section of the paper to see how these are generated.)


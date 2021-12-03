Obj = 0.4
x = [2, 3, 8]
y = [1] ; [1.0]
#Verify
Set of maximally packed x-solution = Any[[2, 10], [2, 3, 5], [2, 3, 8], [2, 5, 8]]
x = [2, 10] ; v** = 0.1984251968503937 ; minCut = [1, 3, 9] ; [0.496063, 0.283465, 0.220472]
x = [2, 3, 5] ; v** = 0.27692307692307694 ; minCut = [1, 9] ; [0.692308, 0.307692]
x = [2, 3, 8] ; v** = 0.4 ; minCut = [1] ; [1.0]
x = [2, 5, 8] ; v** = 0.2545454545454545 ; minCut = [1, 3] ; [0.636364, 0.363636]
#Enumerate Cuts
x = [2, 10] : Size of min Cut on Residual Network = 3.0
x = [2, 3, 5] : Size of min Cut on Residual Network = 2.0
x = [2, 3, 8] : Size of min Cut on Residual Network = 2.0
x = [2, 5, 8] : Size of min Cut on Residual Network = 3.0

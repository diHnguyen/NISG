Obj = 0.9000000000000004
x = [1, 3, 6, 7]
y = [9] ; [1.0]

#Verify
Set of maximally packed x-solution = Any[[2, 6, 7], [6, 7, 8], [6, 7, 9], [1, 3, 5, 6, 7]]
x = [2, 6, 7] ; v** = 0.3428571428571428 ; minCut = [1, 3] ; [0.428571, 0.571429]
x = [6, 7, 8] ; v** = 0.36 ; minCut = [3, 9] ; [0.6, 0.4]
x = [6, 7, 9] ; v** = 0.3428571428571428 ; minCut = [1, 3] ; [0.428571, 0.571429]
x = [1, 3, 5, 6, 7] ; v** = 0.8999999999999999 ; minCut = [9] ; [1.0]
#Enumerate Cuts
x = [2, 6, 7] : Size of min Cut on Residual Network = 1.0
x = [6, 7, 8] : Size of min Cut on Residual Network = 1.0
x = [6, 7, 9] : Size of min Cut on Residual Network = 3.0
x = [1, 3, 5, 6, 7] : Size of min Cut on Residual Network = 1.0

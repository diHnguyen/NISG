function tightenConstraint(lambda, pi_,lambda_pos, P_set,numPaths, numArcs)
#     println("Inside function tightening constraint -- Info:")
#     println("\tlambda = ", lambda)
#     println("\tpi_ = ", pi_)
#     println("\tlambda_pos = ", lambda_pos)
#     println("\tP_set = ", P_set)
#     println("\tnumPaths = ", numPaths)
#     println("\tnumArcs = ", numArcs)
    x_coefs = zeros(numArcs)
    #theta-constraint BEFORE tightening
    #con_ = @build_constraint(theta >= -sum(lambda[i]*sum(x[a] for a in P_set[i]) for i=1:numPaths) + pi_)
    
    X_set = reduce(vcat, P_set[lambda_pos])
#     println("\tX_set = ", X_set)
    unique!(X_set)
#     println("\tUnique X_set = ", X_set)
#     x_coefs[X_set] = pi_
    for i in lambda_pos
        myPath = P_set[i]
        x_coefs[myPath] = x_coefs[myPath].+lambda[i]
#         println("\tx_coefs = ", x_coefs)
    end
    for a in findall(x_coefs.>pi_)
        x_coefs[a] = pi_
    end
#     println("\tx_coefs_adjusted = ", x_coefs)
    return x_coefs, X_set
end
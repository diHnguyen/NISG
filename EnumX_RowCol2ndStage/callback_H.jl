cb_calls = Cint[]
    function my_callback_function(cb_data, cb_where::Cint)
#     function my_callback_function(cb_data::Gurobi.CallbackData, cb_where::Int32)
        push!(cb_calls, cb_where)
        global lambda, pi_, P_set, M_set, R
        println("HERE", cb_where, " vs GRB_CB_MIPSOL ", GRB_CB_MIPSOL)
        println(lambda," ", pi_)
        if cb_where == GRB_CB_MIPSOL #Gurobi.CB_MIPSOL #Gurobi.CB_MIPSOL
            Gurobi.load_callback_variable_primal(cb_data, cb_where)
            x_val = callback_value.(Ref(cb_data), x)
            u_val = callback_value(cb_data, u)
            println("x_val = ", x_val)
            println("u_val = ", u_val)
            TOL = 1e-6
            numPaths = length(P_set)
            numY = length(M_set)
            dual_Obj = sum(lambda[i]*(-R)*sum(x_val[a] for a in P_set[i]) for i=1:numPaths)+pi_
            if u_val+ TOL < dual_Obj#y_val - x_val > 1 + TOL
                con = @build_constraint(u>= sum(lambda[i]*(-R)*sum(x[a] for a in P_set[i]) for i=1:numPaths)+pi_)
                println("Adding ", con)
                MOI.submit(MP, MOI.LazyConstraint(cb_data),con)
            end
        end
    end
    MOI.set(MP, MOI.RawParameter("LazyConstraints"), 1)
    MOI.set(MP, Gurobi.CallbackFunction(), my_callback_function)
using JuMP, Gurobi
using LightGraphs #, LightGraphsFlow
using GraphPlot
const gurobi_env = Gurobi.Env()
model = direct_model(Gurobi.Optimizer())
m = Model(() -> Gurobi.Optimizer(gurobi_env))
    set_optimizer_attribute(m, "OutputFlag", 0)

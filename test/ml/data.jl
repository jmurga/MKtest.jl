# Shared synthetic data of the dfe_ml tests (n = 40): Pn_j = 2000/j, Ps_j = 3000/j.
# Included by engine.jl and extension.jl, which may run in separate processes.
const ML_JS = collect(1.0:39.0)
const ML_SFS = hcat(ML_JS, round.(2000.0 ./ ML_JS), round.(3000.0 ./ ML_JS))
const ML_DIV = [900.0, 1200.0, 1.0e6, 3.0e5]

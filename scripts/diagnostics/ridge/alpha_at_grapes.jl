using MKtest, CSV, DataFrames
tab = CSV.read(joinpath(@__DIR__, "..", "..", "..", "test/fixtures/grapes_parity_tgp.csv"), DataFrame)
sfs, div, _ = MKtest.parse_sfs(joinpath(@__DIR__, "..", "..", "..", "data/MKdata_tgp_may2023.txt"); sample_size = 661)
div = Float64.(vec(div))
rb = MKtest.RombergIntegrator()
models = Dict("GammaZero" => MKtest.GammaZero(), "GammaExpo" => MKtest.GammaExpo(), "GammaGamma" => MKtest.GammaGamma())
function alpha_at(m, s, θ)
    obj = MKtest._objective_for(m, s, div; integrator = rb)
    lay = obj.layout
    x = [MKtest._to_opt(θ[i], lay.log_scale[i]) for i in lay.opt_idx]
    MKtest._set_params!(obj, x); ll, th = MKtest._loglik!(obj)
    f = MKtest._FitOutcome(obj, x, ll, th, true, length(x))
    c = MKtest.FitContext(f; nearly_neutral = 5.0, alpha_ci_threshold = 0.0)
    return ll, MKtest.summarize(MKtest.DFEMAlpha(), c).alpha, (c.exp_div_N1, c.exp_div_N2, c.exp_div_N3, c.totN)
end
for row in eachrow(sort(tab, [:model, :n_ind]))
    haskey(models, row.model) || continue
    m = models[row.model]; s = MKtest.project(sfs, row.n_ind)
    lay = MKtest.build_layout(m, MKtest.PoissonKernel(), true)
    θg = Float64[row["p_" * MKtest.grapes_param_name(nm)] for nm in lay.names]
    llg, ag, _ = alpha_at(m, s, θg)
    jf = MKtest.dfe_ml(s, div, m; integrator = rb)
    θj = [MKtest.param_values(jf, nm)[1] for nm in lay.names]
    println(rpad("n=$(row.n_ind) $(row.model)", 16),
        " GRAPES: lnL ", round(row.lnL; digits = 3), " alpha ", row.alpha,
        " | Julia at GRAPES params: lnL ", round(llg; digits = 3), " alpha ", round(ag; digits = 4),
        " | Julia optimum: lnL ", round(jf.inference.loglik[1]; digits = 3), " alpha ", round(jf.inference.alpha[1]; digits = 4))
    println("    GRAPES params ", Dict(zip(lay.names, round.(θg; sigdigits = 4))))
    println("    Julia  params ", Dict(zip(lay.names, round.(θj; sigdigits = 4))))
end

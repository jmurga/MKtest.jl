# GRAPES parity of dfe_ml (Task 1, D10). Oracle: the GRAPES binary, run by
# test/fixtures/make_grapes_parity_fixture.jl on the tgp SFS projected to 10 and 20
# individuals; its output is test/fixtures/grapes_parity_tgp.csv, so this test needs
# no binary. Every check uses RombergIntegrator(), GRAPES' own integration rule.
#
# Assertions, fixed in the Phase 0 report (section 5) before the fixture existed:
#   objective parity (all models): the Julia log-likelihood AT the binary's
#     parameters is within 1e-3 of the binary's lnL;
#   fit parity: GammaZero lnL and alpha within 1e-3 of the binary; the ridge and
#     stalling models (GammaExpo, GammaGamma, DisplGamma, ScaledBeta) reach a lnL >=
#     the binary's - 1e-3, and their alpha is reported, not asserted.
# A binary row whose lnL equals the Julia lnL at one of GRAPES' start points
# (within 1e-3) is reported as "optimizer stalled".

using CSV
using DataFrames

const PARITY_CSV = joinpath(@__DIR__, "..", "fixtures", "grapes_parity_tgp.csv")
const PARITY_TOL = 1e-3

_model_of(name) = Dict("GammaZero" => MKtest.GammaZero(), "GammaExpo" => MKtest.GammaExpo(),
                       "GammaGamma" => MKtest.GammaGamma(), "DisplGamma" => MKtest.DisplGamma(),
                       "ScaledBeta" => MKtest.ScaledBeta())[name]

@testset "GRAPES binary parity (tgp, n = 10 and 20 individuals)" begin
    @test isfile(PARITY_CSV)
    if isfile(PARITY_CSV)
        tab = CSV.read(PARITY_CSV, DataFrame)
        sfs, div, _ = MKtest.parse_sfs(joinpath(@__DIR__, "..", "..", "data", "MKdata_tgp_may2023.txt");
                                       sample_size = 661)
        div = Float64.(vec(div))
        cache = MKtest.DFECache()
        rb = MKtest.RombergIntegrator()
        report = DataFrame(n_ind = Int[], model = String[], lnL_grapes = Float64[], lnL_julia_at_grapes = Float64[],
                           lnL_julia_fit = Float64[], alpha_grapes = Float64[], alpha_julia = Float64[],
                           stalled = Bool[])
        for row in eachrow(tab)
            m = _model_of(row.model)
            s = MKtest.project(sfs, row.n_ind)
            lay = MKtest.build_layout(m, MKtest.PoissonKernel(), true)
            # the fixture holds GRAPES' own names (p_negGmean, ...); map ours to them
            cols = ["p_" * MKtest.grapes_param_name(nm) for nm in lay.names]
            missing_cols = [c for c in cols if !(c in names(tab)) || ismissing(row[c])]
            @test isempty(missing_cols)
            isempty(missing_cols) || (@warn "GRAPES output lacks $missing_cols for $(row.model)"; continue)
            θ = Float64[row[c] for c in cols]

            # objective parity
            ll_at = MKtest._loglik_at(m, s, div, θ; integrator = rb, cache = cache)
            @test abs(ll_at - row.lnL) < PARITY_TOL

            # fit parity
            jf = MKtest.dfe_ml(s, div, m; integrator = rb, cache = cache)
            if row.model == "GammaZero"
                @test abs(jf.inference.loglik[1] - row.lnL) < PARITY_TOL
                @test abs(jf.inference.alpha[1] - row.alpha) < PARITY_TOL
            else
                @test jf.inference.loglik[1] >= row.lnL - PARITY_TOL
            end

            # did the binary stay at one of its start points?
            starts = [vcat(sp, oer) for sp in MKtest.dfe_start_points(m) for oer in (0.0, 0.01)]
            ll_starts = [MKtest._loglik_at(m, s, div, st; integrator = rb, cache = cache) for st in starts]
            stalled = any(abs.(ll_starts .- row.lnL) .< PARITY_TOL)

            push!(report, (row.n_ind, row.model, row.lnL, ll_at, jf.inference.loglik[1], row.alpha, jf.inference.alpha[1],
                           stalled))
        end
        show(stdout, report; allrows = true, allcols = true)
        println()
        any(report.stalled) && @warn "GRAPES optimizer stalled at a start point for: " *
            join(["n=$(r.n_ind) $(r.model)" for r in eachrow(report) if r.stalled], ", ")
    end
end

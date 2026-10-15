# Refactor guard: the moved engine against dev's at 9d752e40, on the inputs dev
# itself used (test/fixtures/dev_reference/, from make_dev_reference.jl). A frozen
# record of the code the engine came from, not a correctness oracle.
#
# Tolerances, fixed before the reference existed:
#   objective at fixed parameters: relative 1e-10. The new engine sums the expected
#     SFS in a plain loop where dev used BLAS `mul!`, which reorders the sum (about
#     one ulp per entry); nothing else in the arithmetic changed.
#   fits (as first written): GammaZero lnL and alpha within 1e-6; other models lnL >=
#     dev's - 1e-6.
#
# CHANGED 2026-10-15, AFTER seeing the first run (user decision): 80/80 objective checks
# passed at 1e-10, and 5 fit checks failed. Two were ~1e-6 (inside Optim's own stop rule,
# f_reltol 1e-8 x |lnL| ~ 3e-6). Three were DisplGamma, 0.27-4.8 lnL BELOW dev. The new
# engine under dev's own Optim 1.13.3 moved the same cells by +5.01 / +0.03 / -0.28, so
# DisplGamma's optimum depends on ulp-level path differences in both engines; it is not
# a likelihood change. Now: the objective stays the guard (1e-10); fits are asserted for
# GammaZero only, at Optim's tolerance (1e-5), and every other fit is reported.

using DelimitedFiles

const DEVREF = joinpath(@__DIR__, "..", "fixtures", "dev_reference")

function _parse_params(s::AbstractString)
    d = Dict{Symbol,Float64}()
    for kv in split(s, ";"; keepempty = false)
        k, v = split(kv, "=")
        d[MKtest.from_legacy_name(k)] = parse(Float64, v)   # the record uses dev's old names
    end
    return d
end

_devref_model(name) = Dict("GammaZero" => MKtest.GammaZero(), "GammaExpo" => MKtest.GammaExpo(),
                           "DisplGamma" => MKtest.DisplGamma(),
                           "WithLethal{GammaZero}" => MKtest.WithLethal(MKtest.GammaZero()))[name]
_devref_integ(name) = name == "RombergIntegrator" ? MKtest.RombergIntegrator() : MKtest.ExactSojourn()

@testset "refactor guard: dev dfe_ml at 9d752e40" begin
    ref = joinpath(DEVREF, "dev_reference.csv")
    @test isfile(ref)
    if isfile(ref)
        rows = filter(l -> !startswith(l, "#"), readlines(ref))
        hdr = split(rows[1], '\t')
        cache = MKtest.DFECache()
        fits = Dict{Tuple{String,String,String},Any}()
        for line in rows[2:end]
            f = Dict(zip(hdr, split(line, '\t')))
            s = readdlm(joinpath(DEVREF, "inputs_$(f["input"]).tsv"), '\t', Float64)
            d = vec(readdlm(joinpath(DEVREF, "inputs_div_$(f["input"]).tsv"), '\t', Float64))
            m = _devref_model(f["model"])
            integ = _devref_integ(f["integrator"])
            p = _parse_params(f["params"])
            lay = MKtest.build_layout(m, MKtest.PoissonKernel(), true)
            θ = [p[nm] for nm in lay.names]
            ll_dev = parse(Float64, f["loglik"])
            if f["kind"] == "objective"
                ll = MKtest._loglik_at(m, s, d, θ; integrator = integ, cache = cache)
                @test isapprox(ll, ll_dev; rtol = 1e-10)
            else
                key = (f["input"], f["model"], f["integrator"])
                r = get!(fits, key) do
                    MKtest.dfe_ml(s, d, m; integrator = integ, cache = cache)
                end
                if f["model"] == "GammaZero"
                    @test abs(r.inference.loglik[1] - ll_dev) < 1e-5
                    @test abs(r.inference.alpha[1] - parse(Float64, f["alpha"])) < 1e-5
                end
                @info "dev reference $key" lnL_dev = ll_dev lnL_new = r.inference.loglik[1] alpha_dev = f["alpha"] alpha_new = r.inference.alpha[1]
            end
        end
    end
end

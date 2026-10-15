# D11: does dfe_ml give approximately the same answer as the sample size grows?
#
# A CONSISTENCY CHECK, NOT PARITY AND NOT EVIDENCE OF ACCURACY. There is no binary
# beyond n = 20 individuals and no known truth for the tgp data; projecting the SFS
# leaves the expected DFE unchanged, so the estimates should move little and in no
# systematic direction.
#
# Input: data/MKdata_tgp_may2023.txt, parse_sfs(...; sample_size = 661), projected
# with MKtest.project to k = 10, 20, 50, 100, 200 individuals (2k chromosomes).
# Models: GammaZero, GammaExpo, GammaGamma, DisplGamma, ScaledBeta, under the default
# ExactSojourn integrator and every other default.
#
# Tolerance band, FIXED BEFORE THE FIRST RUN (2026-10-15):
#   1. every fit converged;
#   2. GammaZero alpha at each k within +/- 0.05 (absolute) of its k = 200 value;
#   3. GammaZero del_mean (negGmean before 2026-10-15) at each k within a factor of 2 of its
#      k = 200 value;
#   4. each other model's lnL >= the GammaZero lnL at the same k - 1e-6 (the LRT
#      sanity check: GammaZero is nested in GammaExpo, GammaGamma and DisplGamma;
#      ScaledBeta is not nested and is only reported).
#
#   cd ~/.julia/dev/MKtest_grapes && JULIA_NUM_THREADS=8 julia --project=. test/consistency_n.jl
#
# Writes test/fixtures/consistency_n.csv and prints PASS/FAIL per criterion.

using MKtest
using CSV
using DataFrames

const ROOT = normpath(joinpath(@__DIR__, ".."))
const KS = (10, 20, 50, 100, 200)
const MODELS = (MKtest.GammaZero(), MKtest.GammaExpo(), MKtest.GammaGamma(), MKtest.DisplGamma(),
                MKtest.ScaledBeta())
const NESTING_GAMMAZERO = ("GammaExpo", "GammaGamma", "DisplGamma")

sfs, div, _ = MKtest.parse_sfs(joinpath(ROOT, "data", "MKdata_tgp_may2023.txt"); sample_size = 661)
cache = MKtest.DFECache()
fits = DataFrame[]
for k in KS
    sk = MKtest.project(sfs, k)
    for m in MODELS
        t0 = time()
        res = MKtest.dfe_ml(sk, div, m; cache = cache)
        r = innerjoin(select(res.inference, Not(:alpha_ci, :omega_a_ci, :omega_na_ci)),
                      MKtest.params_wide(res); on = [:fit, :dataset, :model])
        r[!, :k] .= k
        r[!, :wall_s] .= time() - t0
        push!(fits, r)
        println("k = $k $(r.model[1]): lnL $(r.loglik[1]) alpha $(r.alpha[1]) converged $(r.converged[1])")
    end
end
all_fits = vcat(fits...; cols = :union)
src_sha = strip(read(pipeline(`bash -c "cd '$ROOT' && find src ext -type f -name '*.jl' -print0 | sort -z | xargs -0 sha256sum | sha256sum | cut -d' ' -f1"`), String))
all_fits[!, :src_sha] .= src_sha
all_fits[!, :threads] .= Threads.nthreads()
CSV.write(joinpath(ROOT, "test", "fixtures", "consistency_n.csv"), all_fits)

g0 = sort(all_fits[all_fits.model .== "GammaZero", :], :k)
ref = g0[g0.k .== 200, :]
c1 = all(all_fits.converged)
c2 = all(abs.(g0.alpha .- ref.alpha[1]) .<= 0.05)
c3 = all(max.(g0.del_mean ./ ref.del_mean[1], ref.del_mean[1] ./ g0.del_mean) .<= 2.0)
c4 = all(begin
             sub = all_fits[(all_fits.k .== k) .& in.(all_fits.model, Ref(NESTING_GAMMAZERO)), :]
             all(sub.loglik .>= g0.loglik[g0.k .== k][1] - 1e-6)
         end for k in KS)
println()
println(select(g0, :k, :n, :loglik, :alpha, :del_shape, :del_mean, :converged))
for (lbl, ok) in (("1 all converged", c1), ("2 GammaZero alpha within 0.05 of k = 200", c2),
                  ("3 GammaZero del_mean within x2 of k = 200", c3),
                  ("4 nested models lnL >= GammaZero lnL", c4))
    println(rpad(lbl, 45), ok ? "PASS" : "FAIL")
end
println("src sha $src_sha, threads $(Threads.nthreads())")

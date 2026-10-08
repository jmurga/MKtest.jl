# make_abcmk_fixture.jl -- reference outputs of the ABC-MK best model, frozen BEFORE the
# dispatch cleanup (branch abcmk-clean), at tag `abcmk-dispatch-archive`.
#
#   julia --project=. test/fixtures/make_abcmk_fixture.jl
#
# The simplified code must reproduce this file bit for bit (testset "best-model fixture").
# The best model is the one validation 4v measured: model_preset(:fixes) with PhiUniform --
# PhiPublished, BeneficialPaper (B level + sojourn at γ·B), PhiUniform, PposGuarded --
# TailWeightProjection, the paper's priors, and IndependentSampling at inference.
# Regenerate ONLY at the archive tag; after the cleanup the old keywords no longer exist.

using MKtest, JLD2, Random, DelimitedFiles

const OUT = joinpath(@__DIR__, "abcmk_best_model.jld2")
const DAC = [2, 4, 5, 10, 20, 50]
mkparam() = MKtest.parameters(N = 500, n = 50, dac = DAC)
const BEST = (phi = PhiPublished(), beneficial_scaling = BeneficialPaper(),
              phi_polymorphism = PhiUniform(), ppos_solver = PposGuarded())
const PRIORS = (alpha = [0.1, 0.9], gam_dfe = [-1000, -200], gam_flanking = [-1000, -500],
                gH = [100, 1000], gL = [1, 10])

isfile(OUT) && rm(OUT)
td = mktempdir()

# 1. rates(): 20 draws x 37 B bins
Random.seed!(20260928)
tab = joinpath(td, "t.jld2")
MKtest.rates(mkparam(); BEST..., PRIORS..., iterations = 20, rho = 0.001,
             projection = TailWeightProjection(), output_file = tab)
r = MKtest._read_rates(JLD2Backend(), tab, mkparam())

# 2. analytical_alpha under the best model, at two B
aa = Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}()
for B in (0.4, 0.999)
    p = MKtest.parameters(N = 500, n = 50, dac = DAC, B = B, al_tot = 0.4, al_low = 0.2,
                          gam_dfe = -457, gam_flanking = -500, gL = 5, gH = 500, ρ = 0.001)
    aa[B] = MKtest.analytical_alpha(p; phi = BEST.phi, scaling = BEST.beneficial_scaling,
                                    phi_poly = BEST.phi_polymorphism)
end

# 3. summary_statistics on a deterministic pseudo-dataset (model row 1's expected spectra
#    scaled to counts), both samplers, fixed seed; the written summaries are the output.
nn = 2 * 50
sfs = Matrix{Float64}(undef, nn - 1, 3)
sfs[:, 1] = 1:(nn - 1)
sfs[:, 2] .= round.(2e6 .* 0.75 ./ (1:(nn - 1)) .* 0.6)     # pn
sfs[:, 3] .= round.(2e6 .* 0.25 ./ (1:(nn - 1)))            # ps
div = [[3000.0, 2000.0, 0.75 * 2e6, 0.25 * 2e6]]
summ = Dict{String,Matrix{Float64}}()
for (tag, samp) in (("indep", IndependentSampling()), ("nested", NestedSampling()))
    of = joinpath(td, tag)
    Random.seed!(7)
    MKtest.summary_statistics(mkparam(), [sfs], div; h5_file = tab, summstat_size = 200,
                              output_folder = of, sampling = samp)
    summ[tag] = readdlm(joinpath(of, "summstat_1.txt"))
end

jldsave(OUT; models = Matrix(r.models), model_names = names(r.models), dsdn = r.dsdn,
        neut = Dict(k => r.neut[k] for k in r.dac), sel = Dict(k => r.sel[k] for k in r.dac),
        dac = r.dac, options = r.options, aa_04 = aa[0.4], aa_0999 = aa[0.999],
        summ_indep = summ["indep"], summ_nested = summ["nested"],
        sfs = sfs, div = div[1],
        commit = readchomp(`git -C $(pkgdir(MKtest)) rev-parse HEAD`))
println("wrote ", OUT, "  (", filesize(OUT), " bytes)")
println("options: ", r.options)
println("rows: ", size(r.models, 1), "  summ sizes: ", size(summ["indep"]), " ", size(summ["nested"]))

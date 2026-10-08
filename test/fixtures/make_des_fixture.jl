# make_des_fixture.jl -- reference outputs of the ABC-MK configuration RjSampling was scored with:
# a rates() table with ρ drawn per model, summarised with RjSampling.
#
# The script runs against the development branch at 8c68fa35, where RjSampling is still
# called DemographyExponentSampling; the code below therefore uses that name.
#
# Generated with the development branch at commit 8c68fa35, the code that produced the
# RjSampling evidence, single-threaded:
#
#   git worktree add --detach /tmp/mk_8c68 8c68fa35
#   julia -t 1 --project=/tmp/mk_8c68 test/fixtures/make_des_fixture.jl
#
# This branch must reproduce the file bit for bit (testset "RjSampling fixture").

using MKtest, JLD2, Random, DelimitedFiles

const OUT = joinpath(@__DIR__, "abcmk_des.jld2")
const DAC = [2, 4, 5, 10, 20, 50]
mkparam() = MKtest.parameters(N = 500, n = 50, dac = DAC)
const BEST = load(joinpath(@__DIR__, "abcmk_best_model.jld2"))

isfile(OUT) && rm(OUT)
td = mktempdir()

# 1. rates() with ρ drawn per model: 20 draws x 37 B bins, defaults otherwise
Random.seed!(20261008)
tab = joinpath(td, "rho.jld2")
MKtest.rates(mkparam(); iterations = 20, rho = [0.0001, 0.004], output_file = tab)
r = MKtest._read_rates(MKtest.JLD2Backend(), tab, mkparam())

# 2. summary_statistics with DES on the best-model fixture's pseudo-dataset: the default
#    w prior, and a fixed w (no w draw, so the composition path is pinned on its own)
summ = Dict{String,Matrix{Float64}}()
for (tag, samp) in (("des", MKtest.DemographyExponentSampling()),
                    ("des_w08", MKtest.DemographyExponentSampling(w = (0.8, 0.8))))
    of = joinpath(td, tag)
    Random.seed!(7)
    MKtest.summary_statistics(mkparam(), [BEST["sfs"]], [BEST["div"]]; h5_file = tab,
                              summstat_size = 200, output_folder = of, sampling = samp)
    summ[tag] = readdlm(joinpath(of, "summstat_1.txt"))
end

jldsave(OUT; models = Matrix(r.models), model_names = names(r.models), dsdn = r.dsdn,
        neut = Dict(k => r.neut[k] for k in r.dac), sel = Dict(k => r.sel[k] for k in r.dac),
        dac = r.dac, options = r.options, summ_des = summ["des"], summ_des_w08 = summ["des_w08"],
        commit = readchomp(`git -C $(pkgdir(MKtest)) rev-parse HEAD`))
println("wrote ", OUT, "  (", filesize(OUT), " bytes)")
println("rows: ", size(r.models, 1), "  summ sizes: ", size(summ["des"]), " ", size(summ["des_w08"]))

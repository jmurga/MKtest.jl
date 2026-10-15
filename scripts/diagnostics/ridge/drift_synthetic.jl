# Is alpha's rise with n a property of GammaZero + projection, or of the tgp data?
# Data drawn from the engine's own GammaZero expectation (shape 0.09, mean 1650, close to
# the tgp n = 200 fit) at 400 chromosomes, Ln/Ls and theta like tgp, Dn/Ds = tgp's.
# A self-consistency smoke check (oracle 4 style), 2 seeds; NOT evidence of accuracy.
using MKtest, Random, Distributions
sfs, div, _ = MKtest.parse_sfs(joinpath(@__DIR__, "..", "..", "..", "data/MKdata_tgp_may2023.txt"); sample_size = 661)
div = Float64.(vec(div)); Ln, Ls = div[3], div[4]
pre = MKtest.get_precomputed(MKtest.DFECache(), 400)
dens = zeros(length(pre.s_grid)); MKtest.dfe_density!(dens, MKtest.GammaZero(), [0.09, 1650.0], pre.s_grid)
sel = [sum(dens .* pre.s_widths .* pre.expected_counts[:, j]) for j in 1:399]
theta = 0.0014
for seed in (1, 2)
    rng = MersenneTwister(seed)
    pn = [rand(rng, Poisson(theta * Ln * sel[j])) for j in 1:399]
    ps = [rand(rng, Poisson(theta * Ls / j)) for j in 1:399]
    s400 = hcat(collect(1.0:399.0), Float64.(pn), Float64.(ps))
    # divergence is tgp's own, so alpha is not a known truth here: only its flatness in n is read
    out = String[]
    for k in (10, 20, 50, 100, 200)
        s = k == 200 ? s400 : MKtest.project(s400, k)
        r = MKtest.dfe_ml(s, div, MKtest.GammaZero())
        push!(out, "k=$k shape $(round(MKtest.param_values(r, :del_shape)[1]; digits=3)) mean $(round(Int, MKtest.param_values(r, :del_mean)[1])) alpha $(round(r.inference.alpha[1]; digits=4))")
    end
    println("seed $seed (truth shape 0.09, mean 1650):\n  ", join(out, "\n  "))
end

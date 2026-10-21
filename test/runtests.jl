using MKtest
using Test
using DataFrames: nrow, names, propertynames, Not
using JLD2
using Random: seed!, MersenneTwister
using DelimitedFiles: readdlm, writedlm
import Logging

@testset "MKtest.jl" begin

    # --- rho prior in rates() ------------------------------------------------
    # `ρ` reaches the model through `set_θ!` (the flanking mutation rate that produces
    # the target B) and `Φ`. `rates(...; ρ = [lo, hi])` draws one value per model and
    # records it in the `rho` column. The default `rho = nothing` consumes NO RNG.
    @testset "rho prior" begin
        param  = MKtest.parameters(N=500, n=100, dac=[2, 4, 5, 10, 20, 50],
                                   cutoff=[0.0, 1.0], B=0.999)
        B_bins = [0.2, 0.5, 0.999]
        nb     = length(B_bins)
        od     = mktempdir()
        counter = Ref(0)
        function run_rates(rho; iterations=60)
            counter[] += 1
            f = joinpath(od, "rates_$(counter[]).jld2")
            seed!(7)
            r = MKtest.rates(param; iterations=iterations, output_file=f, B_bins=B_bins,
                             ρ=rho)
            return r, f
        end

        # one draw per model, shared by that model's B bins, inside the range
        r, f = run_rates([1e-4, 5e-3])
        m = r[:models]
        @test :rho in propertynames(m)
        @test names(m)[9] == "rho"          # appended last: cols 1-8 keep their meaning
        @test all(1e-4 .<= m.rho .<= 5e-3)
        nmodels = size(m, 1) ÷ nb
        @test all(allequal(m.rho[(i-1)*nb+1:i*nb]) for i in 1:nmodels)
        @test length(unique(m.rho)) == nmodels
        @test names(MKtest._read_rates(f, param).models)[9] == "rho"

        # fixed forms
        @test all(==(0.002), run_rates(0.002; iterations=20)[1][:models].rho)
        @test all(==(0.002), run_rates([0.002]; iterations=20)[1][:models].rho)
        @test all(==(param.ρ), run_rates(nothing; iterations=20)[1][:models].rho)

        # invalid ranges must throw rather than silently pick something
        for bad in ([-1.0, 1.0], [5e-3, 1e-4], -0.001, Float64[])
            @test_throws ArgumentError run_rates(bad; iterations=5)
        end

        # the recorded rho is the one set_θ! actually used
        p = MKtest.copy_param(param)
        p.B = m.B[1] / 1000
        p.ρ = m.rho[1]
        p.γ_flanking = Int(m.gam_flanking[1])
        @test MKtest.set_θ!(p) ≈
              -log(p.B) * (p.L_flanking * m.rho[1] - m.gam_flanking[1]) / (2.0 * p.L_flanking)
    end

    # --- the projection the rates table is built with -----------------------
    @testset "tail-weight projection" begin
        param  = MKtest.parameters(N=500, n=100, dac=[2, 4, 5, 10, 20, 50],
                                   cutoff=[0.0, 1.0], B=0.999)
        B_bins = [0.2, 0.5, 0.999]
        od = mktempdir()
        reldiff(a, b) = maximum(abs.(a .- b) ./ max.(abs.(b), 1e-300))

        seed!(21)
        r = MKtest.rates(param; iterations=20, output_file=joinpath(od, "tw.jld2"),
                         B_bins=B_bins)
        m = r[:models]
        # every stored row equals the full-spectrum model at the same parameters
        for i in (1, 2, 3, nrow(m))
            q = MKtest.copy_param(param)
            q.γ_strong, q.γ_weak = Int(m.gH[i]), Int(m.gL[i])
            q.γ_flanking, q.γ_deleterious = Int(m.gam_flanking[i]), Int(m.gam_dfe[i])
            q.shape = m.shape[i]; q.scale = abs(q.shape / q.γ_deleterious)
            q.α, q.α_weak, q.ρ = m.al_tot[i], m.al_low[i], m.rho[i]
            q.B = m.B[i] / 1000; MKtest.set_θ!(q)
            # ppos recovered EXACTLY from the stored fixation rates:
            # f_pos = 0.75 * (ppos/2) * Φ * p_fix
            pp(f, g) = f == 0 ? 0.0 : 2f / (0.75 * MKtest.Φ(q, g) * MKtest.p_fix(q, g))
            q.p_weak = pp(r[:dsdn][i, 3], q.γ_weak)
            q.p_strong = pp(r[:dsdn][i, 4], q.γ_strong)
            bn = MKtest.binom_op(q)
            neut = MKtest.cumulative_vector!(MKtest.sfs_neut(q, bn))
            sel  = MKtest.cumulative_vector!(MKtest.sfs_pos(q, q.γ_strong, q.p_strong, bn) .+
                                             MKtest.sfs_pos(q, q.γ_weak, q.p_weak, bn) .+
                                             MKtest.sfs_neg(q, q.p_weak + q.p_strong, bn))
            @test reldiff([r[:neut][k][i] for k in param.dac], neut[param.dac]) < 1e-12
            @test reldiff([r[:sel][k][i]  for k in param.dac], sel[param.dac])  < 1e-12
        end

        # The tail-weight cache is keyed on the DAC vector as well as B.
        p3 = MKtest.parameters(N=500, n=100, dac=[4, 20, 50], cutoff=[0.0, 1.0], B=0.999)
        seed!(5)
        r6 = MKtest.rates(param; iterations=30, output_file=joinpath(od, "k6.jld2"),
                          B_bins=B_bins)
        seed!(5)
        r3 = MKtest.rates(p3; iterations=30, output_file=joinpath(od, "k3.jld2"),
                          B_bins=B_bins)
        for k in keys(r3[:neut])
            @test reldiff(r3[:neut][k], r6[:neut][k]) < 1e-12
            @test reldiff(r3[:sel][k],  r6[:sel][k])  < 1e-12
        end

        # the identity the whole projection rests on
        binom = MKtest.binom_op(param, B_bins)
        tw = MKtest.tail_weights_op(param, binom, B_bins)
        A = binom[0.999]
        W = tw[0.999].W
        src = rand(size(A, 2))
        ref = MKtest.cumulative_vector!(A * src)[param.dac]
        @test reldiff(W * src, ref) < 1e-12
        NN2 = size(A, 2) - 1
        h = [j == 0 ? 0.0 : 1.0 / j for j in 0:NN2]
        @test reldiff(tw[0.999].neut_tail,
                      MKtest.cumulative_vector!(A * h)[param.dac]) < 1e-12
    end

    # --- set_ppos! guard and rates provenance --------------------------------
    @testset "ppos solver and rates provenance" begin
        mkp() = MKtest.parameters(N = 500, n = 50, dac = [2, 4, 5, 10, 20, 50])

        # A feasible target is solved to the residual tolerance.
        pa = mkp(); pa.α = 0.4; pa.α_weak = 0.2; MKtest.set_θ!(pa)
        MKtest.set_ppos!(pa)
        @test maximum(MKtest._ppos_residual(pa, pa.p_weak, pa.p_strong)) < 1e-6
        @test pa.p_weak > 0 && pa.p_strong > 0

        # al_low == 0 means ppos_l == 0 IS the root, not a failure.
        pz = mkp(); pz.α = 0.4; pz.α_weak = 0.0; MKtest.set_θ!(pz)
        MKtest.set_ppos!(pz)
        @test pz.p_weak < 1e-20
        @test maximum(MKtest._ppos_residual(pz, pz.p_weak, pz.p_strong)) < 1e-6

        # A very peaked deleterious DFE fixes almost nothing, so a low al_tot needs a tiny
        # ppos that nlsolve misses from every start. The system is linear in ppos and its
        # exact root is accepted.
        pj = mkp(); pj.shape = 7.5; pj.γ_deleterious = -500
        pj.scale = abs(pj.shape / pj.γ_deleterious)
        pj.α = 0.06; pj.α_weak = 0.054; MKtest.set_θ!(pj)
        @test MKtest._set_ppos!(pj) === :linear
        @test maximum(MKtest._ppos_residual(pj, pj.p_weak, pj.p_strong)) < 1e-12
        @test 0 < pj.p_weak < 1e-10 && 0 < pj.p_strong < pj.p_weak
        @test MKtest.set_ppos!(pj) == (pj.p_weak, pj.p_strong)
        # Oracle for the closed form: the fixation alpha of its root, recomputed by hand
        # from the three fixation terms.
        fw = MKtest.fix_pos_sim(pj, pj.γ_weak, 0.5pj.p_weak)
        fs = MKtest.fix_pos_sim(pj, pj.γ_strong, 0.5pj.p_strong)
        fn = MKtest.fix_neg(pj, 0.5(pj.p_weak + pj.p_strong))
        @test (fw + fs) / (fw + fs + fn) ≈ 0.06 rtol = 1e-10
        @test fw / (fw + fs + fn) ≈ 0.054 rtol = 1e-10
        # Where nlsolve converges it keeps the root (and its bits); the closed form is the
        # same root to nlsolve's tolerance.
        @test MKtest._set_ppos!(pa) === :nlsolve
        @test sum(MKtest._ppos_linear(pa)) ≈ pa.p_weak + pa.p_strong rtol = 1e-6

        # No finite root: refused. shape == 1 is the pole of the Hurwitz zeta in fix_neg.
        p1 = mkp(); p1.shape = 1.0; p1.γ_deleterious = -500
        p1.scale = abs(p1.shape / p1.γ_deleterious)
        p1.α = 0.06; p1.α_weak = 0.054; MKtest.set_θ!(p1)
        @test_throws ErrorException MKtest.set_ppos!(p1)

        # the model record is written next to the table
        seed!(1)
        f = joinpath(mktempdir(), "prov.jld2")
        MKtest.rates(mkp(); iterations = 8, output_file = f, ρ = 0.001)
        o = MKtest._read_rates(f, mkp()).options
        @test o == MKtest._rates_options()
        @test startswith(o, "model=best-2026-09-28")
    end

    # --- the forward model: what each choice is -------------------------------
    @testset "best-model physics" begin
        p = MKtest.parameters(N = 7310, n = 661, γ_deleterious = -457,
                              γ_flanking = -457, L_flanking = 200_000,
                              γ_weak = 10, γ_strong = 500, ρ = 0.001, B = 0.4,
                              dac = [2, 4, 5, 10, 20, 50, 200, 661, 925])
        MKtest.set_θ!(p)

        # Φ = B^κ follows from set_θ!'s closed form; the published expression carries
        # the factor of two on the beneficial coefficient (κ_L/κ_H = 8.2, not 3.8).
        κ(g) = log(MKtest.Φ(p, g)) / log(p.B)
        @test isapprox(κ(10), 0.5455; atol = 1e-3)
        @test isapprox(κ(500), 0.0662; atol = 1e-3)

        # sfs_pos is the paper's "replace N with Ne in the frequency spectra": a literal
        # transcription -- level B, sojourn exp-form at γ·B, Φ flat.
        q = MKtest.parameters(N = 1000, n = 661, γ_deleterious = -457, γ_flanking = -500,
                              L_flanking = 200_000, γ_weak = 5, γ_strong = 500, ρ = 0.001, B = 0.4,
                              dac = [4, 5, 10, 20, 50, 200, 661, 925])
        q.p_weak, q.p_strong = 0.0077070856, 0.00015930139
        MKtest.set_θ!(q)
        function sfs_pos_ref(q, s, pp, binom)
            NN2 = ceil(Int, q.NN * q.B)
            xs = (0:NN2) ./ NN2
            sc = s * q.B
            g1, g2 = exp(sc * 2), exp(sc * -2)
            v = [pp * 0.5 * (g1 * (1 - g2^(1.0 - x)) / ((g1 - 1.0) * x * (1.0 - x))) for x in xs]
            v[1] = 0.0; v[end] = 0.0
            return q.B * q.θ_coding * MKtest.Φ(q, s) * 0.75 * (binom * (v ./ NN2))
        end
        bn = MKtest.binom_op(q)
        for (s, pp) in ((q.γ_weak, q.p_weak), (q.γ_strong, q.p_strong))
            @test isapprox(MKtest.sfs_pos(q, s, pp, bn), sfs_pos_ref(q, s, pp, bn); rtol = 1e-10)
        end

        # The B level makes beneficial polymorphism relative to neutral FALL as BGS
        # strengthens (without it the ratio rose, 1.486; SLiM measures 0.637).
        function sel_over_neut(B)
            r = MKtest.copy_param(q); r.B = B; MKtest.set_θ!(r)
            b = MKtest.binom_op(r)
            nt = MKtest.cumulative_vector!(MKtest.sfs_neut(r, b))
            sl = MKtest.cumulative_vector!(MKtest.sfs_pos(r, r.γ_weak, r.p_weak, b))
            sl[925] / nt[925]
        end
        ratio = sel_over_neut(0.4) / sel_over_neut(0.999)
        @test 0.3 < ratio < 1.0

        # γ·B against a sojourn at γ: large at low B, vanishing as B -> 1.
        function src_rel(B)
            r = MKtest.copy_param(q); r.B = B
            NN2 = ceil(Int, r.NN * B); sc = zeros(NN2 + 1); sc2 = zeros(NN2 + 1)
            a = MKtest._src_pos!(sc,  r, r.γ_weak * B, r.p_weak, NN2 + 1)
            b = MKtest._src_pos!(sc2, r, r.γ_weak,     r.p_weak, NN2 + 1)
            maximum(abs.(a[2:end-1] .- b[2:end-1]) ./ abs.(b[2:end-1]))
        end
        @test src_rel(0.4) > 0.05
        @test src_rel(0.999) < 1e-2
    end

    # --- syn_nonsyn_ratio reaches the table ------------------------------------
    # `rates()` solves each task on `copy_param(param)`, which lists the fields one by
    # one: a field it misses silently falls back to its default.
    @testset "f_syn" begin
        mkp(x) = MKtest.parameters(N = 500, n = 50, dac = [2, 4, 5, 10, 20, 50],
                                   f_syn = x)
        @test MKtest.copy_param(mkp(0.3)).f_syn == 0.3
        od = mktempdir()
        seed!(3)
        a = MKtest.rates(mkp(0.25); iterations = 5, output_file = joinpath(od, "a.jld2"))
        seed!(3)
        b = MKtest.rates(mkp(0.3); iterations = 5, output_file = joinpath(od, "b.jld2"))
        @test b[:models][:, 1:8] == a[:models][:, 1:8]
        # the neutral SFS is linear in the ratio
        @test all(isapprox(b[:neut][k], a[:neut][k] .* (0.3 / 0.25); rtol = 1e-12) for k in [2, 50])
    end

    # --- the DFE shape prior is derived from param.shape ----------------------
    @testset "DFE shape prior is param-derived" begin
        od = mktempdir()
        lo, hi = 0.184 * 2.0^-2, 0.184 * 2.0^2          # [0.046, 0.736]
        seed!(7)
        out = MKtest.rates(MKtest.parameters(N = 500, n = 50, shape = 0.184,
                                             dac = [2, 4, 5, 10, 20, 50]);
                           iterations = 10, output_file = joinpath(od, "s1.jld2"), ρ = 0.001)
        @test minimum(out[:models].shape) >= lo - 1e-12
        @test maximum(out[:models].shape) <= hi + 1e-12
        seed!(7)
        out2 = MKtest.rates(MKtest.parameters(N = 500, n = 50, shape = 0.4,
                                              dac = [2, 4, 5, 10, 20, 50]);
                            iterations = 10, output_file = joinpath(od, "s2.jld2"), ρ = 0.001)
        @test maximum(out2[:models].shape) > hi
    end

    @testset "shape_factor widens the DFE shape prior" begin
        # f = 2 is the published window exactly (powers of two), so old tables keep bits.
        u2 = MKtest._shape_prior(0.184, 2)
        @test (u2.a, u2.b) == (0.184 * 2^-2, 0.184 * 2^2)
        u3 = MKtest._shape_prior(0.184, 3)
        @test u3.a ≈ 0.184 / 9 && u3.b ≈ 0.184 * 9
        @test_throws ArgumentError MKtest._shape_prior(0.184, 0.5)
        @test_throws ArgumentError MKtest._shape_prior(0.184, Inf)
        od = mktempdir()
        mkp() = MKtest.parameters(N = 500, n = 50, shape = 0.184, dac = [2, 4, 5, 10, 20, 50])
        seed!(11)
        d = MKtest.rates(mkp(); iterations = 10, output_file = joinpath(od, "d.jld2"), ρ = 0.001)
        seed!(11)
        d2 = MKtest.rates(mkp(); iterations = 10, output_file = joinpath(od, "d2.jld2"),
                          ρ = 0.001, shape_factor = 2)
        @test d[:sel] == d2[:sel] && d[:dsdn] == d2[:dsdn]
        # f = 4: [0.0115, 2.944]; draws above shape 1 are solved, not refused
        seed!(11)
        w = MKtest.rates(mkp(); iterations = 40, output_file = joinpath(od, "w.jld2"),
                         ρ = 0.001, shape_factor = 4)
        sh = w[:models].shape
        @test minimum(sh) >= 0.184 / 16 - 1e-12 && maximum(sh) <= 0.184 * 16 + 1e-12
        @test count(>(1), sh) >= 10
        @test size(w[:models], 1) == 40 * 37               # nothing refused
    end

    # --- reference outputs ----------------------------------------------------
    # The fixtures were produced by the code (commit 8c68fa35) that RjSampling was
    # tested with: abcmk_best_model.jld2 for the forward model and the default sampler,
    # abcmk_des.jld2 for a ρ-drawn table and RjSampling (make_des_fixture.jl). Both must
    # be reproduced at one thread.
    function check_table(r, fx, DAC)
        @test Matrix(r.models) == fx["models"]
        @test names(r.models) == fx["model_names"]
        # At ONE thread the rates match to the last bits: another CPU (SIMD width, FMA)
        # changes the summation order and flips the last digit. With more threads, `dsdn` /
        # `sel` drift at the 1e-9 level (the per-task ppos warm start depends on how rows
        # are chunked).
        tol = Threads.nthreads() == 1 ? 1e-10 : 1e-7
        @test r.dsdn ≈ fx["dsdn"] rtol = tol
        @test all(isapprox(r.neut[k], fx["neut"][k]; rtol = tol) &&
                  isapprox(r.sel[k], fx["sel"][k]; rtol = tol) for k in DAC)
    end

    @testset "best-model fixture" begin
        fx = JLD2.load(joinpath(@__DIR__, "fixtures", "abcmk_best_model.jld2"))
        DAC = [2, 4, 5, 10, 20, 50]
        mkparam() = MKtest.parameters(N = 500, n = 50, dac = DAC)
        td = mktempdir()
        seed!(20260928)
        tab = joinpath(td, "t.jld2")
        MKtest.rates(mkparam(); iterations = 20, ρ = 0.001, output_file = tab)
        r = MKtest._read_rates(tab, mkparam())
        check_table(r, fx, DAC)

        for (B, key) in ((0.4, "aa_04"), (0.999, "aa_0999"))
            p = MKtest.parameters(N = 500, n = 50, dac = DAC, B = B, α = 0.4,
                                  α_weak = 0.2, γ_deleterious = -457, γ_flanking = -500,
                                  γ_weak = 5, γ_strong = 500, ρ = 0.001)
            @test MKtest.analytical_alpha(p) == fx[key]
        end

        # the independent sampler with uniform rows (ABC-MK v1.0.0)
        of = joinpath(td, "indep")
        seed!(7)
        MKtest.summary_statistics(mkparam(), [fx["sfs"]], [fx["div"]]; h5_file = tab,
                                  summstat_size = 200, output_folder = of,
                                  sampling = MKtest.IndependentSampling(),
                                  row_prior = MKtest.TableRowPrior())
        @test MKtest._with_gz(readdlm, joinpath(of, "summstat_1.txt.gz")) == fx["summ_indep"]

        # the defaults are RjSampling(w = (0, 2.5), tilt = 0.05) with FlatAlphaPrior(by_B = true)
        summ(; kw...) = begin
            of = mktempdir()
            seed!(7)
            MKtest.summary_statistics(mkparam(), [fx["sfs"]], [fx["div"]]; h5_file = tab,
                                      summstat_size = 200, output_folder = of, kw...)
            MKtest._with_gz(readdlm, joinpath(of, "summstat_1.txt.gz"))
        end
        @test summ() == summ(sampling = MKtest.RjSampling(w = (0.0, 2.5), tilt = 0.05),
                             row_prior = MKtest.FlatAlphaPrior(by_B = true))
        @test summ() != fx["summ_indep"]
    end

    @testset "RjSampling fixture" begin
        fx  = JLD2.load(joinpath(@__DIR__, "fixtures", "abcmk_des.jld2"))
        best = JLD2.load(joinpath(@__DIR__, "fixtures", "abcmk_best_model.jld2"))
        DAC = [2, 4, 5, 10, 20, 50]
        mkparam() = MKtest.parameters(N = 500, n = 50, dac = DAC)
        td = mktempdir()
        seed!(20261008)
        tab = joinpath(td, "rho.jld2")
        MKtest.rates(mkparam(); iterations = 20, ρ = [0.0001, 0.004], output_file = tab)
        r = MKtest._read_rates(tab, mkparam())
        check_table(r, fx, DAC)
        @test r.options == fx["options"]

        # recorded with the untilted sampler, w ~ U(0, 1.5), and uniform rows
        for (tag, samp) in (("summ_des", MKtest.RjSampling(w = (0.0, 1.5), tilt = 0.0)),
                            ("summ_des_w08", MKtest.RjSampling(w = (0.8, 0.8), tilt = 0.0)))
            # the summaries are drawn from the table built here; at one thread it is the
            # fixture's table, so the draws must match exactly
            Threads.nthreads() == 1 || continue
            of = joinpath(td, tag)
            seed!(7)
            MKtest.summary_statistics(mkparam(), [best["sfs"]], [best["div"]]; h5_file = tab,
                                      summstat_size = 200, output_folder = of, sampling = samp,
                                      row_prior = MKtest.TableRowPrior())
            @test MKtest._with_gz(readdlm, joinpath(of, "summstat_1.txt.gz")) == fx[tag]
        end
    end

    # --- parse_sfs ---------------------------------------------------------------
    @testset "parse_sfs" begin
        tgp = joinpath(@__DIR__, "..", "data", "tgp.txt")
        sfs, d, a = MKtest.parse_sfs(tgp; sample_size = 661)
        @test size(sfs) == (1321, 3) && sfs[:, 1] == 1:1321
        @test length(d) == 4                          # Dn, Ds, Ln, Ls
        @test size(a) == (1321, 3) && a[:, 2] ≈ a[:, 1] ./ 1322
        # a cutoff removes rows but keeps the true DAC and frequency
        sfs7, _, a7 = MKtest.parse_sfs(tgp; sample_size = 661, cutoff = [0.0, 0.7])
        @test sfs7[end, 1] == 925 && a7[end, 2] ≈ 925 / 1322

        # control sets (several genes per row) give what get_pol_div gives for the same
        # genes, also below 2n = 100 where frequencies are rounded to 3 digits
        rng = MersenneTwister(1)
        daf(nn, k; d = 4) = "," * join(string.(round.(rand(rng, 1:(nn - 1), k) ./ nn, digits = d)), ",") * ","
        tables = Dict(
            1322 => MKtest.CSV.read(tgp, MKtest.DataFrame; header = false, delim = '\t')[1:30, :],
            60   => MKtest.DataFrame(Column1 = ["g$i" for i in 1:30], Column2 = fill(5, 30),
                                     Column3 = [daf(60, 5; d = 3) for _ in 1:30], Column4 = fill(5, 30),
                                     Column5 = [daf(60, 5; d = 3) for _ in 1:30],
                                     Column6 = rand(rng, 1:9, 30), Column7 = rand(rng, 1:9, 30)))
        for (nn, df) in tables
            genes = string.(df[:, 1])
            rows = [genes[1:10], genes[11:20], genes[21:30]]
            gl = joinpath(mktempdir(), "control.txt")
            write(gl, join(join.(rows, " "), "\n") * "\n")
            s_c, d_c, a_c = MKtest.parse_sfs(df; sample_size = nn, ploidy = 1, gene_list = gl)
            for (i, r) in enumerate(rows)
                sub = df[in(r).(string.(df[:, 1])), :]
                a_r, s_r, d_r = MKtest.get_pol_div(sub, nn, [0.0, 1.0], [3, 5], [6, 7], [8, 9])
                @test s_c[i] == s_r
                @test d_c[i] == d_r
                @test isequal(a_c[i], a_r)
            end
        end
    end

    # --- MK tests -------------------------------------------------------------------
    @testset "MK tests" begin
        data = joinpath(@__DIR__, "..", "data")
        tgp = joinpath(data, "tgp.txt")
        sfs, d, a = MKtest.parse_sfs(tgp; sample_size = 661)
        # parse_sfs rounds alpha(x) to 5 digits
        @test MKtest.α_x(sfs, d)[:, 1:2] == a[:, 1:2]
        @test all(isapprox.(MKtest.α_x(sfs, d)[:, 3], a[:, 3]; atol = 1e-5, nans = true))

        # the model is chosen, never switched: the exponential is the default
        r_exp = MKtest.asymptotic_mk(sfs, d)
        r_lin = MKtest.asymptotic_mk(sfs, d; model = MKtest.LinearFit())
        @test r_exp.model == :exponential && r_lin.model == :linear
        @test r_lin.c == 0.0 && r_lin.asymptotic_mk ≈ r_lin.a + r_lin.b
        for r in (r_exp, r_lin)
            @test r.ci_low < r.asymptotic_mk < r.ci_high
        end
        @test MKtest.asymptotic_mk(a).asymptotic_mk ≈ r_exp.asymptotic_mk atol = 1e-5

        # control set 13 of example_bootstrap.txt: the exponential does not converge, so
        # ExponentialFit returns NaN instead of falling back to the line
        sc, dc, _ = MKtest.parse_sfs(tgp; sample_size = 661,
                                     gene_list = joinpath(data, "example_bootstrap.txt"))
        @test isnan(MKtest.asymptotic_mk(sc[13], dc[13]).asymptotic_mk)
        @test isfinite(MKtest.asymptotic_mk(sc[13], dc[13]; model = MKtest.LinearFit()).asymptotic_mk)
        df = MKtest.asymptotic_mk(sc[1:20], dc[1:20]; model = MKtest.LinearFit())
        @test size(df, 1) == 20 && all(df.model .== :linear)

        # fewer than three alpha(x) values: NaN for either model
        @test isnan(MKtest.asymptotic_mk([0.1, 0.2]; model = MKtest.LinearFit()).asymptotic_mk)

        # the other tests return one row per set, and p-value correction adds a column
        for fn in (MKtest.imputed_mk, MKtest.fww_mk, MKtest.standard_mk, MKtest.ni, MKtest.dos)
            @test size(fn(sc[1:10], dc[1:10]), 1) == 10
        end
        @test :pvalue_adjusted in propertynames(MKtest.standard_mk(sc[1:10], dc[1:10]; correction = :fdr))
        @test MKtest.standard_mk(sfs, d).alpha ≈ 1 - (sum(sfs[:, 2]) / sum(sfs[:, 3])) * (d[2] / d[1]) atol = 1e-3

        # a trimmed SFS (DAC/n <= 0.7) with its n gives what the full SFS gives with the
        # trimmed rows set to 0: the frequencies come from the DAC, not from the row count
        keep = sfs[:, 1] ./ 1322 .<= 0.7
        trimmed, zeroed = sfs[keep, :], copy(sfs)
        zeroed[.!keep, 2:3] .= 0
        for fn in (MKtest.fww_mk, MKtest.imputed_mk)
            @test isequal(fn(trimmed, d; n = 1322), fn(zeroed, d))
            # trimmed at the low end the SFS cannot pass for a full one: n is required (a
            # top-only trim, DAC 1..k, is indistinguishable from a full SFS of k + 1)
            @test_throws ArgumentError fn(sfs[3:end, :], d)
        end
        # a bin exactly at the cutoff is high-frequency only (n = 20: DAC 3 is 0.15)
        s20 = hcat(1:19, fill(2.0, 19), fill(4.0, 19))
        s20[1:3, 2] .= 6.0                       # excess low-frequency Pn, DAC 3 included
        lo = s20[:, 1] ./ 20 .< 0.15
        @test !lo[3]
        pn_lo, ps_lo = sum(s20[lo, 2]), sum(s20[lo, 3])
        pn_hi, ps_hi = sum(s20[.!lo, 2]), sum(s20[.!lo, 3])
        del = pn_lo - pn_hi * ps_lo / ps_hi
        r20 = MKtest.imputed_mk(s20, [30.0, 40.0, 1e5, 3e4])
        @test 0 < del < sum(s20[:, 2])
        @test r20.alpha ≈ round(1 - ((sum(s20[:, 2]) - del) / sum(s20[:, 3])) * (40 / 30), digits = 3)
        # Ps = 0, Dn = 0 and missing site counts give NaN, not an error or an Int 0
        @test isnan(MKtest.imputed_mk(hcat(1:19, ones(19), zeros(19)), d).alpha)
        @test all(isnan, values(MKtest.fww_mk(sfs, d[1:2]))[3:5])
        @test all(isnan, values(MKtest.standard_mk(sfs, d[1:2]))[3:5])
    end

    @testset "reduce_sfs" begin
        # n = 10 chromosomes, 2 bins: DAC 1..5 (x <= 0.5) and 6..9
        s = hcat(1:9, ones(9), 2 .* ones(9))
        @test MKtest.reduce_sfs(s, 2) == [0.5 5 10; 1.0 4 8]
        # empty bins are 0; totals are kept
        r = MKtest.reduce_sfs(s, 20)
        @test size(r) == (20, 3) && sum(r[:, 2]) == 9 && sum(r[:, 3]) == 18
        @test r[1, 2] == 0 && r[2, 2] == 1          # x = 0.1 falls in (0.05, 0.1]
        # a trimmed SFS needs n; the bins are then the same as for the full one
        @test MKtest.reduce_sfs(s[1:4, :], 2; n = 10) == [0.5 4 8; 1.0 0 0]
        @test MKtest.reduce_sfs(s[6:9, :], 2; n = 10) == [0.5 0 0; 1.0 4 8]
        @test_throws ArgumentError MKtest.reduce_sfs(s[6:9, :], 2)
        # the binned table is accepted by the MK tests as frequencies
        @test MKtest.fww_mk(MKtest.reduce_sfs(s, 10), [30.0, 40.0]).alpha ==
              MKtest.fww_mk(s, [30.0, 40.0]).alpha
    end

    # --- ABC output order --------------------------------------------------------
    # Posteriors must come back in the order of the input datasets. `readdir` sorts names
    # as text, so without explicit ordering dataset 10 came second.
    @testset "ABC output order" begin
        @test MKtest._dataset_index("/a/b/summstat_10.txt") == 10
        @test MKtest._dataset_index("out_7.0.post.gz") == 7
        @test MKtest._dataset_index("summstat_10.txt.gz") == 10
        @test MKtest._as_count(1e3, "tol") === 1000
        @test MKtest._as_count(Int32(5), "tol") === 5
        @test_throws ArgumentError MKtest._as_count(0.5, "tol")
        @test_throws ArgumentError MKtest._as_count(0, "tol")
        td = mktempdir()
        for k in (1, 2, 10, 11, 3)
            touch(joinpath(td, "alphas_$k.txt"))
        end
        @test basename.(MKtest._dataset_files(td, "alphas_")) ==
              ["alphas_$k.txt" for k in (1, 2, 3, 10, 11)]

        # end to end: dataset k has a smaller Ds as k grows, so its alpha estimate rises
        best = JLD2.load(joinpath(@__DIR__, "fixtures", "abcmk_best_model.jld2"))
        p() = MKtest.parameters(N = 500, n = 50, dac = [2, 4, 5, 10, 20, 50])
        tab = joinpath(td, "t.jld2")
        seed!(1)
        MKtest.rates(p(); iterations = 300, output_file = tab)
        K = 12
        sfs = [best["sfs"] for _ in 1:K]
        div = [[3000.0, 2000.0 * (1 - 0.06k), best["div"][3], best["div"][4]] for k in 1:K]
        # ABC-MK v1.0.0 settings: this checks the file order, not the estimator. With the
        # defaults the two lowest-α datasets (α near 0) can swap places.
        # summstat_size and tol given as Int for one method and as Float for the other
        for (method, ss, tl) in (("abcreg", 5000, 200), ("abc", 5e3, 200.0))
            seed!(2)
            res = MKtest.abcmk(p(), sfs, div; h5_file = tab, summstat_size = ss,
                               output_folder = joinpath(td, method), tol = tl,
                               abc_method = method, sampling = MKtest.IndependentSampling(),
                               row_prior = MKtest.TableRowPrior(), stat = MKtest.PosteriorMode())
            @test issorted(res[:inference].α)
        end

        # a single set as parse_sfs returns it (Matrix, Vector) equals the wrapped call
        one(s, d) = (seed!(3); MKtest.abcmk(p(), s, d; h5_file = tab, summstat_size = 2000,
                                            output_folder = joinpath(td, "one"), tol = 100))
        @test one(sfs[1], div[1])[:inference] == one([sfs[1]], [div[1]])[:inference]

        # --- alpha composition: α == α_weak + α_strong through the ABC ---------------
        # α_strong is α - α_weak per Poisson draw (summary_statistics.jl), rounded to 5
        # digits (gap <= 1.5e-5 per row); the ABCreg regression is linear per column, so
        # each posterior row keeps its rounding gap, and the mean is additive. The median
        # and the KDE mode are per-column statistics and do not add up, so Composed()
        # derives the strong component of the point estimates. The driver
        # stuck_on_cell/scripts/src/alpha_composition.jl ran the same checks at full scale.
        @testset "alpha composition" begin
            cn = [:α_weak, :α_strong, :α, :ωₐ_weak, :ωₐ_strong, :ωₐ, :ωₙₐ,
                  :γ_deleterious, :β, :γ_weak, :γ_strong, :B]
            gap3(m) = maximum(abs.(m[:, 3] .- m[:, 1] .- m[:, 2]))
            # the stages of one dataset: summaries file, ABCreg posterior, point estimates
            function stages(s, d, tag)
                of = joinpath(td, tag)
                seed!(5)
                MKtest.summary_statistics(p(), [s], [d]; h5_file = tab, summstat_size = 3000,
                                          output_folder = of)
                m = MKtest.read_gz(joinpath(of, "summstat_1.txt.gz"))
                target = vec(MKtest.read_gz(joinpath(of, "alphas_1.txt.gz")))
                P = size(m, 2) - length(target)
                @test P in (8, 12)
                @test gap3(m) <= 1.5e-5
                post = MKtest.abcreg_loclinear(target, m[:, 1:P], m[:, (P + 1):end]; tol = 200)
                @test size(post, 1) == 200
                @test gap3(post) <= 2e-5
                pdf = [MKtest.DataFrame(post, P == 12 ? cn : cn[[1, 2, 3, 8, 9, 10, 11, 12]])]
                mdf = Matrix(pdf[1])
                for st in (MKtest.PosteriorMean(), MKtest.PosteriorMedian(), MKtest.PosteriorMode(),
                           MKtest.PosteriorQuantile(0.975))
                    # Composed: exact composition, α_strong derived
                    r = MKtest.summary_abc(pdf; stat = st, components = MKtest.Composed())[:inference]
                    @test r.α[1] - r.α_weak[1] - r.α_strong[1] == 0
                    @test isfinite(r.α_strong[1])
                    # Marginal (default): α_strong is the statistic of its own column
                    r = MKtest.summary_abc(pdf; stat = st)[:inference]
                    @test r.α_strong[1] == MKtest.point_estimate(st, mdf[:, 2])
                end
                # the marginal mean composes to rounding without Composed
                r = MKtest.summary_abc(pdf; stat = MKtest.PosteriorMean())[:inference]
                @test abs(r.α[1] - r.α_weak[1] - r.α_strong[1]) <= 2e-5
                return post
            end
            stages(sfs[1], div[1], "comp_fixture")

            # a synthetic set with a known composition: the forward model's expected spectra
            # and fixation rates at α 0.4, α_weak 0.2 (as analytical_alpha computes them),
            # scaled and Poisson-sampled; truth α_weak + α_strong == α by construction
            q = p(); q.α = 0.4; q.α_weak = 0.2
            binom = MKtest.binom_op(q)
            MKtest.set_θ!(q); MKtest.set_ppos!(q)
            fN   = q.B * MKtest.fix_neut(q)
            fNeg = q.B * MKtest.fix_neg(q, 0.5 * q.p_strong + 0.5 * q.p_weak)
            fL   = MKtest.fix_pos_sim(q, q.γ_weak, 0.5 * q.p_weak)
            fH   = MKtest.fix_pos_sim(q, q.γ_strong, 0.5 * q.p_strong)
            neut = MKtest.sfs_neut(q, binom)
            sel  = MKtest.sfs_pos(q, q.γ_strong, q.p_strong, binom) .+
                   MKtest.sfs_pos(q, q.γ_weak, q.p_weak, binom) .+
                   MKtest.sfs_neg(q, q.p_strong + q.p_weak, binom)
            dn, ds = fNeg + fL + fH, fN
            truth = (α_weak = fL / dn, α_strong = fH / dn, α = (fL + fH) / dn)
            @test truth.α_weak + truth.α_strong == truth.α
            # the label is reached to ~B: B multiplies the deleterious fixations here (as in
            # rates and analytical_alpha) but not inside the ppos solve
            @test isapprox(truth.α, 0.4; atol = 1e-3) && isapprox(truth.α_weak, 0.2; atol = 1e-3)
            @test isapprox(1 .- (ds / dn) .* (sel ./ neut),
                           MKtest.analytical_alpha(MKtest.copy_param(q))[1]; rtol = 1e-10)
            seed!(9)
            sc_p, sc_d = 2e4 / sum(neut), 5e3 / (dn + ds)
            syn_sfs = hcat(1:(q.nn - 1), MKtest.pois_rand.(sel .* sc_p), MKtest.pois_rand.(neut .* sc_p))
            syn_div = [Float64(MKtest.pois_rand(dn * sc_d)), Float64(MKtest.pois_rand(ds * sc_d)),
                       best["div"][3], best["div"][4]]
            @test only(MKtest.valid_data(p(), [syn_sfs], [syn_div]))   # main: a Bool mask
            post = stages(Float64.(syn_sfs), syn_div, "comp_synthetic")
            # the posterior holds the truth of each component and of the composition
            for (j, tv) in enumerate(truth)
                @test minimum(post[:, j]) <= tv <= maximum(post[:, j])
            end
            # abcmk passes both keywords on: Composed composes exactly with every stat, and
            # the default Marginal is the same posterior summarised column by column
            for st in (MKtest.PosteriorMedian(), MKtest.PosteriorMode(), MKtest.PosteriorMean())
                seed!(10)
                rc = MKtest.abcmk(p(), syn_sfs, syn_div; h5_file = tab, summstat_size = 2000,
                                  output_folder = joinpath(td, "comp_abcmk"), tol = 100,
                                  stat = st, components = MKtest.Composed())[:inference]
                @test rc.α[1] - rc.α_weak[1] - rc.α_strong[1] == 0
                @test rc.ωₐ[1] - rc.ωₐ_weak[1] - rc.ωₐ_strong[1] == 0
                seed!(10)
                rm = MKtest.abcmk(p(), syn_sfs, syn_div; h5_file = tab, summstat_size = 2000,
                                  output_folder = joinpath(td, "comp_abcmk"), tol = 100,
                                  stat = st)[:inference]
                @test rm.α[1] == rc.α[1] && rm.α_weak[1] == rc.α_weak[1]
                @test rm.γ_weak[1] == rc.γ_weak[1]
            end
        end
    end

    # --- summary_abc --------------------------------------------------------------
    @testset "summary_abc" begin
        rng = MersenneTwister(11)
        cn = [:α_weak, :α_strong, :α, :γ_deleterious, :β, :γ_weak, :γ_strong, :B]
        post = [MKtest.DataFrame(rand(rng, 500, 8), cn) for _ in 1:3]
        q(v, p) = MKtest.quantile(v, p)
        # default: PosteriorMedian, Marginal; :stats are the seven marginal statistics
        for (k, P) in ((1, post[1:1]), (3, post), (3, Matrix.(post)))
            r = MKtest.summary_abc(P)
            @test nrow(r[:inference]) == k && nrow(r[:quantiles]) == k
            @test nrow(r[:stats]) == 7k
            @test r[:stats].Stats[1:7] ==
                  ["Min.", "2.5% Perc.", "Median", "Mean", "Mode", "97.5% Perc.", "Max"]
            for i in 1:k, c in cn
                v = post[i][:, c]
                blk = r[:stats][(7i - 6):7i, c]
                @test blk ≈ [minimum(v), q(v, 0.025), MKtest.median(v), MKtest.mean(v),
                             blk[5], q(v, 0.975), maximum(v)]
                @test blk[5] == vec(MKtest.get_mode(post[i]))[findfirst(==(c), cn)]
                @test r[:inference][i, c] == MKtest.median(v)
            end
        end
        # each estimator is its :stats row; PosteriorQuantile(0)/(1) are min/max
        for (st, row) in ((MKtest.PosteriorMedian(), 3), (MKtest.PosteriorMean(), 4),
                          (MKtest.PosteriorMode(), 5), (MKtest.PosteriorQuantile(0.025), 2),
                          (MKtest.PosteriorQuantile(0.975), 6), (MKtest.PosteriorQuantile(0), 1),
                          (MKtest.PosteriorQuantile(1), 7))
            r = MKtest.summary_abc(post; stat = st)
            for i in 1:3
                @test Vector(r[:inference][i, :]) == Vector(r[:stats][7(i - 1) + row, cn])
            end
        end
        # the mode row is get_mode; :quantiles starts with the point estimate
        r = MKtest.summary_abc(post[1:1]; stat = MKtest.PosteriorMode())
        @test vec(Matrix(r[:inference])) == vec(MKtest.get_mode(post[1]))
        @test startswith(r[:quantiles][1, :α], string(round(r[:inference][1, :α], digits = 3)) * " [")
        # the median type and the 0.5 quantile agree to rounding
        @test isapprox(Matrix(MKtest.summary_abc(post; stat = MKtest.PosteriorQuantile(0.5))[:inference]),
                       Matrix(MKtest.summary_abc(post)[:inference]); rtol = 1e-14)
        # Composed: α_strong = α - α_weak in :inference and :quantiles; :stats unchanged
        for st in (MKtest.PosteriorMedian(), MKtest.PosteriorMode(), MKtest.PosteriorQuantile(0.9))
            rm = MKtest.summary_abc(post; stat = st)
            rc = MKtest.summary_abc(post; stat = st, components = MKtest.Composed())
            @test rc[:stats] == rm[:stats]
            @test rc[:inference].α_strong == rc[:inference].α .- rc[:inference].α_weak
            @test rc[:inference].α_strong != rm[:inference].α_strong
            @test rc[:inference][:, Not(:α_strong)] == rm[:inference][:, Not(:α_strong)]
            @test rc[:quantiles][:, Not(:α_strong)] == rm[:quantiles][:, Not(:α_strong)]
            @test startswith(rc[:quantiles][1, :α_strong],
                             string(round(rc[:inference][1, :α_strong], digits = 3)) * " [")
        end
        # the ω block composes the same way; a NaN ω stays NaN and α still composes
        c12 = [:α_weak, :α_strong, :α, :ωₐ_weak, :ωₐ_strong, :ωₐ, :ωₙₐ,
               :γ_deleterious, :β, :γ_weak, :γ_strong, :B]
        p12 = [MKtest.DataFrame(rand(rng, 300, 12), c12)]
        r12 = MKtest.summary_abc(p12; components = MKtest.Composed())[:inference]
        @test r12.ωₐ[1] - r12.ωₐ_weak[1] - r12.ωₐ_strong[1] == 0
        @test r12.ωₙₐ[1] == MKtest.median(p12[1].ωₙₐ)
        p12[1][:, [:ωₐ_weak, :ωₐ_strong, :ωₐ, :ωₙₐ]] .= NaN
        r12 = MKtest.summary_abc(p12; components = MKtest.Composed())[:inference]
        @test isnan(r12.ωₐ_strong[1]) && r12.α[1] - r12.α_weak[1] - r12.α_strong[1] == 0
        # types only: strings are rejected; PosteriorQuantile checks p
        @test_throws TypeError MKtest.summary_abc(post; stat = "median")
        @test_throws TypeError MKtest.summary_abc(post; components = "composed")
        @test_throws ArgumentError MKtest.PosteriorQuantile(1.5)
        @test_throws ArgumentError MKtest.PosteriorQuantile(-0.1)
        @test MKtest.PosteriorQuantile(1).p === 1.0
        # a user-defined estimator needs only a type and a point_estimate method
        @eval struct TrimmedMean10 <: MKtest.AbstractPointEstimate end
        @eval MKtest.point_estimate(::TrimmedMean10, v::Vector{Float64}) =
            MKtest.mean(sort(v)[(length(v) ÷ 10 + 1):(length(v) - length(v) ÷ 10)])
        r = MKtest.summary_abc(post; stat = TrimmedMean10())[:inference]
        @test r[1, :β] == MKtest.mean(sort(post[1].β)[51:450])
    end

    # --- case/control empirical test ----------------------------------------------
    @testset "empirical_test" begin
        cn = [:α_weak, :α_strong, :α]
        ctrl = MKtest.DataFrame(hcat(1.0:99, 1.0:99, 1.0:99), cn)
        mk(x, a) = MKtest.case_control_table(
            MKtest.OrderedDict(:inference => MKtest.DataFrame([x x x], cn)), ctrl; analysis = a)
        tbl = vcat(mk(95.0, "a"), mk(100.0, "b"), mk(50.0, "c"))
        @test nrow(tbl) == 300
        @test names(tbl)[1:3] == ["analysis", "type", "set"]
        r = MKtest.empirical_test(tbl)
        @test r.param[1:3] == ["α", "α_strong", "α_weak"]
        pa = r.p[r.param .== "α"]
        # 95 vs 1..99: 5 controls >= 95 -> 6 / 100; 100: none -> 1 / 100; 50: 50 -> 51 / 100
        @test pa ≈ [0.06, 0.01, 0.51]
        @test r.stars[r.param .== "α"] == ["", "*", ""]
        @test r.control_median[1] == 50.0 && r.n_control[1] == 99
        rl = MKtest.empirical_test(tbl; alternative = :less)
        @test rl.p[rl.param .== "α"] ≈ [96 / 100, 100 / 100, 51 / 100]
        rt = MKtest.empirical_test(tbl; alternative = :two_sided)
        @test rt.p[rt.param .== "α"] ≈ [0.12, 0.02, 1.0]
        @test MKtest._stars.([0.0005, 0.005, 0.03, 0.05]) == ["***", "**", "*", ""]
        # per analysis: the controls of one analysis never enter another's test
        tbl2 = vcat(mk(95.0, "a"), MKtest.case_control_table(
            MKtest.DataFrame([95.0 95.0 95.0], cn), ctrl .* 0 .+ 1000; analysis = "z"))
        @test MKtest.empirical_test(tbl2).p[[1, 4]] ≈ [0.06, 1.0]
        @test_throws ArgumentError MKtest.empirical_test(tbl; omega = true)
        @test_throws ArgumentError MKtest.case_control_table(ctrl, ctrl; analysis = "x")
    end

    # --- RjSampling: analytic properties ----------------------
    @testset "RjSampling" begin
        dac = [4, 5, 10, 20, 50, 200, 661]
        n   = length(dac)
        rng = MersenneTwister(3)
        # a row's TAIL sums: positive, strictly decreasing
        ips = sort(rand(rng, n) .+ 0.1; rev = true) .* 100     # interval λps (any positive)
        ipn = rand(rng, n) .* 80 .+ 1                          # interval λpn
        tail(v) = reverse(cumsum(reverse(v)))
        λps_col, λpn_col = tail(ips), tail(ipn)
        band = [false, false, false, false, false, true, false] # the 200-660 interval
        tps = zeros(n); tpn = zeros(n); r = zeros(n)

        # r ≡ 1: observed Ps per interval proportional to the model's -> independent composition
        s1 = MKtest.RjSampling(0.0, 1.5, 0.15, 0.70, 3.7 .* ips, band)
        for w in (0.0, 0.4, 1.3)
            MKtest._rj_tails!(tps, tpn, r, λps_col, λpn_col, s1, w)
            @test tps ./ (tps .+ tpn) ≈ λps_col ./ (λps_col .+ λpn_col) rtol = 1e-12
        end

        # w = 1: GRAPES form, independent of the anchor
        obs = rand(rng, n) .* 50 .+ 5
        s2  = MKtest.RjSampling(1.0, 1.0, 0.15, 0.70, obs, band)
        MKtest._rj_tails!(tps, tpn, r, λps_col, λpn_col, s2, 1.0)
        @test tpn ./ tps ≈ tail(obs .* ipn ./ ips) ./ tail(obs) rtol = 1e-12

        # anchor: r is normalised so its λps-weighted mean on the band is 1
        s3 = MKtest.RjSampling(0.8, 0.8, 0.15, 0.70, obs, band)
        MKtest._rj_tails!(tps, tpn, r, λps_col, λpn_col, s3, 0.8)
        c = r[6]                      # one interval in the band -> r_band / c == 1
        @test tps[6] - tps[7] ≈ ips[6] * r[6] / c rtol = 1e-12
        @test (tps[6] - tps[7]) / ips[6] ≈ 1.0 rtol = 1e-12

        # _prepare_sampling: fills per-interval Ps from a RAW [DAC pn ps] SFS; others untouched
        raw = hcat(1.0:1321.0, fill(3.0, 1321), fill(2.0, 1321))
        p = MKtest._prepare_sampling(MKtest.RjSampling(), raw, dac, 1322)
        @test p.ps_obs ≈ 2.0 .* [1, 5, 10, 30, 150, 461, 661] rtol = 1e-12
        @test p.band == band
        @test MKtest._prepare_sampling(MKtest.IndependentSampling(), raw, dac, 1322) ===
              MKtest.IndependentSampling()
        @test_throws ArgumentError MKtest._prepare_sampling(
            MKtest.RjSampling(anchor = (0.75, 0.80)), raw, dac, 1322)
        @test_throws ArgumentError MKtest._prepare_sampling(
            MKtest.RjSampling(), raw, reverse(dac), 1322)
        rawz = copy(raw); rawz[200:660, 3] .= 0.0                 # no Ps in one interval
        @test_throws ArgumentError MKtest._prepare_sampling(
            MKtest.RjSampling(), rawz, dac, 1322)
        # the anchor is a FREQUENCY band resolved with the sample size
        dac90 = [4, 5, 10, 20, 50, 200, 400, 661, 800, 925, 1100, 1190]
        @test MKtest._prepare_sampling(MKtest.RjSampling(), raw, dac90, 1322).band ==
              [false, false, false, false, false, true, true, true, true, false, false, false]
        dac50 = [2, 4, 5, 10, 20, 50, 70]
        raw50 = hcat(1.0:99.0, fill(3.0, 99), fill(2.0, 99))
        p50 = MKtest._prepare_sampling(MKtest.RjSampling(), raw50, dac50, 100)
        @test p50.band == [false, false, false, false, true, true, false]
        @test p50.ps_obs ≈ 2.0 .* [2, 1, 5, 10, 30, 20, 30] rtol = 1e-12
        @test_throws ArgumentError MKtest.RjSampling(anchor = (200, 925))  # DAC form
        @test_throws ArgumentError MKtest.RjSampling(w = (1.0, 0.5))
        @test_throws ArgumentError MKtest.poisson_polymorphism!(zeros(n, 2), zeros(n, 2), ones(n),
            ones(n, 2), ones(n, 2), MKtest.RjSampling())    # not prepared

        # sampler means: fixed w, many columns of one row -> E[ps_k] = P_k O_k
        seed!(20260929)
        N = 100_000
        O = 1000.0 .* collect(n:-1:1)
        sps = zeros(n, N); spn = zeros(n, N)
        MKtest.poisson_polymorphism!(spn, sps, O, repeat(λps_col, 1, N), repeat(λpn_col, 1, N), s3)
        MKtest._rj_tails!(tps, tpn, r, λps_col, λpn_col, s3, 0.8)
        P = tps ./ (tps .+ tpn)
        @test all(isapprox.(vec(sum(sps; dims = 2)) ./ N, P .* O; rtol = 2e-3))
        @test all(isapprox.(vec(sum(spn; dims = 2)) ./ N, (1 .- P) .* O; rtol = 2e-3))

        # fixations: RjSampling draws IndependentSampling's, bit for bit
        M = 500
        λds = rand(rng, M) .+ 0.5; λw = rand(rng, M) .* 0.2; λs = rand(rng, M) .* 0.2
        λdn = λw .+ λs .+ rand(rng, M) .+ 0.5
        obs_d = [5000.0, 1.5e6, 0.5e6]
        function fixdraw(smp)
            αm = zeros(M, 3); ωm = zeros(M, 4); dn = zeros(M); ds = zeros(M)
            seed!(99)
            MKtest.poisson_fixation!(αm, ωm, dn, ds, obs_d, λds, λdn, λw, λs, smp)
            return (αm, ωm, dn, ds)
        end
        @test fixdraw(MKtest.RjSampling()) == fixdraw(MKtest.IndependentSampling())

        # tilt: τ = 0 is the untilted sampler; τ ≠ 0 moves Pn mass along ln f, band mass fixed
        pt = MKtest._prepare_sampling(MKtest.RjSampling(tilt = 0.05), raw, dac, 1322)
        @test pt.tilt == 0.05 && pt.band == band
        @test pt.lnf ≈ log.(sqrt.(dac .* [4, 9, 19, 49, 199, 660, 1321]) ./ 1322) rtol = 1e-12
        a0 = zeros(n); b0 = zeros(n); a1 = zeros(n); b1 = zeros(n)
        MKtest._rj_tails!(a0, b0, r, λps_col, λpn_col, pt, 0.7)
        MKtest._rj_tails!(a1, b1, r, λps_col, λpn_col, pt, 0.7, 0.0)
        @test (a0, b0) == (a1, b1)
        MKtest._rj_tails!(a1, b1, r, λps_col, λpn_col, pt, 0.7, 0.3)
        iv(t) = [i < n ? t[i] - t[i+1] : t[n] for i in 1:n]
        @test a1 == a0                                        # Ps untouched
        @test sum(iv(b1)[band]) ≈ sum(iv(b0)[band]) rtol = 1e-12   # no level on the band
        q = iv(b1) ./ iv(b0)
        @test issorted(q)                                     # τ > 0: Pn toward high frequency
        @test q ./ q[1] ≈ exp.(0.3 .* (pt.lnf .- pt.lnf[1])) rtol = 1e-10
        @test_throws ArgumentError MKtest.RjSampling(tilt = -0.1)
        # tilt = 0 draws exactly the untilted sampler (no extra random numbers)
        p0 = MKtest._prepare_sampling(MKtest.RjSampling(0.0, 1.5, 0.15, 0.70, Float64[], Bool[]),
                                      raw, dac, 1322)
        pz = MKtest._prepare_sampling(MKtest.RjSampling(w = (0.0, 1.5), tilt = 0.0), raw, dac, 1322)
        function polydraw(smp)
            a = zeros(n, 50); b = zeros(n, 50)
            seed!(5)
            MKtest.poisson_polymorphism!(a, b, O, repeat(λps_col, 1, 50), repeat(λpn_col, 1, 50), smp)
            return a, b
        end
        @test polydraw(p0) == polydraw(pz)
        @test polydraw(pt) != polydraw(p0)
        # defaults
        d = MKtest.RjSampling()
        @test (d.w_lo, d.w_hi, d.anchor_lo, d.anchor_hi, d.tilt) == (0.0, 2.5, 0.15, 0.70, 0.05)
    end

    # --- FlatAlphaPrior -----------------------------------------------------------
    @testset "FlatAlphaPrior" begin
        rng = MersenneTwister(11)
        M = 200_000
        α = rand(rng, M) .^ 2 .* 0.5                 # realised α piled near 0
        λdn = rand(rng, M) .+ 1.0
        dsdn = hcat(ones(M), λdn, 0.4 .* α .* λdn, 0.6 .* α .* λdn)
        idx = collect(1:M)
        @test MKtest._row_weights(MKtest.TableRowPrior(), dsdn, idx) === nothing
        w = MKtest._row_weights(MKtest.FlatAlphaPrior(floor = 0.0, by_B = false), dsdn, idx)
        mass = [sum(w[floor.(Int, α ./ 0.01) .== k]) for k in 0:49]
        @test maximum(mass) / minimum(mass) < 1.0 + 1e-9  # flat: every bin carries the same weight
        # the floor caps the weight of sparse bins
        wf = MKtest._row_weights(MKtest.FlatAlphaPrior(floor = 0.05, by_B = false), dsdn, idx)
        @test maximum(wf) <= 1 / (0.05 * M) + 1e-15
        seed!(3)
        s = MKtest._sample_rows(MKtest.FlatAlphaPrior(floor = 0.0, by_B = false), dsdn, idx, 400_000)
        h = [count(==(k), floor.(Int, α[s] ./ 0.01)) for k in 0:49] ./ 400_000
        @test maximum(abs.(h .- 0.02)) < 0.003
        # a non-finite α (λdn = 0) is never drawn
        d2 = copy(dsdn); d2[1, 2] = 0.0
        @test MKtest._row_weights(MKtest.FlatAlphaPrior(by_B = false), d2, idx)[1] == 0.0
        # by_B: α flat within each B value, and each B keeps its share of rows
        Bc = ifelse.(1:M .<= M ÷ 4, 400.0, 999.0)               # 25% / 75% of rows
        αB = ifelse.(Bc .== 400.0, rand(rng, M) .* 0.2, 0.2 .+ rand(rng, M) .* 0.3)
        dB = hcat(ones(M), λdn, 0.5 .* αB .* λdn, 0.5 .* αB .* λdn)
        wB = MKtest._row_weights(MKtest.FlatAlphaPrior(floor = 0.0, by_B = true), dB, idx; B = Bc)
        @test sum(wB[Bc .== 400.0]) / sum(wB) ≈ 0.25 rtol = 1e-9
        mB = [sum(wB[(Bc .== 999.0) .& (floor.(Int, αB ./ 0.01) .== k)]) for k in 20:49]
        @test maximum(mB) / minimum(mB) < 1.0 + 1e-9
        # by_B = false ignores B; the default (by_B = true) needs it
        @test MKtest._row_weights(MKtest.FlatAlphaPrior(by_B = false), dsdn, idx) ==
              MKtest._row_weights(MKtest.FlatAlphaPrior(by_B = false), dsdn, idx; B = Bc)
        @test MKtest.FlatAlphaPrior().by_B
        @test_throws ArgumentError MKtest._row_weights(MKtest.FlatAlphaPrior(), dB, idx)
        @test_throws ArgumentError MKtest.FlatAlphaPrior(bin = 0)
        @test_throws ArgumentError MKtest.FlatAlphaPrior(floor = 1.0)
    end

    # --- valid_data ----------------------------------------------------------------
    @testset "valid_data" begin
        dac = [4, 5, 10, 20, 50, 200, 661]
        param = MKtest.parameters(n = 661, dac = dac)
        raw = hcat(1.0:1321.0, fill(3.0, 1321), fill(2.0, 1321))
        gap = copy(raw); gap[200:660, 3] .= 0.0      # no Ps in one interval; every tail > 0
        top = copy(raw); top[661:end, 2] .= 0.0      # no Pn at DAC >= 661
        d = [3000.0, 2000.0, 1e6, 3e5]
        sfs = [raw, gap, top, raw, raw]
        div = [d, d, d, [0.0, 2000.0, 1e6, 3e5], permutedims(d)]   # Dn = 0; 1 x 4 form
        @test MKtest.valid_data(param, sfs, div) == [true, false, false, false, true]
        @test MKtest.valid_data(param, sfs, div; sampling = MKtest.IndependentSampling()) ==
              [true, true, false, false, true]
        # the same intervals as the sampler: what valid_data rejects, RjSampling throws for
        @test MKtest._interval_ps(raw, dac) == MKtest._prepare_sampling(MKtest.RjSampling(), raw, dac, 1322).ps_obs
        @test_throws ArgumentError MKtest._prepare_sampling(MKtest.RjSampling(), gap, dac, 1322)
        @test_throws ArgumentError MKtest.valid_data(param, sfs, div[1:3])
    end

    # --- matched control sets ---------------------------------------------------
    @testset "bootstrap" begin
        rng = MersenneTwister(5)
        od = mktempdir()
        # 3 chromosomes of 200 genes; chr3 has no case gene
        genes = ["G$(lpad(i, 4, '0'))" for i in 1:600]
        chrom = repeat(["1", "2", "3"], inner = 200)
        start = repeat(collect(0:199) .* 10_000, 3)
        bed = joinpath(od, "genes.bed")
        open(bed, "w") do io
            for i in 1:600
                println(io, join((chrom[i], start[i], start[i] + 3000, genes[i], "+"), '\t'))
            end
        end
        case = genes[vcat(5:5:100, 205:5:300)]               # 40 case genes on chr1 and chr2
        cfile = joinpath(od, "case.txt")
        open(io -> foreach(g -> println(io, g), case), cfile, "w")
        ff = joinpath(od, "factors.txt")
        open(ff, "w") do io
            for g in genes
                println(io, g, " ", round(10 + randn(rng), digits = 4), " ", round(-5 + randn(rng), digits = 4))
            end
        end

        # get_distance against brute force over the midpoints
        dist = MKtest.get_distance(cfile, bed)
        mid = div.(2 .* start .+ 3000, 2)
        brute = map(1:600) do i
            same = [mid[j] for j in 1:600 if genes[j] in case && chrom[j] == chrom[i]]
            isempty(same) ? typemax(Int) : minimum(abs.(same .- mid[i]))
        end
        @test dist.gene == genes && dist.distance == brute
        @test all(dist.distance[401:600] .== typemax(Int))   # chr3 kept, not skipped

        prefix = joinpath(od, "out")
        bp = MKtest.bootstrap_parameters(data = cfile, annotations = bed, confounding_factors = ff,
                                         output = prefix, distance = 20_000, iterations = 2,
                                         tolerance = 0.05, filter_hla_hist = false, seed = 1)
        res = MKtest.bootstrap(bp)
        lines = readlines(prefix * "_control.txt")
        @test lines == res.control && length(lines) == 200
        @test [split(l)[1] for l in lines] == ["sample_$i" for i in 1:200]   # unique labels
        sets = [split(l)[2:end] for l in lines]
        @test all(length.(sets) .== 40)
        far = Set(dist.gene[dist.distance .>= 20_000])
        @test all(g -> g in far && !(g in case), reduce(vcat, sets))
        @test any(g -> g in genes[401:600], reduce(vcat, sets))           # chr3 genes are used
        @test Set(readlines(prefix * "_case.txt")) == Set(case)
        # the pooled controls match the case averages within the tolerance
        fac = Dict(r[1] => parse.(Float64, r[2:3]) for r in split.(readlines(ff)))
        case_avg = sum(fac[g] for g in case) ./ 40
        pool_avg = sum(fac[g] for g in reduce(vcat, sets)) ./ (200 * 40)
        @test all(abs.(pool_avg .- case_avg) .<= 0.05 .* abs.(case_avg))
        # a seed reproduces the sets
        @test MKtest.bootstrap(bp).control == lines
        # an impossible tolerance fails loudly
        fac_bad = Dict(g => (g in case ? [100.0, 1.0] : [1.0, 1.0]) for g in genes)
        @test isempty(MKtest.get_bootstrap(case, setdiff(genes, case); factors = fac_bad,
                                           max_reps = 3, tolerance = 0.01, rng = MersenneTwister(1)))
    end

    # --- small methods: NI_TG, alpha_b, the Julia ABC -------------------------
    @testset "ni_cmh" begin
        # hand-computed NI_TG (Stoletzki & Eyre-Walker 2011): sum(Ds Pn / (Ps + Ds)) /
        # sum(Dn Ps / (Ps + Ds)) over genes
        sfs = [[1.0 4.0 6.0], [1.0 2.0 2.0]]
        div = [[3.0, 5.0], [1.0, 4.0]]
        out = MKtest.ni_cmh(sfs, div)
        @test out isa NamedTuple
        @test out.ni ≈ (5 * 4 / 11 + 4 * 2 / 6) / (3 * 6 / 11 + 1 * 2 / 6)
        # an Int zero in the divergence is caught like a Float one
        @test isnan(MKtest.ni([1.0 4.0 6.0], [0, 5]).ni)
        # 39 bins = 40 chromosomes: 20 individuals is the most it can be projected to
        sfs40 = hcat(1:39, ones(39), ones(39))
        @test size(MKtest.project(sfs40, 20), 1) == 39
        @test_throws AssertionError MKtest.project(sfs40, 30)
    end

    @testset "α_b" begin
        seed!(3)
        r_s = [[40, 100, 60, 100] .+ i for i in 1:12]
        out = MKtest.α_b(r_s)
        @test nrow(out) == 12              # every dataset, not the first 10
        @test all(out.α_low .<= out.αᵦ .<= out.α_high)
        # 1 - 1/Z with Z = (sN/sS) / (rN/rS)
        @test out.αᵦ[1] ≈ 1 - 1 / ((61 / 101) / (41 / 101))
    end

    @testset "abc_loclinear transformations" begin
        # oracle: the tangent transformation and its inverse are a round trip
        for x in (0.1, 0.5, 2.0)
            @test MKtest.undo_tangent_transfromation(MKtest.tangent_transfromation(x, 0.1, 2.0), 0.1, 2.0) ≈ x
        end
        # synthetic table: 2 positive parameters, 3 statistics linear in them
        rng = MersenneTwister(11)
        nsim = 5000
        θ = hcat(rand(rng, nsim) .* 0.8 .+ 0.1, rand(rng, nsim) .* 9 .+ 1)
        stat(t, e) = [t[1] + 0.1 * t[2], 2 * t[1] - 0.05 * t[2], 0.3 * t[2]] .+ e
        summ = vcat([vcat(θ[i, :], stat(θ[i, :], 0.02 .* randn(rng, 3)))' for i in 1:nsim]...)
        truth = [0.5, 5.0]
        od = mktempdir()
        tfile, sfile = joinpath(od, "target.tsv"), joinpath(od, "summ.tsv")
        writetab(f, m) = open(io -> foreach(r -> println(io, join(r, '\t')), eachrow(m)), f, "w")
        writetab(tfile, permutedims(stat(truth, zeros(3))))
        writetab(sfile, summ)
        for tr in ("none", "log", "tan"), k in ("epanechnikov", "gaussian")
            acc, adj = MKtest.abc_loclinear(tfile, sfile; P = 2, tol = 250,
                                            transformation = tr, kernel = k)
            @test size(acc) == (250, 2) && all(isfinite, adj)
            @test all(abs.(vec(sum(adj, dims = 1)) ./ 250 .- truth) .< [0.05, 0.5])
            tr == "log" && @test all(>(0), adj)
            if tr == "tan"
                @test all(minimum(θ, dims = 1) .- 1e-3 .<= adj .<= maximum(θ, dims = 1) .+ 1e-3)
            end
        end
        neg = copy(summ); neg[:, 1] .-= 1.0
        nfile = joinpath(od, "neg.tsv")
        writetab(nfile, neg)
        @test_throws ArgumentError MKtest.abc_loclinear(tfile, nfile; P = 2, tol = 250,
                                                        transformation = "log")
    end

    include("abc_in_memory.jl")

end

# dfe_ml engine (src/ml/)
include("ml/engine.jl")
include("ml/extension.jl")
include("ml/grapes_parity.jl")

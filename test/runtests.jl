using MKtest
using Test
using DataFrames: nrow, names, propertynames
using JLD2
using Random: seed!, MersenneTwister
using DelimitedFiles: readdlm

@testset "MKtest.jl" begin

    # --- ABC forward model: parameter solve ---------------------------------
    # `≈` rather than `==`: `set_θ!` is a closed-form solve (it was `Roots.find_zero`),
    # and exact equality tested the arithmetic path rather than the answer.
    @testset "parameters" begin
        adap = MKtest.parameters(N=1000, n=661, gam_dfe=-457, gam_flanking=-1000,
                                 al_low=0.2, B=0.999)
        binom = MKtest.binom_op(adap)
        MKtest.set_θ!(adap)
        MKtest.set_ppos!(adap)

        # 0.00020733640485526956 with the previous Φ (no factor 2 on s). Regression
        # value, not an independent check.
        @test adap.ppos_h ≈ 0.00020731419036651257
        @test adap.θ_flanking ≈ 3.001501000750676e-6
    end

    # --- rho prior in rates() ------------------------------------------------
    # `ρ` reaches the model through `set_θ!` (the flanking mutation rate that produces
    # the target B) and `Φ`. `rates(...; rho = [lo, hi])` draws one value per model and
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
                             rho=rho)
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
        p.gam_flanking = Int(m.gam_flanking[1])
        @test MKtest.set_θ!(p) ≈
              -log(p.B) * (p.Lf * m.rho[1] - m.gam_flanking[1]) / (2.0 * p.Lf)
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
            q.gH, q.gL = Int(m.gH[i]), Int(m.gL[i])
            q.gam_flanking, q.gam_dfe = Int(m.gam_flanking[i]), Int(m.gam_dfe[i])
            q.shape = m.shape[i]; q.scale = abs(q.shape / q.gam_dfe)
            q.al_tot, q.al_low, q.ρ = m.al_tot[i], m.al_low[i], m.rho[i]
            q.B = m.B[i] / 1000; MKtest.set_θ!(q)
            # ppos recovered EXACTLY from the stored fixation rates:
            # f_pos = 0.75 * (ppos/2) * Φ * p_fix
            pp(f, g) = f == 0 ? 0.0 : 2f / (0.75 * MKtest.Φ(q, g) * MKtest.p_fix(q, g))
            q.ppos_l = pp(r[:dsdn][i, 3], q.gL)
            q.ppos_h = pp(r[:dsdn][i, 4], q.gH)
            bn = MKtest.binom_op(q)
            neut = MKtest.cumulative_vector!(MKtest.sfs_neut(q, bn))
            sel  = MKtest.cumulative_vector!(MKtest.sfs_pos(q, q.gH, q.ppos_h, bn) .+
                                             MKtest.sfs_pos(q, q.gL, q.ppos_l, bn) .+
                                             MKtest.sfs_neg(q, q.ppos_l + q.ppos_h, bn))
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
        pa = mkp(); pa.al_tot = 0.4; pa.al_low = 0.2; MKtest.set_θ!(pa)
        MKtest.set_ppos!(pa)
        @test maximum(MKtest._ppos_residual(pa, pa.ppos_l, pa.ppos_h)) < 1e-6
        @test pa.ppos_l > 0 && pa.ppos_h > 0

        # al_low == 0 means ppos_l == 0 IS the root, not a failure.
        pz = mkp(); pz.al_tot = 0.4; pz.al_low = 0.0; MKtest.set_θ!(pz)
        MKtest.set_ppos!(pz)
        @test pz.ppos_l < 1e-20
        @test maximum(MKtest._ppos_residual(pz, pz.ppos_l, pz.ppos_h)) < 1e-6

        # An INFEASIBLE target is refused, not written: a very peaked deleterious DFE
        # produces almost no deleterious fixations, so a low al_tot cannot be reached.
        pj = mkp(); pj.shape = 7.5; pj.gam_dfe = -500
        pj.scale = abs(pj.shape / pj.gam_dfe)
        pj.al_tot = 0.06; pj.al_low = 0.054; MKtest.set_θ!(pj)
        @test_throws ErrorException MKtest.set_ppos!(pj)

        # the model record is written next to the table
        seed!(1)
        f = joinpath(mktempdir(), "prov.jld2")
        MKtest.rates(mkp(); iterations = 8, output_file = f, rho = 0.001)
        o = MKtest._read_rates(f, mkp()).options
        @test o == MKtest._rates_options()
        @test startswith(o, "model=best-2026-09-28")
    end

    # --- the forward model: what each choice is -------------------------------
    @testset "best-model physics" begin
        p = MKtest.parameters(N = 7310, n = 661, gam_dfe = -457,
                              gam_flanking = -457, Lf = 200_000,
                              gL = 10, gH = 500, ρ = 0.001, B = 0.4,
                              dac = [2, 4, 5, 10, 20, 50, 200, 661, 925])
        MKtest.set_θ!(p)

        # Φ = B^κ follows from set_θ!'s closed form; the published expression carries
        # the factor of two on the beneficial coefficient (κ_L/κ_H = 8.2, not 3.8).
        κ(g) = log(MKtest.Φ(p, g)) / log(p.B)
        @test isapprox(κ(10), 0.5455; atol = 1e-3)
        @test isapprox(κ(500), 0.0662; atol = 1e-3)

        # sfs_pos is the paper's "replace N with Ne in the frequency spectra": a literal
        # transcription -- level B, sojourn exp-form at γ·B, Φ flat.
        q = MKtest.parameters(N = 1000, n = 661, gam_dfe = -457, gam_flanking = -500,
                              Lf = 200_000, gL = 5, gH = 500, ρ = 0.001, B = 0.4,
                              dac = [4, 5, 10, 20, 50, 200, 661, 925])
        q.ppos_l, q.ppos_h = 0.0077070856, 0.00015930139
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
        for (s, pp) in ((q.gL, q.ppos_l), (q.gH, q.ppos_h))
            @test isapprox(MKtest.sfs_pos(q, s, pp, bn), sfs_pos_ref(q, s, pp, bn); rtol = 1e-10)
        end

        # The B level makes beneficial polymorphism relative to neutral FALL as BGS
        # strengthens (without it the ratio rose, 1.486; SLiM measures 0.637).
        function sel_over_neut(B)
            r = MKtest.copy_param(q); r.B = B; MKtest.set_θ!(r)
            b = MKtest.binom_op(r)
            nt = MKtest.cumulative_vector!(MKtest.sfs_neut(r, b))
            sl = MKtest.cumulative_vector!(MKtest.sfs_pos(r, r.gL, r.ppos_l, b))
            sl[925] / nt[925]
        end
        ratio = sel_over_neut(0.4) / sel_over_neut(0.999)
        @test 0.3 < ratio < 1.0

        # γ·B against a sojourn at γ: large at low B, vanishing as B -> 1.
        function src_rel(B)
            r = MKtest.copy_param(q); r.B = B
            NN2 = ceil(Int, r.NN * B); sc = zeros(NN2 + 1); sc2 = zeros(NN2 + 1)
            a = MKtest._src_pos!(sc,  r, r.gL * B, r.ppos_l, NN2 + 1)
            b = MKtest._src_pos!(sc2, r, r.gL,     r.ppos_l, NN2 + 1)
            maximum(abs.(a[2:end-1] .- b[2:end-1]) ./ abs.(b[2:end-1]))
        end
        @test src_rel(0.4) > 0.05
        @test src_rel(0.999) < 1e-2
    end

    # --- syn_nonsyn_ratio reaches the table ------------------------------------
    # `rates()` solves each task on `copy_param(param)`, which lists the fields one by
    # one: a field it misses silently falls back to its default.
    @testset "syn_nonsyn_ratio" begin
        mkp(x) = MKtest.parameters(N = 500, n = 50, dac = [2, 4, 5, 10, 20, 50],
                                   syn_nonsyn_ratio = x)
        @test MKtest.copy_param(mkp(0.3)).syn_nonsyn_ratio == 0.3
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
                           iterations = 10, output_file = joinpath(od, "s1.jld2"), rho = 0.001)
        @test minimum(out[:models].shape) >= lo - 1e-12
        @test maximum(out[:models].shape) <= hi + 1e-12
        seed!(7)
        out2 = MKtest.rates(MKtest.parameters(N = 500, n = 50, shape = 0.4,
                                              dac = [2, 4, 5, 10, 20, 50]);
                            iterations = 10, output_file = joinpath(od, "s2.jld2"), rho = 0.001)
        @test maximum(out2[:models].shape) > hi
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
        MKtest.rates(mkparam(); iterations = 20, rho = 0.001, output_file = tab)
        r = MKtest._read_rates(tab, mkparam())
        check_table(r, fx, DAC)

        for (B, key) in ((0.4, "aa_04"), (0.999, "aa_0999"))
            p = MKtest.parameters(N = 500, n = 50, dac = DAC, B = B, al_tot = 0.4,
                                  al_low = 0.2, gam_dfe = -457, gam_flanking = -500,
                                  gL = 5, gH = 500, ρ = 0.001)
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
        MKtest.rates(mkparam(); iterations = 20, rho = [0.0001, 0.004], output_file = tab)
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
                               row_prior = MKtest.TableRowPrior(), stat = "mode")
            @test issorted(res[:inference].α)
        end

        # a single set as parse_sfs returns it (Matrix, Vector) equals the wrapped call
        one(s, d) = (seed!(3); MKtest.abcmk(p(), s, d; h5_file = tab, summstat_size = 2000,
                                            output_folder = joinpath(td, "one"), tol = 100))
        @test one(sfs[1], div[1])[:inference] == one([sfs[1]], [div[1]])[:inference]
    end

    # --- summary_abc --------------------------------------------------------------
    @testset "summary_abc" begin
        rng = MersenneTwister(11)
        cn = [:α_weak, :α_strong, :α, :γ₋, :β, :γ₊, :γ₊₊, :B]
        post = [MKtest.DataFrame(rand(rng, 500, 8), cn) for _ in 1:3]
        q(v, p) = MKtest.quantile(v, p)
        for (k, P) in ((1, post[1:1]), (3, post), (3, Matrix.(post)))
            r = MKtest.summary_abc(P; stat = "median")
            @test nrow(r[:inference]) == k && nrow(r[:quantiles]) == k
            @test nrow(r[:stats]) == 7k
            @test r[:stats].Stats[1:7] ==
                  ["Min.", "2.5% Perc.", "Median", "Mean", "Mode", "97.5% Perc.", "Max"]
            for i in 1:k, c in cn
                v = post[i][:, c]
                blk = r[:stats][(7i - 6):7i, c]
                @test blk ≈ [minimum(v), q(v, 0.025), MKtest.median(v), MKtest.mean(v),
                             blk[5], q(v, 0.975), maximum(v)]
                @test r[:inference][i, c] == MKtest.median(v)
            end
        end
        # the mode row is get_mode, the inference row follows `stat`
        r = MKtest.summary_abc(post[1:1]; stat = "Mode")
        @test vec(Matrix(r[:inference])) == vec(MKtest.get_mode(post[1]))
        @test startswith(r[:quantiles][1, :α], string(round(r[:inference][1, :α], digits = 3)) * " [")
        @test_throws ArgumentError MKtest.summary_abc(post; stat = "foo")
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

end

# The ABC-MK forward model and its reference table: priors, the tail-weight projection,
# the ppos solver, the physics of the chosen model, and the two frozen fixtures.
# Group "rates" of runtests.jl.

@testset "ABC-MK rates" begin

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
        m = r.models
        @test :ρ in propertynames(m)
        @test names(m)[9] == "ρ"            # appended last: cols 1-8 keep their meaning
        @test all(1e-4 .<= m.ρ .<= 5e-3)
        nmodels = size(m, 1) ÷ nb
        @test all(allequal(m.ρ[(i-1)*nb+1:i*nb]) for i in 1:nmodels)
        @test length(unique(m.ρ)) == nmodels
        @test names(m) == string.(MKtest._MODEL_COLUMNS) && !("ne_ratio" in names(m))
        @test names(MKtest._read_rates(f, param).models)[9] == "ρ"
        # the return value is the table load_rates reads back, field by field
        @test r isa MKtest.RatesTable
        t = MKtest.load_rates(f, param)
        @test (r.models, r.neut, r.sel, r.dsdn, r.dac, r.options, r.file, r.N, r.n, r.cutoff) ==
              (t.models, t.neut, t.sel, t.dsdn, t.dac, t.options, t.file, t.N, t.n, t.cutoff)
        @test typeof(r.models) == typeof(t.models) && typeof(r.neut) == typeof(t.neut)
        # on disk the block holds the same fields as plain arrays, and the block is unique
        grp = MKtest._jld2_group_key(param.N, param.n, param.cutoff[1], param.cutoff[end])
        jldopen(f, "r") do h
            g = h[grp]
            @test g["models"] isa Matrix{Float64} && g["models"] == Matrix(r.models)
            @test g["model_names"] == names(r.models)
            @test g["neut"] isa Matrix{Float64} && size(g["neut"]) == (size(r.dsdn, 1), length(r.dac))
            @test g["neut"][:, 2] == r.neut[r.dac[2]] && g["sel"][:, end] == r.sel[r.dac[end]]
            @test g["dsdn"] == r.dsdn && g["dac"] == r.dac && g["options"] == r.options
        end
        @test_throws ArgumentError MKtest._write_rates(f, r)
        @test_throws ArgumentError MKtest.load_rates(f, MKtest.parameters(N = 500, n = 10, dac = [2, 4]))
        # a file of the layout before 2026-10-22 (DataFrame with the old column names and
        # ne_ratio, OrderedDict{Int,Array} tails) reads to the same table
        f_old = joinpath(od, "old_layout.jld2")
        old_names = ["B", "al_low", "al_tot", "gam_flanking", "gL", "gH", "shape", "gam_dfe", "rho"]
        jldopen(f_old, "w") do h
            h[grp * "/models"] = hcat(MKtest.rename(r.models, old_names),
                                      MKtest.DataFrame(ne_ratio = ones(size(r.models, 1))))
            h[grp * "/neut"]   = MKtest.OrderedDict{Int,Array}(k => r.neut[k] for k in r.dac)
            h[grp * "/sel"]    = MKtest.OrderedDict{Int,Array}(k => r.sel[k] for k in r.dac)
            h[grp * "/dsdn"]   = r.dsdn
            h[grp * "/dac"]    = r.dac
            h[grp * "/options"] = r.options
        end
        t_old = MKtest.load_rates(f_old, param)
        @test (t_old.models, t_old.neut, t_old.sel, t_old.dsdn, t_old.dac, t_old.options) ==
              (r.models, r.neut, r.sel, r.dsdn, r.dac, r.options)
        @test typeof(t_old.neut) == typeof(r.neut)

        # fixed forms: a Real, a one-value range and nothing fix ρ and draw no RNG; one of
        # them through rates() to the stored column
        @test MKtest._rho_prior(0.002, param) === (nothing, 0.002)
        @test MKtest._rho_prior([0.002], param) === (nothing, 0.002)
        @test MKtest._rho_prior(nothing, param) === (nothing, param.ρ)
        @test all(==(0.002), run_rates(0.002; iterations=10)[1].models.ρ)

        # invalid ranges must throw rather than silently pick something
        for bad in ([-1.0, 1.0], [5e-3, 1e-4], -0.001, Float64[])
            @test_throws ArgumentError run_rates(bad; iterations=5)
        end

        # the recorded rho is the one set_θ! actually used
        p = MKtest.copy_param(param)
        p.B = m.B[1] / 1000
        p.ρ = m.ρ[1]
        p.γ_flanking = Int(m.γ_flanking[1])
        @test MKtest.set_θ!(p) ≈
              -log(p.B) * (p.L_flanking * m.ρ[1] - m.γ_flanking[1]) / (2.0 * p.L_flanking)
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
        m = r.models
        # every stored row equals the full-spectrum model at the same parameters
        for i in (1, 2, 3, nrow(m))
            q = MKtest.copy_param(param)
            q.γ_strong, q.γ_weak = Int(m.γ_strong[i]), Int(m.γ_weak[i])
            q.γ_flanking, q.γ_deleterious = Int(m.γ_flanking[i]), Int(m.γ_deleterious[i])
            q.shape = m.shape[i]; q.scale = abs(q.shape / q.γ_deleterious)
            q.α, q.α_weak, q.ρ = m.α[i], m.α_weak[i], m.ρ[i]
            q.B = m.B[i] / 1000; MKtest.set_θ!(q)
            # ppos recovered EXACTLY from the stored fixation rates:
            # f_pos = 0.75 * (ppos/2) * Φ * p_fix
            pp(f, g) = f == 0 ? 0.0 : 2f / (0.75 * MKtest.Φ(q, g) * MKtest.p_fix(q, g))
            q.p_weak = pp(r.dsdn[i, 3], q.γ_weak)
            q.p_strong = pp(r.dsdn[i, 4], q.γ_strong)
            bn = MKtest.binom_op(q)
            neut = MKtest.cumulative_vector!(MKtest.sfs_neut(q, bn))
            sel  = MKtest.cumulative_vector!(MKtest.sfs_pos(q, q.γ_strong, q.p_strong, bn) .+
                                             MKtest.sfs_pos(q, q.γ_weak, q.p_weak, bn) .+
                                             MKtest.sfs_neg(q, q.p_weak + q.p_strong, bn))
            @test reldiff([r.neut[k][i] for k in param.dac], neut[param.dac]) < 1e-12
            @test reldiff([r.sel[k][i]  for k in param.dac], sel[param.dac])  < 1e-12
        end

        # The tail-weight cache is keyed on the DAC vector as well as B.
        p3 = MKtest.parameters(N=500, n=100, dac=[4, 20, 50], cutoff=[0.0, 1.0], B=0.999)
        seed!(5)
        r6 = MKtest.rates(param; iterations=30, output_file=joinpath(od, "k6.jld2"),
                          B_bins=B_bins)
        seed!(5)
        r3 = MKtest.rates(p3; iterations=30, output_file=joinpath(od, "k3.jld2"),
                          B_bins=B_bins)
        for k in keys(r3.neut)
            @test reldiff(r3.neut[k], r6.neut[k]) < 1e-12
            @test reldiff(r3.sel[k],  r6.sel[k])  < 1e-12
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

        # α_weak == 0 means p_weak == 0 IS the root, not a failure.
        pz = mkp(); pz.α = 0.4; pz.α_weak = 0.0; MKtest.set_θ!(pz)
        MKtest.set_ppos!(pz)
        @test pz.p_weak < 1e-20
        @test maximum(MKtest._ppos_residual(pz, pz.p_weak, pz.p_strong)) < 1e-6

        # A very peaked deleterious DFE fixes almost nothing, so a low α needs a tiny
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
        @test b.models[:, 1:8] == a.models[:, 1:8]
        # the neutral SFS is linear in the ratio
        @test all(isapprox(b.neut[k], a.neut[k] .* (0.3 / 0.25); rtol = 1e-12) for k in [2, 50])
    end

    # --- the DFE shape prior is derived from param.shape ----------------------
    @testset "DFE shape prior is param-derived" begin
        od = mktempdir()
        lo, hi = 0.184 * 2.0^-2, 0.184 * 2.0^2          # [0.046, 0.736]
        seed!(7)
        out = MKtest.rates(MKtest.parameters(N = 500, n = 50, shape = 0.184,
                                             dac = [2, 4, 5, 10, 20, 50]);
                           iterations = 10, output_file = joinpath(od, "s1.jld2"), ρ = 0.001)
        @test minimum(out.models.shape) >= lo - 1e-12
        @test maximum(out.models.shape) <= hi + 1e-12
        seed!(7)
        out2 = MKtest.rates(MKtest.parameters(N = 500, n = 50, shape = 0.4,
                                              dac = [2, 4, 5, 10, 20, 50]);
                            iterations = 10, output_file = joinpath(od, "s2.jld2"), ρ = 0.001)
        @test maximum(out2.models.shape) > hi
    end

    @testset "shape_factor widens the DFE shape prior" begin
        # f = 2 is the published window exactly (powers of two), so old tables keep bits.
        u2 = MKtest._shape_prior(0.184, 2)
        @test (u2.a, u2.b) === (0.184 / 4, 0.184 * 4)
        u3 = MKtest._shape_prior(0.184, 3)
        @test (u3.a, u3.b) === (0.184 / 8, 0.184 * 8)
        @test MKtest._shape_prior(0.184, 2.5).b ≈ 0.184 * 2^2.5
        @test_throws ArgumentError MKtest._shape_prior(0.184, 0)
        @test_throws ArgumentError MKtest._shape_prior(0.184, -1)
        @test_throws ArgumentError MKtest._shape_prior(0.184, Inf)
        od = mktempdir()
        mkp() = MKtest.parameters(N = 500, n = 50, shape = 0.184, dac = [2, 4, 5, 10, 20, 50])
        seed!(11)
        d = MKtest.rates(mkp(); iterations = 10, output_file = joinpath(od, "d.jld2"), ρ = 0.001)
        seed!(11)
        d2 = MKtest.rates(mkp(); iterations = 10, output_file = joinpath(od, "d2.jld2"),
                          ρ = 0.001, shape_factor = 2)
        @test d.sel == d2.sel && d.dsdn == d2.dsdn
        # f = 4: [0.0115, 2.944]; draws above shape 1 are solved, not refused
        seed!(11)
        w = MKtest.rates(mkp(); iterations = 40, output_file = joinpath(od, "w.jld2"),
                         ρ = 0.001, shape_factor = 4)
        sh = w.models.shape
        @test minimum(sh) >= 0.184 / 16 - 1e-12 && maximum(sh) <= 0.184 * 16 + 1e-12
        @test count(>(1), sh) >= 10
        @test size(w.models, 1) == 40 * 37               # nothing refused
    end

    # --- reference outputs ----------------------------------------------------
    # The fixtures were produced by the code (commit 8c68fa35) that RjSampling was
    # tested with: abcmk_best_model.jld2 for the forward model and the default sampler,
    # abcmk_des.jld2 for a ρ-drawn table and RjSampling (make_des_fixture.jl). Both must
    # be reproduced at one thread.
    function check_table(r, fx, DAC)
        # the fixtures hold the pre-2026-10-22 column names and a constant ne_ratio column
        fxnames = [get(MKtest._LEGACY_MODEL_NAMES, n, n) for n in fx["model_names"]]
        keep = findall(!=("ne_ratio"), fxnames)
        @test all(==(1.0), fx["models"][:, findfirst(==("ne_ratio"), fxnames)])
        @test Matrix(r.models) == fx["models"][:, keep]
        @test names(r.models) == fxnames[keep]
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
        fx = BEST
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

        # the defaults are RjSampling(w = (0.5, 1.5), tilt = 0.05) with FlatAlphaPrior(by_B = true)
        summ(; kw...) = begin
            of = mktempdir()
            seed!(7)
            MKtest.summary_statistics(mkparam(), [fx["sfs"]], [fx["div"]]; h5_file = tab,
                                      summstat_size = 200, output_folder = of, kw...)
            MKtest._with_gz(readdlm, joinpath(of, "summstat_1.txt.gz"))
        end
        @test summ() == summ(sampling = MKtest.RjSampling(w = (0.5, 1.5), tilt = 0.05),
                             row_prior = MKtest.FlatAlphaPrior(by_B = true))
        @test summ() != fx["summ_indep"]
    end

    @testset "RjSampling fixture" begin
        fx  = JLD2.load(joinpath(@__DIR__, "fixtures", "abcmk_des.jld2"))
        best = BEST
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
end

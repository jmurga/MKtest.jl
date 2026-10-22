# The Poisson samplers, the row priors, valid_data, the ABC (file and in-memory routes),
# the posterior summaries and the case/control test. Group "abc" of runtests.jl.

# One 300-draw table (seed 1) for every testset that needs a table: the same draws each of
# them built for itself before 2026-10-22.
seed!(1)
const TAB300 = quiet(() -> MKtest.rates(P500(); iterations = 300,
                                        output_file = joinpath(mktempdir(), "t300.jld2")))

@testset "ABC-MK inference" begin

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
        best = BEST
        p = P500
        tab = TAB300.file
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
            @test only(MKtest.valid_mask(p(), [syn_sfs], [syn_div]))   # main: a Bool mask
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
                                  tol = 100, in_memory = true,
                                  stat = st, components = MKtest.Composed())[:inference]
                @test rc.α[1] - rc.α_weak[1] - rc.α_strong[1] == 0
                @test rc.ωₐ[1] - rc.ωₐ_weak[1] - rc.ωₐ_strong[1] == 0
                seed!(10)
                rm = MKtest.abcmk(p(), syn_sfs, syn_div; h5_file = tab, summstat_size = 2000,
                                  tol = 100, in_memory = true,
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
        @test (d.w_lo, d.w_hi, d.anchor_lo, d.anchor_hi, d.tilt) == (0.5, 1.5, 0.15, 0.70, 0.05)
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
        @test MKtest.valid_mask(param, sfs, div) == [true, false, false, false, true]
        @test MKtest.valid_mask(param, sfs, div; sampling = MKtest.IndependentSampling()) ==
              [true, true, false, false, true]
        # valid_data applies the mask and returns the pair to pass on
        s2, d2 = @test_logs (:warn, r"dropping 3 of 5") MKtest.valid_data(param, sfs, div)
        @test s2 == [raw, raw] && d2 == [d, permutedims(d)]
        @test MKtest.valid_data(param, sfs[[1, 5]], div[[1, 5]]) == (sfs[[1, 5]], div[[1, 5]])
        # the same intervals as the sampler: what valid_data rejects, RjSampling throws for
        @test MKtest._interval_ps(raw, dac) == MKtest._prepare_sampling(MKtest.RjSampling(), raw, dac, 1322).ps_obs
        @test_throws ArgumentError MKtest._prepare_sampling(MKtest.RjSampling(), gap, dac, 1322)
        @test_throws ArgumentError MKtest.valid_data(param, sfs, div[1:3])
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

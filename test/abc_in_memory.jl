# In-memory ABC (`abcmk(...; in_memory = true)`) and the Julia port of ABCreg.
#
# Oracles: the ABCreg binary (scripts/reg) and the file `abc` route, on the very same
# summaries. Thresholds fixed before running: in-memory "abc" equals the file route
# exactly; `abcreg_loclinear` equals the binary to its printed precision (6 significant
# digits, so |a - b| <= 1e-5 |b| + 1e-9) with the same accepted rows.
#
# Included by test/abcmk_abc.jl (group "abc": `Pkg.test(test_args = ["abc"])`), which
# provides `quiet`, `BEST`, `P500` and the shared 300-draw table `TAB300`.

@testset "ABCreg in memory" begin
    reg = joinpath(pkgdir(MKtest), "scripts", "reg")
    have_reg = Sys.islinux() && isfile(reg)
    close6(a, b) = size(a) == size(b) && all(abs.(a .- b) .<= 1e-5 .* abs.(b) .+ 1e-9)

    # -- direct: an all-positive fixture, so ABCreg's -L is defined everywhere ------
    rng = MersenneTwister(20261018)
    Nf, S = 2000, 6
    par = hcat(rand(rng, Nf, 3) .* 0.8 .+ 0.05, rand(rng, Nf) .* 900 .+ 100,
               ceil.(rand(rng, Nf) .* 9), rand(rng, Nf) .* 0.8 .+ 0.2)
    ss = par[:, 3] .+ par[:, 6] .* 0.2 .+ randn(rng, Nf, S) .* 0.05
    target = ss[7, :] .+ 0.01
    tl = 150
    @test size(MKtest.abcreg_loclinear(target, par, ss; tol = tl)) == (tl, 6)
    @test_throws ArgumentError MKtest.abcreg_loclinear(target, par, ss; tol = Nf + 1)
    @test_throws ArgumentError MKtest.abcreg_loclinear(target, par, ss; tol = tl, transformation = "sqrt")
    @test_throws DimensionMismatch MKtest.abcreg_loclinear(target[1:end-1], par, ss; tol = tl)
    if have_reg
        td = mktempdir()
        writedlm(joinpath(td, "s.txt"), hcat(par, ss), '\t')
        writedlm(joinpath(td, "a.txt"), reshape(target, 1, S), '\t')
        for (t, flag, suffix) in (("none", String[], ""), ("log", ["-L"], ".log"),
                                  ("tan", ["-T"], ".tangent"))
            out = joinpath(td, "o_$t")
            run(pipeline(`$reg -d $(joinpath(td, "a.txt")) -p $(joinpath(td, "s.txt")) -P 6 -S $S -t $(tl / Nf) -b $out $flag`,
                         stdout = devnull, stderr = devnull))
            bin = quiet(() -> MKtest.read_gz(out * ".0" * suffix * ".post.gz"))
            jl = MKtest.abcreg_loclinear(target, par, ss; tol = tl, transformation = t)
            @test close6(jl, bin)
        end
        # -L with a negative parameter: the binary aborts; the port gives NaN in that
        # column only, and the other columns equal the binary run without the negatives
        parn = copy(par); parn[1:40, 2] .= -parn[1:40, 2]
        jl = MKtest.abcreg_loclinear(target, parn, ss; tol = tl, transformation = "log")
        bin = quiet(() -> MKtest.read_gz(joinpath(td, "o_log.0.log.post.gz")))
        @test all(isnan, jl[:, 2])
        @test close6(jl[:, [1; 3:6]], bin[:, [1; 3:6]])
    else
        @test_skip "ABCreg binary not available"
    end

    # -- end to end: the same table rows and Poisson draws on disk and in memory ----
    best = BEST
    p = P500
    td = mktempdir()
    tab = TAB300.file
    K = 4
    sfs = [best["sfs"] for _ in 1:K]
    div = [[3000.0, 2000.0 * (1 - 0.06k), best["div"][3], best["div"][4]] for k in 1:K]
    St = length(p().dac)
    kw = (h5_file = tab, summstat_size = 5000, alpha = nothing, B_bins = nothing,
          sampling = MKtest.RjSampling(), row_prior = MKtest.FlatAlphaPrior())
    files(folder; o...) = (seed!(2); quiet(() -> MKtest.summary_statistics(p(), sfs, div;
        h5_file = tab, summstat_size = 5000, output_folder = folder, o...)))
    mem(method; o...) = (seed!(2); quiet(() -> MKtest._abc_in_memory(p(), sfs, div;
        abc_method = method, tol = 200, merge(kw, o)...)))
    Ms(x) = Matrix.(x)

    for (lab, o) in (("default", (;)),
                     ("v1.0.0 sampling", (sampling = MKtest.IndependentSampling(),
                                          row_prior = MKtest.TableRowPrior())))
        fo = joinpath(td, replace(lab, " " => "_"))
        files(fo; o...)
        pf = quiet(() -> MKtest.abc(output_folder = fo, S = St, tol = 200))["loclinear"]
        pm = mem("abc"; o...)
        @test length(pm) == length(pf) == K
        @test Ms(pm) == Ms(pf)
        @test names(pm[1]) == names(pf[1])
        if have_reg
            pr = quiet(() -> MKtest.ABCreg(output_folder = fo, S = St, tol = 200))
            mr = mem("abcreg"; o...)
            @test length(mr) == length(pr) == K
            @test all(close6(a, b) for (a, b) in zip(Ms(mr), Ms(pr)))
            @test names(mr[1]) == names(pr[1])
        end
    end
    have_reg || @test_skip "ABCreg binary not available"

    # the matrix method of abc_loclinear is the file method on the same numbers
    fd = joinpath(td, "default")
    f1 = MKtest._dataset_files(fd, "alphas_")[1]
    s1 = MKtest._dataset_files(fd, "summstat_")[1]
    tbl = quiet(() -> MKtest.read_gz(s1))
    Pc = size(tbl, 2) == 12 + St ? 12 : 8
    @test MKtest.abc_loclinear(f1, s1; P = Pc, tol = 200) ==
          MKtest.abc_loclinear(vec(MKtest.read_gz(f1)), tbl[:, 1:Pc], tbl[:, (Pc+1):end]; tol = 200)

    # GSL's quantile from two order statistics equals it on the sorted data
    for x in (rand(MersenneTwister(4), 1001), [3.0, 1.0, 1.0, 2.0, 2.0], [5.0])
        for f in (0.0, 0.01, 0.5, 0.9999, 1.0)
            @test MKtest._gsl_quantile_unsorted(x, f) == MKtest._gsl_quantile(sort(x), f)
        end
    end

    # summary buffers are reused across datasets: none allocated once the pool has one,
    # and a reused buffer gives the summaries of a fresh one
    st = quiet(() -> (seed!(2); MKtest._summaries_setup(p(), sfs, div; kw...)))
    pool = MKtest._BufferPool()
    asm(buf, i) = Matrix(MKtest._assemble_summaries!(buf, st.sfs_p[i], st.divergence_p[i], st.neut,
        st.sel, st.dsdn, st.gn, st.sh, st.gL, st.gH, st.B, st.summstat_size, st.samplers[i]))
    seed!(9); fresh = [asm(MKtest._dataset_buffer(st), i) for i in 1:K]
    seed!(9); @test [MKtest._with_buffer(b -> asm(b, i), pool, st) for i in 1:K] == fresh
    @test length(pool.free) == 1
    reuse() = MKtest._with_buffer(b -> nothing, pool, st)
    reuse()
    @test @allocated(reuse()) < 10_000

    # the in-memory ABCreg port asks to cite ABCreg, as the binary route does
    seed!(2)
    @test_logs (:info, r"10\.1186/1471-2156-10-35") match_mode = :any MKtest._abc_in_memory(
        p(), sfs, div; abc_method = "abcreg", tol = 200, kw...)

    # through abcmk: in memory gives the file route's summary
    run_abcmk(method, im) = (seed!(2); quiet(() -> MKtest.abcmk(p(), sfs, div; h5_file = tab,
        summstat_size = 5000, tol = 200, abc_method = method, in_memory = im,
        output_folder = im ? "" : joinpath(td, "abcmk_$method"))))
    @test run_abcmk("abc", true)[:inference] == run_abcmk("abc", false)[:inference]
    # a single set as parse_sfs returns it (Matrix, Vector) equals the wrapped call
    one(s, d) = (seed!(3); quiet(() -> MKtest.abcmk(p(), s, d; h5_file = tab, summstat_size = 2000,
                                                    tol = 100, in_memory = true)))
    @test one(sfs[1], div[1])[:inference] == one([sfs[1]], [div[1]])[:inference]
    if have_reg
        a, b = run_abcmk("abcreg", true)[:inference], run_abcmk("abcreg", false)[:inference]
        @test isapprox(Matrix(a), Matrix(b); rtol = 1e-4)
    end
    # the file route still needs a folder
    @test_throws ArgumentError MKtest.abcmk(p(), sfs, div; h5_file = tab, summstat_size = 5000,
                                            in_memory = false)
end

# A table loaded once (`load_rates`) and passed as `rates = t`. Oracles: the `h5_file` route
# on the same table and seed (identical files and posteriors), and the Dict-counting
# `_row_weights` that the array version replaced (copied below, identical weights).
@testset "load_rates" begin
    function weights_dict(p::MKtest.FlatAlphaPrior, dsdn, idx; B = nothing)
        key = Vector{Tuple{Float64,Int}}(undef, length(idx))
        cnt = Dict{Tuple{Float64,Int},Int}(); nB = Dict{Float64,Int}()
        for (k, i) in enumerate(idx)
            α = (dsdn[i, 3] + dsdn[i, 4]) / dsdn[i, 2]
            g = p.by_B ? Float64(B[i]) : 0.0
            key[k] = (g, isfinite(α) ? floor(Int, α / p.bin) : -1)
            nB[g] = get(nB, g, 0) + 1
            key[k][2] >= 0 && (cnt[key[k]] = get(cnt, key[k], 0) + 1)
        end
        nbin = Dict{Float64,Int}()
        for k in keys(cnt); nbin[k[1]] = get(nbin, k[1], 0) + 1; end
        share(g) = p.by_B ? nB[g] / (length(idx) * nbin[g]) : 1.0
        return [x[2] >= 0 ? share(x[1]) / max(cnt[x], p.floor * nB[x[1]]) : 0.0 for x in key]
    end
    rng = MersenneTwister(5)
    R = 20_000
    dsdn = hcat(rand(rng, R) .+ 0.5, rand(rng, R) .+ 1.0, rand(rng, R) .* 0.4, rand(rng, R) .* 0.4)
    dsdn[1:50, 2] .= 0.0                               # α = Inf -> weight 0
    dsdn[51:80, 3] .= -2.0                             # α < 0 -> weight 0
    Bc = rand(rng, [200.0, 500.0, 999.0, -0.0, 0.0], R)
    idx = sort(unique(rand(rng, 1:R, 15_000)))
    for pr in (MKtest.FlatAlphaPrior(), MKtest.FlatAlphaPrior(by_B = false),
               MKtest.FlatAlphaPrior(floor = 0.05), MKtest.FlatAlphaPrior(bin = 0.05, floor = 0.0))
        @test MKtest._row_weights(pr, dsdn, idx; B = Bc) == weights_dict(pr, dsdn, idx; B = Bc)
    end
    @test all(iszero, MKtest._row_weights(MKtest.FlatAlphaPrior(), dsdn, 1:50; B = Bc))

    best = BEST
    p = P500
    td = mktempdir()
    tab = TAB300.file
    sfs = [best["sfs"] for _ in 1:3]
    div = [[3000.0, 2000.0 * (1 - 0.06k), best["div"][3], best["div"][4]] for k in 1:3]

    t = MKtest.load_rates(tab, p())
    @test t isa MKtest.RatesTable
    @test occursin("RatesTable(", sprint(show, t)) && occursin("rows", sprint(show, t))
    @test t.dac == MKtest.load_rates(tab, p()).dac
    @test isempty(t.weights)

    # summary files: h5_file, rates = t (twice: the second takes the cached weights) and
    # rates = path are the same files
    function sumfiles(src)
        f = mktempdir()
        seed!(2)
        quiet(() -> MKtest.summary_statistics(p(), sfs, div; src..., summstat_size = 4000, output_folder = f))
        [MKtest.read_gz(joinpath(f, "$(pre)$k.txt.gz")) for k in 1:3 for pre in ("summstat_", "alphas_")]
    end
    ref = sumfiles((h5_file = tab,))
    @test sumfiles((rates = t,)) == ref
    @test length(t.weights) == 1
    w1 = only(values(t.weights))
    @test sumfiles((rates = t,)) == ref
    @test only(values(t.weights)) === w1              # reused, not recomputed
    @test sumfiles((rates = tab,)) == ref

    # a filtered prior gets its own weights
    seed!(2)
    quiet(() -> MKtest.summary_statistics(p(), sfs, div; rates = t, summstat_size = 2000,
                                          output_folder = mktempdir(), B_bins = [0.3, 0.999]))
    @test length(t.weights) == 2

    # abcmk, in memory and on files: the same posteriors as h5_file
    am(src; o...) = (seed!(3); quiet(() -> MKtest.abcmk(p(), sfs, div; src..., summstat_size = 4000,
                                                       tol = 150, abc_method = "abc", o...)))
    @test am((rates = t,); in_memory = true)[:inference] == am((h5_file = tab,); in_memory = true)[:inference]
    @test am((rates = t,); output_folder = mktempdir(), in_memory = false)[:inference] ==
          am((h5_file = tab,); output_folder = mktempdir(), in_memory = false)[:inference]

    # exactly one source, and a loaded table must match param
    @test_throws ArgumentError MKtest.abcmk(p(), sfs, div; summstat_size = 4000, in_memory = true)
    @test_throws ArgumentError MKtest.abcmk(p(), sfs, div; h5_file = tab, rates = t,
                                            summstat_size = 4000, in_memory = true)
    @test_throws ArgumentError quiet(() -> MKtest.summary_statistics(p(), sfs, div; summstat_size = 4000,
                                                                     output_folder = mktempdir()))
    p_n = MKtest.parameters(N = 500, n = 40, dac = [2, 4, 5, 10, 20, 40])
    @test_throws ArgumentError quiet(() -> MKtest.abcmk(p_n, sfs, div; rates = t, summstat_size = 4000,
                                                        in_memory = true))
    p_dac = MKtest.parameters(N = 500, n = 50, dac = [2, 4, 5, 10, 20, 30])
    @test_throws ArgumentError quiet(() -> MKtest.abcmk(p_dac, sfs, div; rates = t, summstat_size = 4000,
                                                        in_memory = true))
end

# The ω block is a property of the data: present when Ln > 0 and Ls > 0, absent (NaN ω,
# ABC on 8 parameters) without site counts. A draw with Ds = 0 leaves ω undefined for that
# row only, and the row is dropped like any non-finite row, as alpha's are (dn = 0).
# Three datasets on the same draws: one with ample divergence, one with Dn + Ds = 10 (many
# Ds = 0 draws, keeps its ω block with fewer rows) and one without site counts. Oracles:
# the ABC on the last one's 8-parameter table (the regression of every other column is
# independent of ω), the file route with the binary (6 digits) and with `abc` (exact), on
# the same draws.
@testset "the ω block follows the site counts, not the draws" begin
    closeNaN(a, b) = size(a) == size(b) &&
        all((isnan(x) && isnan(y)) || abs(x - y) <= 1e-5 * abs(y) + 1e-9 for (x, y) in zip(a, b))
    reg = joinpath(pkgdir(MKtest), "scripts", "reg")
    have_reg = Sys.islinux() && isfile(reg)
    best = BEST
    p = P500
    td = mktempdir()
    tab = TAB300.file
    sfs = [best["sfs"], best["sfs"], best["sfs"]]
    div = [best["div"],
           [6.0, 4.0, best["div"][3], best["div"][4]],        # Dn + Ds = 10: Ds = 0 draws
           [best["div"][1], best["div"][2], 0.0, 0.0]]         # no site counts: no ω
    S = length(p().dac)
    kw = (h5_file = tab, summstat_size = 5000, alpha = nothing, B_bins = nothing,
          sampling = MKtest.RjSampling(), row_prior = MKtest.FlatAlphaPrior())

    # summaries, row by row: the Ds = 0 rows leave, the ω columns stay
    seed!(2)
    st = quiet(() -> MKtest._summaries_setup(p(), sfs, div; kw...))
    buf = MKtest._dataset_buffer(st)
    seed!(3)
    out2 = MKtest._assemble_summaries!(buf, st.sfs_p[2], st.divergence_p[2], st.neut, st.sel, st.dsdn,
                                       st.gn, st.sh, st.gL, st.gH, st.B, st.summstat_size, st.samplers[2])
    finite_row(i) = all(isfinite, view(buf.α_mat, i, :)) && all(isfinite, view(buf.ω_mat, i, :)) &&
                    all(isfinite, view(buf.α_summaries, i, :))
    nds0 = count(==(0), buf.sampled_ds)
    @test nds0 > 0 && any(i -> any(isnan, view(buf.ω_mat, i, :)), 1:st.N)   # the case under test
    @test size(out2, 2) == 12 + S
    @test size(out2, 1) == count(finite_row, 1:st.N) == count(buf.mask)
    @test all(i -> buf.mask[i] == finite_row(i), 1:st.N)
    @test !any(buf.mask[i] for i in 1:st.N if buf.sampled_ds[i] == 0)
    @test all(isfinite, out2)
    # the α columns of the kept rows are the draws themselves, untouched
    @test out2[:, 1:3] == round.(buf.α_mat[buf.mask, :], digits = 5)
    # without site counts: 8 columns, every finite-α row kept
    buf3 = MKtest._dataset_buffer(st)
    seed!(3)
    out3 = MKtest._assemble_summaries!(buf3, st.sfs_p[3], st.divergence_p[3], st.neut, st.sel, st.dsdn,
                                       st.gn, st.sh, st.gL, st.gH, st.B, st.summstat_size, st.samplers[3])
    @test size(out3, 2) == 8 + S
    @test size(out3, 1) == count(i -> all(isfinite, view(buf3.α_mat, i, :)) &&
                                      all(isfinite, view(buf3.α_summaries, i, :)), 1:st.N)

    mem(m) = (seed!(2); quiet(() -> MKtest._abc_in_memory(p(), sfs, div; abc_method = m, tol = 200, kw...)))
    fo = joinpath(td, "files")
    seed!(2)
    quiet(() -> MKtest.summary_statistics(p(), sfs, div; h5_file = tab, summstat_size = 5000,
                                          output_folder = fo))
    files = MKtest._dataset_files(fo, "summstat_")
    @test MKtest._count_cols.(files) == [12 + S, 12 + S, 8 + S]
    n_rows = [quiet(() -> MKtest._with_gz(countlines, f)) for f in files]
    @test n_rows[2] < 5000 && !any(isnan, quiet(() -> MKtest.read_gz(files[2])))
    t3 = quiet(() -> MKtest.read_gz(files[3]))
    a3 = vec(quiet(() -> MKtest.read_gz(MKtest._dataset_files(fo, "alphas_")[3])))
    ω = 4:7
    rest = [1:3; 8:12]

    names12 = MKtest._abc_param_names(12)
    for (m, direct) in (("abcreg", MKtest.abcreg_loclinear(a3, t3[:, 1:8], t3[:, 9:end]; tol = 200)),
                        ("abc", MKtest.abc_loclinear(a3, t3[:, 1:8], t3[:, 9:end]; tol = 200)[2]))
        post = mem(m)
        @test length(post) == 3
        @test names(post[1]) == names(post[2]) == names(post[3]) == string.(names12)
        @test !any(isnan, Matrix(post[1])) && !any(isnan, Matrix(post[2]))
        @test all(isnan, Matrix(post[3][:, ω])) && !any(isnan, Matrix(post[3][:, rest]))
        @test Matrix(post[3][:, rest]) == direct
    end

    # file route: the same posteriors, NaN ω included
    fa = quiet(() -> MKtest.abc(output_folder = fo, S = S, tol = 200))["loclinear"]
    @test all(isequal(Matrix(a), Matrix(b)) for (a, b) in zip(fa, mem("abc")))
    @test names(fa[3]) == string.(names12)
    if have_reg
        fr = quiet(() -> MKtest.ABCreg(output_folder = fo, S = S, tol = 200))
        @test all(closeNaN(Matrix(a), Matrix(b)) for (a, b) in zip(fr, mem("abcreg")))
    else
        @test_skip "ABCreg binary not available"
    end

    # summaries: NaN statistics for the unestimated ω, the usual ones elsewhere
    post = mem("abcreg")
    res = quiet(() -> MKtest.summary_abc(post; stat = MKtest.PosteriorMedian()))
    @test all(isnan, Matrix(res[:inference][3:3, ω])) && !any(isnan, Matrix(res[:inference][:, rest]))
    @test !any(isnan, Matrix(res[:inference][1:2, :]))
    @test isequal(Matrix(res[:inference][1:1, :]), Matrix(quiet(() -> MKtest.summary_abc(post[1:1]; stat = MKtest.PosteriorMedian()))[:inference]))
    @test occursin("NaN", res[:quantiles][3, :ωₐ])
    # and through abcmk, with summary_abc's statistic
    r = quiet(() -> (seed!(2); MKtest.abcmk(p(), sfs, div; h5_file = tab, summstat_size = 5000,
                                            tol = 200, in_memory = true,
                                            stat = MKtest.PosteriorMedian(), components = MKtest.Marginal())))
    @test isequal(Matrix(r[:inference]), Matrix(res[:inference]))
end

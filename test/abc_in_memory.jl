# In-memory ABC (`abcmk(...; in_memory = true)`) and the Julia port of ABCreg.
#
# Oracles: the ABCreg binary (scripts/reg) and the file `abc` route, on the very same
# summaries. Thresholds fixed before running: in-memory "abc" equals the file route
# exactly; `abcreg_loclinear` equals the binary to its printed precision (6 significant
# digits, so |a - b| <= 1e-5 |b| + 1e-9) with the same accepted rows.
#
# Standalone: julia --project=. -e 'using MKtest, Test, JLD2, Logging; using Random: seed!,
#   MersenneTwister; using DelimitedFiles: writedlm; include("test/abc_in_memory.jl")'

@testset "ABCreg in memory" begin
    reg = joinpath(pkgdir(MKtest), "scripts", "reg")
    have_reg = Sys.islinux() && isfile(reg)
    close6(a, b) = size(a) == size(b) && all(abs.(a .- b) .<= 1e-5 .* abs.(b) .+ 1e-9)
    quiet(f) = Logging.with_logger(f, Logging.NullLogger())

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
    best = JLD2.load(joinpath(@__DIR__, "fixtures", "abcmk_best_model.jld2"))
    p() = MKtest.parameters(N = 500, n = 50, dac = [2, 4, 5, 10, 20, 50])
    td = mktempdir()
    tab = joinpath(td, "t.jld2")
    seed!(1)
    quiet(() -> MKtest.rates(p(); iterations = 300, output_file = tab))
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
    @test_throws ArgumentError MKtest.abcmk(p(), sfs, div; h5_file = tab, summstat_size = 5000)
end

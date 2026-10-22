# parse_sfs, the classic MK tests, the SFS binning, matched control sets and the small
# methods. Group "mk" of runtests.jl.

@testset "MK tests and data" begin

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
end

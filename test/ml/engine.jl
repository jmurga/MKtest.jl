# dfe_ml engine: the testsets ported from dev (9d752e40), rewritten for the new
# interface with their meaning kept, plus the output and axis checks of the move.
#
# Shared synthetic data (n = 40): Pn_j = 2000/j, Ps_j = 3000/j.

using MKtest
using Test
using DataFrames: DataFrame, nrow, names
using JLD2
using SpecialFunctions: loggamma
import Distributions
using Distributions: logpdf

import MKtest: GammaZero, GammaExpo, GammaGamma, ExpoZero, ExpoExpo, ExpoGamma, DisplGamma,
               ScaledBeta, LogNormalZero, NeuGamma, FGMBesselK, FixedSpikes, WithLethal,
               ExactSojourn, GaussLegendreSojourn, RombergIntegrator, LogGammaGrid, GrapesGrid,
               PoissonKernel, NegativeBinomialKernel, MidpointDensity, IntervalCDF, DFECache

const ML_JS = collect(1.0:39.0)
const ML_SFS = hcat(ML_JS, round.(2000.0 ./ ML_JS), round.(3000.0 ./ ML_JS))
const ML_DIV = [900.0, 1200.0, 1.0e6, 3.0e5]

# A DFE quadrature defined outside the package (the axis is open): the density
# averaged at the two cell bounds.
struct TrapezoidDensity <: MKtest.AbstractDFEQuadrature end
function MKtest.fill_dfe_vals!(out, ::TrapezoidDensity, m::MKtest.AbstractDFEModel, p, pre)
    lo, hi = similar(out), similar(out)
    MKtest.dfe_density!(lo, m, p, pre.s_lo)
    MKtest.dfe_density!(hi, m, p, pre.s_hi)
    @inbounds for i in eachindex(out)
        out[i] = (lo[i] + hi[i]) / 2
    end
    return out
end

@testset "dfe_ml engine" begin

    # --- integrators ----------------------------------------------------------
    # RombergIntegrator transliterates GRAPES and keeps its defect (the truncated
    # interval drops ~n*QUASI_ZERO of bin 1); both sides are pinned.
    @testset "integrators" begin
        logbinom(n, j) = sum(log, 1:n) - sum(log, 1:j) - sum(log, 1:(n-j))
        ex = ExactSojourn(); gl = GaussLegendreSojourn(); rb = RombergIntegrator()

        # exact identity counts(S -> 0, j, n) = 1/j
        for (n, j) in ((20, 1), (200, 7), (1322, 3))
            lc = logbinom(n, j)
            @test MKtest.expected_count(ex, 1e-12, j, n, lc) ≈ 1/j rtol=1e-12
            @test MKtest.expected_count(gl, 1e-12, j, n, lc) ≈ 1/j rtol=1e-11
        end

        # two independent methods agree to 1e-10 (measured 1.4e-11 on dev)
        for n in (20, 100), j in (1, 2, n÷2, n-1), S in (-2000.0, -500.0, -5.0, 5.0, 500.0, 10000.0)
            lc = logbinom(n, j)
            e = MKtest.expected_count(ex, S, j, n, lc)
            e > 1e-6 * (n / (j * (n - j))) || continue
            @test abs(MKtest.expected_count(gl, S, j, n, lc) - e) / e < 1e-10
        end

        @test isnan(MKtest._sojourn_density(0.01, -2000.0))
        @test isfinite(MKtest.expected_count(ex, -2000.0, 3, 100, logbinom(100, 3)))
        @test isfinite(MKtest.expected_count(gl, -2000.0, 3, 100, logbinom(100, 3)))

        # the Romberg hook is the historical expression, ulp for ulp
        for (n, j, S) in ((40, 1, 5.0), (40, 20, -50.0), (100, 3, 500.0))
            lc = logbinom(n, j)
            f = x -> MKtest._sojourn_density(x, S) * exp(lc + j * log(x) + (n - j) * log(1.0 - x))
            @test MKtest.expected_count(rb, S, j, n, lc) ===
                  MKtest.integrate(rb, f, MKtest.QUASI_ZERO, 1.0 - MKtest.QUASI_ZERO) / 2.0
        end
        # ...and still wrong by the documented amount
        @test abs(MKtest.expected_count(rb, 5.0, 1, 661, logbinom(661, 1)) /
                  MKtest.expected_count(ex, 5.0, 1, 661, logbinom(661, 1)) - 1) > 0.05

        # the cache is keyed by value on the integrator
        c = DFECache()
        a = MKtest.get_precomputed(c, 20; integrator = rb)
        b = MKtest.get_precomputed(c, 20; integrator = ex)
        @test length(c) == 2
        @test a.expected_counts != b.expected_counts
        @test MKtest.get_precomputed(c, 20; integrator = RombergIntegrator()) === a
        @test MKtest.dn_correction(a, 100.0) != MKtest.dn_correction(b, 100.0)

        # the S < -1500 refusal belongs to Romberg
        wide = LogGammaGrid(; int_bounds = (1e-4, 1000.0), Npts = 40)
        @test_throws ErrorException MKtest.get_precomputed(c, 20; integrator = rb, grid = wide)
        @test all(isfinite, MKtest.get_precomputed(c, 20; integrator = ex, grid = wide).expected_counts)
    end

    # --- the GRAPES integrator must not move ----------------------------------
    # test/fixtures/integrator_baseline.jld2 is dev's data/reg/integrator_baseline.jld2
    # (commit 2daa9d1e, 2026-09-08, scripts/gen_integrator_baseline.jl; sha256
    # bad7ca5f...d3bd203). A missing file is a failure.
    @testset "regression: RombergIntegrator is frozen" begin
        base = joinpath(@__DIR__, "..", "fixtures", "integrator_baseline.jld2")
        @test isfile(base)
        if isfile(base)
            d = JLD2.load(base)
            rb = RombergIntegrator()
            # FNV-1a over the raw bits (Base.hash on arrays is not stable across versions)
            function bits_hash(m)
                h = 0xcbf29ce484222325
                @inbounds for x in m
                    u = reinterpret(UInt64, x)
                    for k in 0:7
                        h = (h ⊻ ((u >> (8k)) & 0xff)) * 0x100000001b3
                    end
                end
                return h
            end
            c = DFECache()
            for n in (20, 40, 100)
                m = MKtest.get_precomputed(c, n; integrator = rb).expected_counts
                @test bits_hash(m) == d["expected_counts_n$(n)_hash"]
                ri = d["expected_counts_n$(n)_probe_rows"]
                ci = d["expected_counts_n$(n)_probe_cols"]
                @test m[ri, ci] == d["expected_counts_n$(n)_probe"]
            end
            for n in (20, 200)
                t = MKtest.DnCorrectionTable(n, rb)
                @test [MKtest.dn_correction(t, S) for S in d["dn_correction_S"]] == d["dn_correction_n$(n)"]
            end
        end
    end

    # --- likelihood kernel ----------------------------------------------------
    @testset "likelihood kernel" begin
        pk = PoissonKernel(); nb = NegativeBinomialKernel()
        for (k, mu, lf) in ((0.0, 3.0, 0.0), (7.0, 6.5, loggamma(8.0)), (16745.0, 16000.0, loggamma(16746.0)))
            @test MKtest.log_kernel(pk, k, mu, lf, NaN) === k * log(mu) - mu - lf
        end
        @test MKtest.log_kernel(pk, 7.0, 6.5, 0.0, 1.0) === MKtest.log_kernel(pk, 7.0, 6.5, 0.0, 1e9)

        # against Distributions.jl (measured worst 3.7e-12 on dev)
        for k in (0, 1, 7, 1000, 16745), mu in (0.1, 3.0, 950.0, 16000.0), a in (1e-4, 1.0, 100.0, 1000.0)
            ours = MKtest.log_kernel(nb, Float64(k), mu, loggamma(k + 1), a)
            ref  = logpdf(Distributions.NegativeBinomial(a, a / (a + mu)), k)
            @test abs(ours - ref) <= 1e-10 * max(abs(ref), 1e-12)
        end
        # Poisson limit at rate O(1/a)
        gap(a) = maximum(abs(MKtest.log_kernel(nb, Float64(k), mu, loggamma(k + 1), a) -
                             MKtest.log_kernel(pk, Float64(k), mu, loggamma(k + 1), NaN)) /
                         abs(MKtest.log_kernel(pk, Float64(k), mu, loggamma(k + 1), NaN))
                         for k in (0, 1, 7, 1000, 16745), mu in (0.1, 3.0, 950.0, 16000.0))
        g6, g7, g8 = gap(1e6), gap(1e7), gap(1e8)
        @test g8 < g7 < g6 < 0.05
        @test 5.0 < g6 / g7 < 20.0
        @test 5.0 < g7 / g8 < 20.0

        # the kernel's parameter is appended after the DFE and the orientation error
        lp = MKtest.build_layout(GammaZero(), pk, true)
        ln = MKtest.build_layout(GammaZero(), nb, true)
        @test lp.names == [:negGshape, :negGmean, :syno_orientation_error]
        @test ln.names == [:negGshape, :negGmean, :syno_orientation_error, :a]
        @test length(ln.opt_idx) == length(lp.opt_idx) + 1
        @test lp.kernel_idx == 0 && ln.kernel_idx == 4
        @test ln.n_spec == lp.n_spec == 2
        @test ln.log_scale == [false, true, false, true]
        # `fixed` reaches the kernel parameter like any other
        lf = MKtest.build_layout(GammaZero(), nb, true; fixed = (a = 500.0,))
        @test !(4 in lf.opt_idx)
        # unknown names are errors, not silent no-ops
        @test_throws ArgumentError MKtest.build_layout(GammaZero(), pk, true; fixed = (foo = 1.0,))
        @test_throws ArgumentError MKtest.build_layout(GammaZero(), pk, true; bounds = (negGmean = (10.0, 1.0),))
        @test_throws ArgumentError MKtest.build_layout(GammaZero(), pk, false; fixed = (syno_orientation_error = 0.0,))
    end

    # --- discrete DFE components and the extra shapes ---------------------------
    @testset "discrete components and DFE shapes" begin
        for m in (GammaZero(), GammaExpo(), GammaGamma(), ExpoZero(), ExpoExpo(), ExpoGamma(),
                  DisplGamma(), ScaledBeta(), LogNormalZero())
            @test isempty(collect(MKtest.discrete_components(m, first(MKtest.dfe_start_points(m)))))
        end

        # a point mass at S = 0 is the neutral spectrum, 1/j exactly
        pre = MKtest.get_precomputed(DFECache(), 40)
        row = MKtest.discrete_row(pre, 0.0)
        @test length(row) == 39
        @test maximum(abs(row[j] - 1 / j) / (1 / j) for j in 1:39) < 1e-12
        # a spike far outside the grid is evaluated at its true value
        deep = MKtest.discrete_row(pre, -10000.0)
        @test all(isfinite, deep) && all(>(0.0), deep)
        @test deep[1] / deep[20] > 100

        @test MKtest.dfe_param_names(NeuGamma()) == [:p_neu, :negGshape, :negGmean]
        @test collect(MKtest.discrete_components(NeuGamma(), [0.3, 0.2, 1000.0])) == [(0.0, 0.3)]
        sgrid = collect(range(-500.0, 100.0; length = 401))
        a, b = zeros(401), zeros(401)
        MKtest.dfe_density!(a, GammaZero(), [0.2, 1000.0], sgrid)
        MKtest.dfe_density!(b, NeuGamma(), [0.0, 0.2, 1000.0], sgrid)
        @test a == b
        @test MKtest.dfe_mass_below(NeuGamma(), [0.0, 0.2, 1000.0], 5.0) ==
              MKtest.dfe_mass_below(GammaZero(), [0.2, 1000.0], 5.0)
        for bb in (1e-9, 0.5, 5.0, 500.0)
            @test MKtest.dfe_mass_below(NeuGamma(), [0.3, 0.2, 1000.0], bb) ≈
                  0.3 + 0.7 * MKtest.dfe_mass_below(GammaZero(), [0.2, 1000.0], bb)
            @test MKtest.dfe_mass_above(NeuGamma(), [0.3, 0.2, 1000.0], bb) ≈
                  0.7 * MKtest.dfe_mass_above(GammaZero(), [0.2, 1000.0], bb)
        end

        m6 = FixedSpikes()
        @test length(m6.positions) == 6
        @test MKtest.dfe_param_names(m6) == [:w1, :w2, :w3, :w4, :w5]
        for v in ([0.5, 0.5, 0.5, 0.5, 0.5], [0.01, 0.99, 0.5, 0.2, 0.8],
                  [1e-6, 1e-6, 1e-6, 1e-6, 1e-6], [1 - 1e-6, 0.3, 0.3, 0.3, 0.3])
            w = MKtest.spike_weights(m6, v)
            @test length(w) == 6 && all(>=(0.0), w)
            @test sum(w) ≈ 1.0 atol=1e-12
        end
        @test all(≈(1 / 6; atol = 1e-12), MKtest.spike_weights(m6, first(MKtest.dfe_start_points(m6))))
        z = ones(401)
        MKtest.dfe_density!(z, m6, [0.5, 0.5, 0.5, 0.5, 0.5], sgrid)
        @test all(iszero, z)
        @test sum(w for (_, w) in MKtest.discrete_components(m6, [0.5, 0.5, 0.5, 0.5, 0.5])) ≈ 1.0
        m3 = FixedSpikes([0.0, -50.0, -5000.0])
        @test MKtest.dfe_param_names(m3) == [:w1, :w2]
        @test MKtest.dfe_name(m3) == "FixedSpikes"

        # dev raised a MethodError here: FixedSpikes' density was typed ::Vector, and
        # WithLethal hands it a view.
        wf = WithLethal(FixedSpikes())
        z2 = ones(401)
        MKtest.dfe_density!(z2, wf, [0.5, 0.5, 0.5, 0.5, 0.5, 0.1], sgrid)
        @test all(iszero, z2)
        @test sum(w for (_, w) in MKtest.discrete_components(wf, [0.5, 0.5, 0.5, 0.5, 0.5, 0.1])) ≈ 0.9

        # FGMBesselK: GRAPES' bounds, a normalised density, GRAPES' S = 0 convention
        @test MKtest.dfe_param_names(FGMBesselK()) == [:BKm, :BKnsz2, :BKsigma, :BKscale]
        @test MKtest.dfe_parameters(FGMBesselK()).BKm == (1.0, 10.0)
        @test MKtest.dfe_log_scale(FGMBesselK()) == (false, false, false, true)
        fine = collect(range(-1000.0, 1000.0; length = 40001))
        dd = zeros(length(fine))
        MKtest.dfe_density!(dd, FGMBesselK(), [2.0, 100.0, 0.3, 100.0], fine)
        @test all(isfinite, dd) && all(>=(0.0), dd)
        @test sum(dd) * step(range(-1000.0, 1000.0; length = 40001)) ≈ 1.0 atol=2e-3
        z0 = zeros(1); MKtest.dfe_density!(z0, FGMBesselK(), [2.0, 100.0, 0.3, 100.0], [0.0])
        zm = zeros(2); MKtest.dfe_density!(zm, FGMBesselK(), [2.0, 100.0, 0.3, 100.0],
                                           [-MKtest.QUASI_ZERO, MKtest.QUASI_ZERO])
        @test z0[1] ≈ (zm[1] + zm[2]) / 2
    end

    # --- WithLethal and the two alpha definitions -------------------------------
    @testset "WithLethal and the alpha split" begin
        wl = WithLethal(GammaZero())
        @test MKtest.dfe_param_names(wl) == [:negGshape, :negGmean, :p_lethal]
        @test MKtest.dfe_parameters(wl).p_lethal == (0.0, 0.99)
        @test MKtest.dfe_log_scale(wl) == (false, true, false)
        @test MKtest.dfe_name(wl) == "WithLethal{GammaZero}"
        @test all(length(sp) == 3 for sp in MKtest.dfe_start_points(wl))

        sgrid = collect(range(-500.0, 100.0; length = 401))
        base, wrapped = zeros(401), zeros(401)
        MKtest.dfe_density!(base, GammaZero(), [0.2, 1000.0], sgrid)
        MKtest.dfe_density!(wrapped, wl, [0.2, 1000.0, 0.0], sgrid)
        @test base == wrapped
        @test MKtest.dfe_mass_above(wl, [0.2, 1000.0, 0.0], 500.0) ==
              MKtest.dfe_mass_above(GammaZero(), [0.2, 1000.0], 500.0)
        MKtest.dfe_density!(wrapped, wl, [0.2, 1000.0, 0.35], sgrid)
        @test wrapped ≈ 0.65 .* base
        @test MKtest.dfe_mass_above(wl, [0.2, 1000.0, 0.35], 500.0) ≈
              0.65 * MKtest.dfe_mass_above(GammaZero(), [0.2, 1000.0], 500.0)
        @test MKtest.dfe_mass_below(wl, [0.2, 1000.0, 0.35], 5.0) ≈
              0.65 * MKtest.dfe_mass_below(GammaZero(), [0.2, 1000.0], 5.0)

        wn = WithLethal(NeuGamma())
        @test MKtest.dfe_param_names(wn) == [:p_neu, :negGshape, :negGmean, :p_lethal]
        @test collect(MKtest.discrete_components(wn, [0.4, 0.2, 1000.0, 0.25])) == [(0.0, 0.75 * 0.4)]
        # wrapping twice would give two parameters named p_lethal: refused
        ww = WithLethal(WithLethal(GammaZero()))
        @test_throws ArgumentError MKtest.dfe_param_names(ww)

        # the alpha split: with a free Dn, alpha_div IS alpha
        r = MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero())
        @test r.alpha_div[1] ≈ r.alpha[1]
        # no beneficial class, so the DFE predicts no adaptive substitution
        @test r.alpha_dfe[1] == 0.0
        r2 = MKtest.dfe_ml(ML_SFS, ML_DIV, GammaExpo())
        @test r2.alpha_dfe[1] > 0.0

        cmp = MKtest.dfe_model_comparison(vcat(r, r2; cols = :union))
        @test sort(cmp.aic.model) == ["GammaExpo", "GammaZero"]
        @test nrow(cmp.lrt) <= 1
        if nrow(cmp.lrt) == 1
            @test cmp.lrt.model_null[1] == "GammaZero" && cmp.lrt.model_alt[1] == "GammaExpo"
            @test cmp.lrt.df[1] == 2
        end

        # at p_lethal = 0 (fixed) the wrapper is its base, fit included
        rl = MKtest.dfe_ml(ML_SFS, ML_DIV, wl; fixed = (p_lethal = 0.0,))
        @test rl.loglik[1] ≈ r.loglik[1] atol=1e-9
        @test rl.n_params[1] == r.n_params[1]

        # no divergence data: nothing throws, alpha_div is undefined
        r3 = MKtest.dfe_ml(ML_SFS, [0.0, 0.0, 1.0e6, 3.0e5], GammaZero(); div_model = MKtest.NoDivergence())
        @test isnan(r3.alpha_div[1])
        @test isfinite(r3.loglik[1]) && r3.loglik[1] < 0
        # without divergence the SFS alone is scored: the saturated Dn/Ds terms are absent
        @test r3.loglik[1] != MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero()).loglik[1]
    end

    # --- IntervalCDF quadrature -------------------------------------------------
    @testset "IntervalCDF quadrature" begin
        pre = MKtest.get_precomputed(DFECache(), 40; grid = GrapesGrid())
        @test length(pre.s_lo) == length(pre.s_hi) == length(pre.s_grid)
        @test all(pre.s_lo .<= pre.s_grid .<= pre.s_hi)
        # GrapesGrid cells tile exactly
        @test pre.s_lo[2:end] ≈ pre.s_hi[1:(end - 1)]
        @test pre.s_lo[1] ≈ -MKtest.LARGE_S_VALUE
        @test pre.s_hi[end] ≈ MKtest.VERY_LARGE_S_VALUE

        lg = LogGammaGrid(; int_bounds = (1e-4, 1000.0), Npts = 200)
        gpts = MKtest.s_grid_points(lg)
        gw = MKtest.s_grid_widths(lg, gpts)
        llo, lhi = MKtest.s_cell_bounds(lg, gpts, gw)
        @test all(llo .<= gpts .<= lhi)
        @test lhi .- llo ≈ gw
        ineg = findlast(<(0), gpts)
        @test lhi[ineg] < 0 < llo[ineg + 1]

        # every CDF is a CDF, reaching its own total mass
        for (m, total) in ((GammaZero(), 1.0), (GammaExpo(), 1.0), (GammaGamma(), 1.0), (ExpoZero(), 1.0),
                           (LogNormalZero(), 1.0), (NeuGamma(), 1.0), (DisplGamma(), 1.0))
            p = first(MKtest.dfe_start_points(m))
            disc = sum(w for (_, w) in MKtest.discrete_components(m, p); init = 0.0)
            @test MKtest.dfe_cdf_continuous(m, p, 1.0e8) + disc ≈ total atol=1e-6
            @test MKtest.dfe_cdf_continuous(m, p, -1.0e8) ≈ 0.0 atol=1e-8
            @test issorted([MKtest.dfe_cdf_continuous(m, p, S)
                            for S in (-5000.0, -500.0, -50.0, -5.0, -0.5, 0.5, 5.0, 500.0)])
        end
        pb = first(MKtest.dfe_start_points(ScaledBeta()))
        @test MKtest.dfe_cdf_continuous(ScaledBeta(), pb, 30.0) ≈ pb[3]
        @test MKtest.dfe_cdf_continuous(ScaledBeta(), pb, -30.0) == 0.0
        @test issorted([MKtest.dfe_cdf_continuous(ScaledBeta(), pb, S) for S in (-25.0, -24.0, -10.0, 0.0, 10.0, 25.0)])

        # the CDF agrees with a trapezoid of the model's own density
        function trap(m, p, a, b, npts = 100001)
            xs = collect(range(a, b; length = npts)); o = zeros(npts)
            MKtest.dfe_density!(o, m, p, xs)
            (sum(o) - (o[1] + o[end]) / 2) * step(range(a, b; length = npts))
        end
        for m in (GammaZero(), GammaExpo(), LogNormalZero(), DisplGamma())
            p = first(MKtest.dfe_start_points(m))
            for S in (-2000.0, -200.0, -20.0)
                @test MKtest.dfe_cdf_continuous(m, p, S) ≈ trap(m, p, -2.0e5, S) atol=1e-5
            end
        end

        # the midpoint rule loses mass where the density is singular; the CDF rule
        # does not (dev: 21.8% lost at shape 0.1, 5.7% at 0.186)
        function total_weight(quad, p)
            v = similar(pre.s_grid)
            MKtest.fill_dfe_vals!(v, quad, GammaZero(), p, pre)
            sum(v .* pre.s_widths)
        end
        analytic(p) = MKtest.dfe_cdf_continuous(GammaZero(), p, pre.s_hi[end]) -
                      MKtest.dfe_cdf_continuous(GammaZero(), p, pre.s_lo[1])
        for sh in (0.1, 0.186, 0.3)
            p = [sh, 2000.0]
            @test total_weight(IntervalCDF(), p) ≈ analytic(p) rtol=1e-10
            @test total_weight(MidpointDensity(), p) < analytic(p)
        end
        for sh in (0.5, 1.0, 2.0)
            p = [sh, 2000.0]
            @test total_weight(MidpointDensity(), p) ≈ analytic(p) rtol=1e-4
        end

        # MidpointDensity is exactly the density
        for m in (GammaZero(), GammaExpo())
            p = first(MKtest.dfe_start_points(m))
            a, b = similar(pre.s_grid), similar(pre.s_grid)
            MKtest.fill_dfe_vals!(a, MidpointDensity(), m, p, pre)
            MKtest.dfe_density!(b, m, p, pre.s_grid)
            @test a == b
        end

        r0 = MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero())
        rm = MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero(); quadrature = MidpointDensity())
        rc = MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero(); quadrature = IntervalCDF())
        @test r0.loglik[1] === rm.loglik[1]
        @test isfinite(rc.loglik[1])
    end

    @testset "a user-defined quadrature works end to end" begin
        pre = MKtest.get_precomputed(DFECache(), 40; grid = GrapesGrid())
        v = similar(pre.s_grid)
        MKtest.fill_dfe_vals!(v, TrapezoidDensity(), GammaZero(), [0.3, 2000.0], pre)
        @test all(isfinite, v) && all(>=(0.0), v)
        mass(q, p) = (w = similar(pre.s_grid); MKtest.fill_dfe_vals!(w, q, GammaZero(), p, pre); sum(w .* pre.s_widths))
        p = [0.2, 2000.0]
        exact = MKtest.dfe_cdf_continuous(GammaZero(), p, pre.s_hi[end]) -
                MKtest.dfe_cdf_continuous(GammaZero(), p, pre.s_lo[1])
        @test mass(IntervalCDF(), p) ≈ exact rtol=1e-10
        @test mass(MidpointDensity(), p) < exact
        @test mass(TrapezoidDensity(), p) < exact
        @test mass(TrapezoidDensity(), p) != mass(MidpointDensity(), p)
        r = MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero(); quadrature = TrapezoidDensity())
        @test isfinite(r.loglik[1])
        @test r.loglik[1] != MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero()).loglik[1]
    end

    # --- AIC and LRT core -------------------------------------------------------
    @testset "model comparison: AIC selection vs loglik ordering" begin
        nm = ["simple", "complex"]
        t = MKtest._aic_lrt_tables(nm, [1, 5], [-100.0, -99.0], (i, j) -> i == 1 && j == 2)
        @test t.aic.model[1] == "complex"                 # loglik order
        @test t.aic.model[argmin(t.aic.aic)] == "simple"  # AIC order, the real answer
        @test t.aic.aic == [2*5 - 2*(-99.0), 2*1 - 2*(-100.0)]
        @test nrow(t.lrt) == 1
        @test t.lrt.df[1] == 4
        @test t.lrt.D[1] ≈ 2.0
        @test 0.7 < t.lrt.p_value[1] < 0.8
        t2 = MKtest._aic_lrt_tables(nm, [1, 5], [-100.0, -99.0], (i, j) -> false; extra_cols = (:count => [3, 7],))
        @test eltype(t2.aic.count) == Int
    end

    # --- the output DataFrame ---------------------------------------------------
    @testset "dfe_ml output" begin
        c = DFECache()
        r = MKtest.dfe_ml([ML_SFS, ML_SFS], [ML_DIV, ML_DIV], GammaZero(); cache = c)
        @test nrow(r) == 2 && r.dataset == [1, 2]
        @test length(c) == 1                       # one table for both datasets
        for col in (:model, :n, :converged, :loglik, :aic, :n_params, :free_params, :theta,
                    :negGshape, :negGmean, :syno_orientation_error, :alpha, :alpha_ci_low,
                    :alpha_ci_high, :omega_a, :omega_na, :omega_obs, :alpha_div, :alpha_dfe,
                    :integrator, :kernel, :quadrature, :s_grid, :ne_model, :div_model,
                    :dominance, :demography, :expectation, :optimizer)
            @test col in propertynames(r)
        end
        @test r.n[1] == 40
        @test r.free_params[1] == "negGshape,negGmean,syno_orientation_error"
        @test r.n_params[1] == 3
        @test r.aic[1] ≈ 2 * 3 - 2 * r.loglik[1]
        @test r.loglik[1] == r.loglik[2]           # same data, same answer
        # GRAPES' likelihood CI is on by default (threshold 1.92), and 0 switches it off
        @test r.alpha_ci_low[1] < r.alpha[1] < r.alpha_ci_high[1]
        @test !(:alpha_fww in propertynames(r))    # FWW is MKtest.fww_mk, not part of dfe_ml
        r0ci = MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero(); alpha_ci_threshold = 0, cache = c)
        @test isnan(r0ci.alpha_ci_low[1]) && isnan(r0ci.alpha_ci_high[1])
        # the CI is GRAPES' rule: one lnL threshold drop in E[Dn]; a larger threshold is wider
        rw = MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero(); alpha_ci_threshold = 4.0, cache = c)
        @test rw.alpha_ci_low[1] < r.alpha_ci_low[1] && rw.alpha_ci_high[1] > r.alpha_ci_high[1]

        # fixed parameters keep their column and leave n_params
        rf = MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero(); fixed = (negGshape = 0.3,), cache = c)
        @test rf.negGshape[1] == 0.3 && rf.n_params[1] == 2
        @test_throws ArgumentError MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero(); fixed = (foo = 1.0,))

        # Route A
        ra = MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero(); ne_model = MKtest.AncestralNeRatio(:free), cache = c)
        @test 1.0 <= ra.Ne_ratio[1] <= 50.0
        @test ra.n_params[1] == 4
        @test occursin("Ne_ratio", ra.free_params[1])
        @test_throws ArgumentError MKtest.dfe_ml(ML_SFS, ML_DIV, GammaZero();
            ne_model = MKtest.AncestralNeRatio(:free), div_model = MKtest.DFEConstrainedDivergence())
        @test_throws ArgumentError MKtest.AncestralNeRatio(:fixed)
        @test_throws ArgumentError MKtest.AncestralNeRatio(-1.0)

        # Dominance(0.5) is Additive() exactly
        a = MKtest.get_precomputed(c, 20; dominance = MKtest.Additive())
        b = MKtest.get_precomputed(c, 20; dominance = MKtest.Dominance(0.5))
        @test a !== b && a.expected_counts == b.expected_counts

        # the comparison on a stacked table, and the fitting form
        both = MKtest.dfe_model_comparison(ML_SFS, ML_DIV; models = [GammaZero(), GammaExpo()])
        @test nrow(both.fits) == 2 && nrow(both.aic) == 2
        @test all(both.aic.dataset .== 1)
    end
end

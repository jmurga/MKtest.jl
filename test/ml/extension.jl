# The extension examples of the brief, as written (with the `MKtest.` qualification a
# method added from outside the package needs), one user subtype per axis, and the
# type-stability and allocation checks of the hot path.

module UserDFEs
using MKtest
using Distributions
import MKtest: AbstractDFEModel

# --- the brief's two examples -------------------------------------------------
struct NormalDFE <: AbstractDFEModel end
MKtest.dfe_parameters(::NormalDFE) = (μ = (-1e4, 1e5), σ = (1, 1e5))   # name = (lower, upper)
MKtest.distribution(::NormalDFE, μ, σ) = truncated(Normal(μ, σ), 0, Inf)

struct WithLethal{M} <: AbstractDFEModel
    dfe::M
end
MKtest.dfe_parameters(m::WithLethal) = (MKtest.dfe_parameters(m.dfe)..., p_lethal = (0.0, 0.99))
MKtest.dfe_density!(out, m::WithLethal, p, s) =
    (MKtest.dfe_density!(out, m.dfe, p[1:end-1], s); out .*= 1 - p[end])

# --- one user subtype per axis -----------------------------------------------
struct SameAsPoisson <: MKtest.AbstractLikelihoodKernel end
MKtest.log_kernel(::SameAsPoisson, k, mu, lf, d) = MKtest.log_kernel(MKtest.PoissonKernel(), k, mu, lf, d)

# implements only `integrate`; the hooks' fallbacks are GRAPES' expressions
struct MyRomberg <: MKtest.AbstractIntegrator end
MKtest.integrate(::MyRomberg, f, a::Float64, b::Float64) = MKtest.integrate(MKtest.RombergIntegrator(), f, a, b)

struct UniformGrid <: MKtest.AbstractSGrid end
MKtest.s_grid_points(::UniformGrid) = collect(range(-499.75, 1999.75; length = 5000))
MKtest.s_grid_widths(::UniformGrid, g::Vector{Float64}) = fill(0.5, length(g))

struct RatioTwo <: MKtest.AbstractNeModel end
MKtest.ne_ratio(::RatioTwo) = 2.0

struct PredictedDn <: MKtest.AbstractDivergenceModel end
MKtest.divergence_free(::PredictedDn) = false

struct HalfH <: MKtest.AbstractDominance end
MKtest.effective_s(::HalfH, S::Float64) = 2.0 * 0.5 * S

struct NoRj <: MKtest.AbstractDemography end
MKtest.free_rj(::NoRj) = false

struct WrappedNM <: MKtest.AbstractOptimizer end
MKtest.run_optimizer(::WrappedNM, f, b, x0) = MKtest.run_optimizer(MKtest.NelderMeadOptimizer(), f, b, x0)

struct ThetaTwice <: MKtest.AbstractAlphaSummary end
MKtest.summarize(::ThetaTwice, c::MKtest.FitContext) = (theta_x2 = 2 * c.theta,)

# the PRF sojourn with fitdadi's exterior terms switched on
struct PRFExterior <: MKtest.AbstractSFSExpectation end
MKtest.expectation_counts(::PRFExterior, s, n, integ) = MKtest.expectation_counts(MKtest.PRFSojourn(), s, n, integ)
MKtest.uses_exterior(::PRFExterior) = true
end # module UserDFEs

using .UserDFEs: NormalDFE
using Distributions: truncated, Normal, pdf, cdf, Poisson
using Random: MersenneTwister

# Poisson data from the engine's own expected SFS (theta per site 0.002, n = 40).
function _simulated_sfs(model, p, seed; theta = 0.002, Ln = 1.0e6, Ls = 3.0e5)
    pre = MKtest.get_precomputed(MKtest.DFECache(), 40)
    dens = zeros(length(pre.s_grid))
    MKtest.dfe_density!(dens, model, p, pre.s_grid)
    sel = [sum(dens .* pre.s_widths .* pre.expected_counts[:, j]) for j in 1:39]
    rng = MersenneTwister(seed)
    pn = [rand(rng, Poisson(theta * Ln * sel[j])) for j in 1:39]
    ps = [rand(rng, Poisson(theta * Ls / j)) for j in 1:39]
    return hcat(collect(1.0:39.0), Float64.(pn), Float64.(ps))
end

@testset "dfe_ml extension" begin

    @testset "NormalDFE as written in the brief" begin
        m = NormalDFE()
        # bounds and scales derived from dfe_parameters alone
        lay = MKtest.build_layout(m, MKtest.PoissonKernel(), true)
        @test lay.names == [:μ, :σ, :syno_orientation_error]
        @test lay.bounds[1:2] == [(-1e4, 1e5), (1.0, 1e5)]
        @test lay.log_scale[1:2] == [false, true]
        @test length(MKtest.dfe_start_points(m)) == 8
        @test all(all(lay.bounds[i][1] <= sp[i] <= lay.bounds[i][2] for i in 1:2) for sp in MKtest.dfe_start_points(m))

        # the generic density is the hand-written pdf on |S| (oracle 1)
        p = [100.0, 300.0]
        d = truncated(Normal(100.0, 300.0), 0, Inf)
        sg = collect(range(-2000.0, 50.0; length = 2051))
        out = zeros(length(sg))
        MKtest.dfe_density!(out, m, p, sg)
        @test out == [s < 0 ? pdf(d, -s) : 0.0 for s in sg]

        # IntervalCDF's cell masses are differences of `cdf` (oracle 1)
        pre = MKtest.get_precomputed(MKtest.DFECache(), 40)
        v = similar(pre.s_grid)
        MKtest.fill_dfe_vals!(v, MKtest.IntervalCDF(), m, p, pre)
        mass = v .* pre.s_widths
        neg = findall(pre.s_hi .<= 0)
        @test mass[neg] ≈ [cdf(d, -pre.s_lo[i]) - cdf(d, -pre.s_hi[i]) for i in neg] rtol=1e-10
        @test sum(mass) ≈ cdf(d, 500.0) rtol=1e-10
        # masses and the tail terms come from the distribution too
        @test MKtest.dfe_mass_below(m, p, 50.0) ≈ cdf(d, 50.0)
        @test MKtest.dfe_mass_above(m, p, 50.0) ≈ 1 - cdf(d, 50.0)

        # a fit runs on the shared synthetic data
        r = MKtest.dfe_ml(ML_SFS, ML_DIV, NormalDFE())
        @test r.converged[1] && isfinite(r.loglik[1])
        @test r.model[1] == "NormalDFE"

        # Smoke test, NOT evidence of accuracy: data drawn from the engine's own
        # expectation under NormalDFE(300, 400). The optimizer must reach a
        # likelihood at least as high as the truth's on each of 3 seeds.
        ptrue = [300.0, 400.0]
        for seed in (1, 2, 3)
            sfs = _simulated_sfs(m, ptrue, seed)
            div = [900.0, 1200.0, 1.0e6, 3.0e5]
            rs = MKtest.dfe_ml(sfs, div, m)
            ll_true = maximum(MKtest._loglik_at(m, sfs, div, [ptrue..., oer]) for oer in (0.0, 0.01, 0.05))
            @test rs.loglik[1] >= ll_true - 1e-6
            @info "NormalDFE smoke seed $seed: truth (300, 400), fit ($(round(rs.μ[1]; digits = 1)), $(round(rs.σ[1]; digits = 1))), lnL $(round(rs.loglik[1]; digits = 3)) vs truth $(round(ll_true; digits = 3))"
        end
    end

    @testset "WithLethal as written in the brief" begin
        uw = UserDFEs.WithLethal(MKtest.GammaZero())
        bw = MKtest.WithLethal(MKtest.GammaZero())
        @test MKtest.dfe_param_names(uw) == [:negGshape, :negGmean, :p_lethal]
        # same likelihood as its base at p_lethal = 0, and as the built-in wrapper
        for θ in ([0.2, 1000.0, 0.01], [0.5, 30000.0, 0.0])
            ll_base = MKtest._loglik_at(MKtest.GammaZero(), ML_SFS, ML_DIV, θ)
            # `≈`, not `==`: the user model's default log scale on negGshape sends the
            # value through exp(log(x)), which is not always the identity in floating point
            @test MKtest._loglik_at(uw, ML_SFS, ML_DIV, [θ[1:2]..., 0.0, θ[3]]) ≈ ll_base rtol=1e-12
            for pl in (0.0, 0.3)
                @test MKtest._loglik_at(uw, ML_SFS, ML_DIV, [θ[1:2]..., pl, θ[3]]) ≈
                      MKtest._loglik_at(bw, ML_SFS, ML_DIV, [θ[1:2]..., pl, θ[3]]) rtol=1e-12
            end
        end
        # the fit at p_lethal = 0 equals GammaZero's (its default log scale differs,
        # so the optimizer path differs; same optimum)
        r0 = MKtest.dfe_ml(ML_SFS, ML_DIV, MKtest.GammaZero())
        ru = MKtest.dfe_ml(ML_SFS, ML_DIV, uw; fixed = (p_lethal = 0.0,))
        @test ru.loglik[1] ≈ r0.loglik[1] atol=1e-4
        # free p_lethal fits and stays in its bounds
        rf = MKtest.dfe_ml(ML_SFS, ML_DIV, uw)
        @test isfinite(rf.loglik[1]) && 0.0 <= rf.p_lethal[1] <= 0.99
        @test rf.loglik[1] >= ru.loglik[1] - 1e-6
    end

    @testset "one user subtype per axis" begin
        g0 = MKtest.GammaZero()
        r0 = MKtest.dfe_ml(ML_SFS, ML_DIV, g0)

        # kernel
        rk = MKtest.dfe_ml(ML_SFS, ML_DIV, g0; kernel = UserDFEs.SameAsPoisson())
        @test rk.loglik[1] == r0.loglik[1]

        # integrator: only `integrate` defined, the table equals Romberg's
        c = MKtest.DFECache()
        a = MKtest.get_precomputed(c, 20; integrator = UserDFEs.MyRomberg())
        b = MKtest.get_precomputed(c, 20; integrator = MKtest.RombergIntegrator())
        @test a.expected_counts == b.expected_counts

        # S grid
        rg = MKtest.dfe_ml(ML_SFS, ML_DIV, g0; s_grid = UserDFEs.UniformGrid())
        # a different grid is a different quadrature: on this data negGshape sits at its
        # 0.05 bound, a density far too spiky at S -> 0- for 0.5-wide cells, so alpha
        # moves; the check is that the grid reaches the fit
        @test isfinite(rg.loglik[1]) && rg.converged[1]
        @test rg.loglik[1] != r0.loglik[1]

        # Ne model and divergence model
        rn = MKtest.dfe_ml(ML_SFS, ML_DIV, MKtest.GammaExpo(); ne_model = UserDFEs.RatioTwo())
        rn2 = MKtest.dfe_ml(ML_SFS, ML_DIV, MKtest.GammaExpo(); ne_model = MKtest.AncestralNeRatio(2.0))
        @test rn.alpha[1] == rn2.alpha[1]
        rd = MKtest.dfe_ml(ML_SFS, ML_DIV, g0; div_model = UserDFEs.PredictedDn())
        rd2 = MKtest.dfe_ml(ML_SFS, ML_DIV, g0; div_model = MKtest.DFEConstrainedDivergence())
        @test rd.loglik[1] == rd2.loglik[1]

        # dominance
        h1 = MKtest.get_precomputed(c, 20; dominance = UserDFEs.HalfH())
        h2 = MKtest.get_precomputed(c, 20; dominance = MKtest.Additive())
        @test h1.expected_counts == h2.expected_counts

        # demography
        re = MKtest.dfe_ml(ML_SFS, ML_DIV, g0; demography = UserDFEs.NoRj())
        re2 = MKtest.dfe_ml(ML_SFS, ML_DIV, g0; demography = MKtest.Equilibrium())
        @test re.loglik[1] == re2.loglik[1]
        @test re.loglik[1] <= r0.loglik[1] + 1e-8  # fewer free r_j cannot fit better

        # optimizer
        ro = MKtest.dfe_ml(ML_SFS, ML_DIV, g0; optimizer = UserDFEs.WrappedNM())
        ro2 = MKtest.dfe_ml(ML_SFS, ML_DIV, g0; optimizer = MKtest.NelderMeadOptimizer())
        @test ro.loglik[1] == ro2.loglik[1]
        @test abs(ro.loglik[1] - r0.loglik[1]) < 1e-2

        # summaries
        rs = MKtest.dfe_ml(ML_SFS, ML_DIV, g0; summaries = (UserDFEs.ThetaTwice(),))
        @test rs.theta_x2[1] == 2 * rs.theta[1]
        @test !(:alpha in propertynames(rs))       # only the summaries asked for

        # expectation with exterior terms: lethal_mass is the mass beyond the outer
        # cell EDGE (dev evaluated it at the outer grid point)
        lg = MKtest.LogGammaGrid(; Npts = 200)
        rx = MKtest.dfe_ml(ML_SFS, ML_DIV, g0; expectation = UserDFEs.PRFExterior(), s_grid = lg)
        pre = MKtest.get_precomputed(MKtest.DFECache(), 40; grid = lg, expectation = UserDFEs.PRFExterior())
        p = [rx.negGshape[1], rx.negGmean[1]]
        @test rx.lethal_mass[1] == MKtest.dfe_mass_above(g0, p, pre.s_ext_hi)
        @test pre.s_ext_hi > abs(pre.s_grid[1])
    end

    # --- type stability and allocations of the hot path -------------------------
    # Measured 2026-10-15 on the laptop (Julia 1.13.1, 1 thread, src sha in the handoff):
    #   objective obj(x), GammaZero / GammaExpo / NormalDFE ............ 0 bytes
    #   fill_dfe_vals!, MidpointDensity and IntervalCDF (GammaZero) .... 0 bytes
    #   generic dfe_density! through `distribution` (NormalDFE) ........ 0 bytes
    #   expected_count, ExactSojourn / GaussLegendreSojourn ............ 0 bytes
    #   integrate(f), GaussLegendreSojourn / ExactSojourn .............. 0 bytes
    #     (160 bytes before `integrate` was made to specialise on `f`)
    #   expected_count / integrate, RombergIntegrator ......... 3776 / 2240 bytes
    #     (GRAPES' transliteration, kept as is; it only runs at table build)
    @testset "hot path: @inferred and @allocated" begin
        alloc_call(f, x) = (f(x); @allocated f(x))
        for m in (MKtest.GammaZero(), MKtest.GammaExpo(), NormalDFE())
            obj = MKtest._objective_for(m, ML_SFS, ML_DIV)
            lay = obj.layout
            full = vcat(first(MKtest.dfe_start_points(m)), 0.01)
            x = [MKtest._to_opt(full[i], lay.log_scale[i]) for i in lay.opt_idx]
            @test (@inferred obj(x)) isa Float64
            @test alloc_call(obj, x) == 0
        end

        pre = MKtest.get_precomputed(MKtest.DFECache(), 40)
        v = similar(pre.s_grid)
        scratch = similar(pre.s_edges)
        p = [0.2, 2000.0]
        @inferred MKtest.fill_dfe_vals!(v, MKtest.MidpointDensity(), MKtest.GammaZero(), p, pre, scratch)
        @inferred MKtest.fill_dfe_vals!(v, MKtest.IntervalCDF(), MKtest.GammaZero(), p, pre, scratch)
        fdv(q) = (MKtest.fill_dfe_vals!(v, q, MKtest.GammaZero(), p, pre, scratch);
                  @allocated MKtest.fill_dfe_vals!(v, q, MKtest.GammaZero(), p, pre, scratch))
        @test fdv(MKtest.MidpointDensity()) == 0
        @test fdv(MKtest.IntervalCDF()) == 0

        # the generic `distribution` path is type stable and allocation free
        sg = pre.s_grid
        pn = [100.0, 300.0]
        @inferred MKtest.dfe_density!(v, NormalDFE(), pn, sg)
        dd() = (MKtest.dfe_density!(v, NormalDFE(), pn, sg); @allocated MKtest.dfe_density!(v, NormalDFE(), pn, sg))
        @test dd() == 0

        lc = log(binomial(BigInt(40), 3))
        lcf = Float64(lc)
        for integ in (MKtest.ExactSojourn(), MKtest.GaussLegendreSojourn(), MKtest.RombergIntegrator())
            @test (@inferred MKtest.expected_count(integ, -5.0, 3, 40, lcf)) isa Float64
        end
        ec(integ) = (MKtest.expected_count(integ, -5.0, 3, 40, lcf); @allocated MKtest.expected_count(integ, -5.0, 3, 40, lcf))
        @test ec(MKtest.ExactSojourn()) == 0
        @test ec(MKtest.GaussLegendreSojourn()) == 0
        f = x -> x^2
        gi(g) = (MKtest.integrate(g, f, 0.0, 1.0); @allocated MKtest.integrate(g, f, 0.0, 1.0))
        @test gi(MKtest.GaussLegendreSojourn()) == 0
        @test gi(MKtest.ExactSojourn()) == 0
        @info "hot-path allocations (bytes)" objective_GammaZero = alloc_call(MKtest._objective_for(MKtest.GammaZero(), ML_SFS, ML_DIV), [0.2, log(2000.0), 0.01]) romberg_expected_count = ec(MKtest.RombergIntegrator())
    end
end

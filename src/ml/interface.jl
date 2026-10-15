# =============================================================================
# dfe_ml engine: the axes and their interfaces
#
# `dfe_ml` is the GRAPES model (Galtier 2016) written as a set of independent
# axes. Each axis is an abstract type with a small interface; `dfe_ml` takes one
# instance per axis as a keyword and calls only the interface methods, so a new
# subtype of one axis needs no edit to any other. The defaults reproduce GRAPES,
# except the integrator (`ExactSojourn`, see its docstring).
#
#   axis                     abstract type              default
#   SFS expectation          AbstractSFSExpectation     PRFSojourn()
#   dominance                AbstractDominance          Additive()
#   demography (SFS shape)   AbstractDemography         NuisanceRj()
#   ancestral Ne             AbstractNeModel            ConstantNe()
#   DFE                      AbstractDFEModel           (required argument)
#   likelihood kernel        AbstractLikelihoodKernel   PoissonKernel()
#   role of divergence       AbstractDivergenceModel    SaturatedDivergence()
#   sojourn integrals        AbstractIntegrator         ExactSojourn()
#   DFE weight per S cell    AbstractDFEQuadrature      MidpointDensity()
#   S grid                   AbstractSGrid              GrapesGrid()
#   optimizer                AbstractOptimizer          LBFGSOptimizer()
#   reported summaries       AbstractAlphaSummary       (DFEMAlpha(), FWWAlpha())
#
# Selection coefficients are S = 4 Ne s (GRAPES' convention) throughout.
# =============================================================================

# Constants, GRAPES' where it has one.
const QUASI_ZERO = 1e-4
const LARGE_S_VALUE = 500.0
const VERY_LARGE_S_VALUE = 10000.0
const GAMMA_INTEGRAL_BREAK = 0.5
const S_GRID_K1 = 1000
const S_GRID_K2 = 1000
const SCALED_BETA_BOUNDARY = 25.0

# -----------------------------------------------------------------------------
# SFS expectation
# -----------------------------------------------------------------------------

"""
    AbstractSFSExpectation

Who produces the expected SFS of a mutation with selection coefficient S, at
theta = 1. Everything downstream (the DFE integral, theta, the likelihood, alpha)
is shared, so a new producer (a diffusion solve under a fitted demography, a
Wright-Fisher transition matrix) plugs in here.

A subtype implements:
- `expectation_counts(e, s_eff, n, integ) -> Matrix{Float64}` (`length(s_eff)` x `n-1`)

and may override:
- `neutral_shape(e, n) -> Vector{Float64}`: per-bin factor of the S = 0 spectrum
  relative to `1/j` (all ones by default);
- `uses_exterior(e) -> Bool`: add the analytic mass outside the S grid (fitdadi's
  exterior terms). `false` by default, as in GRAPES;
- `default_s_grid(e)`, `default_demography(e)`.

The value of `e` is part of the `DFECache` key, so a subtype must be immutable or
define `==` and `hash`.
"""
abstract type AbstractSFSExpectation end

"""
    PRFSojourn()

The closed-form Sawyer-Hartl Poisson random field sojourn density at equilibrium,
GRAPES' model. Demographic distortion of the SFS is absorbed by the per-bin `r_j`
nuisance parameters (`NuisanceRj()`), not modelled.
"""
struct PRFSojourn <: AbstractSFSExpectation end

"""
    expectation_counts(e, s_eff, n, integ) -> Matrix{Float64}

Expected SFS bin `j` (columns, `j = 1..n-1`) for each effective selection
coefficient `s_eff[i]` (rows), at theta = 1, integrated with `integ`.
"""
function expectation_counts end

"""
    neutral_shape(e, n) -> Vector{Float64}

Per-bin factor by which the S = 0 expected SFS of `e` departs from the
equilibrium `1/j`. All ones for `PRFSojourn`. Multiplying by exactly 1.0 is exact,
so the default path is unchanged by this hook.
"""
neutral_shape(::AbstractSFSExpectation, n::Int) = ones(n - 1)

"""
    uses_exterior(e) -> Bool

Whether the DFE integral adds two analytic terms for the mass outside the S grid:
mass below the least deleterious cell, applied to the neutral spectrum, and mass
beyond the most deleterious cell, applied to that cell's row (fitdadi). `false` for
`PRFSojourn`, which reproduces GRAPES, and whose `GrapesGrid` tiles |S| down to 0.
"""
uses_exterior(::AbstractSFSExpectation) = false

"""
    default_s_grid(e) -> AbstractSGrid

The S grid used when `dfe_ml` gets `s_grid = nothing`. `GrapesGrid()` by default.
"""
default_s_grid(::AbstractSFSExpectation) = GrapesGrid()

"""
    default_demography(e) -> AbstractDemography

The demography axis used when `dfe_ml` gets `demography = nothing`. `NuisanceRj()`
by default. A producer that carries a demographic history itself should return
`Equilibrium()`: free `r_j` on top of a modelled history count it twice.
"""
default_demography(::AbstractSFSExpectation) = NuisanceRj()

# -----------------------------------------------------------------------------
# Dominance
# -----------------------------------------------------------------------------

"""
    AbstractDominance

Dominance of new mutations, through the fast heterozygote approximation: a mutation
with selection coefficient S and dominance h is treated as additive with
`S_eff = 2 h S` (exact for rare alleles, the additive model at h = 0.5). The DFE is
still a density over the real S; only the sojourn and fixation use `S_eff`.

A subtype implements `effective_s(d, S) -> Float64`. Its value is part of the
`DFECache` key.
"""
abstract type AbstractDominance end

"""
    Additive()

h = 0.5, `S_eff = S`. The default, and GRAPES' model.
"""
struct Additive <: AbstractDominance end

"""
    Dominance(h)

A fixed dominance coefficient: `S_eff = 2 h S`. `Dominance(0.5)` equals `Additive()`
exactly (`2 * 0.5 * S == S` in floating point).
"""
struct Dominance <: AbstractDominance
    h::Float64
end

"""
    DominanceFunction(f)

Dominance as a function of S, `h = f(S)`: `S_eff = 2 f(S) S`. The cache key uses the
value of `f`, so two different anonymous functions never share a cache entry, and
the same function object always does.
"""
struct DominanceFunction{F} <: AbstractDominance
    f::F
end

"""
    effective_s(d::AbstractDominance, S) -> Float64

Effective selection coefficient of a mutation with coefficient `S` under `d`.
"""
function effective_s end

@inline effective_s(::Additive, S::Float64) = S
@inline effective_s(d::Dominance, S::Float64) = 2.0 * d.h * S
@inline effective_s(d::DominanceFunction, S::Float64) = 2.0 * Float64(d.f(S)) * S

# -----------------------------------------------------------------------------
# Demography (how the SFS shape is treated)
# -----------------------------------------------------------------------------

"""
    AbstractDemography

How demographic distortion of the SFS enters the likelihood. GRAPES estimates one
nuisance factor `r_j` per frequency bin, shared by synonymous and nonsynonymous
sites (Eyre-Walker et al. 2006). This axis does NOT fit a demographic history.

A subtype implements `free_rj(d) -> Bool`.
"""
abstract type AbstractDemography end

"""
    NuisanceRj()

The `r_j` are profiled out of the likelihood at every evaluation (GRAPES). Default.
"""
struct NuisanceRj <: AbstractDemography end

"""
    Equilibrium()

`r_j = 1`: the population is assumed to be at equilibrium, which removes `n - 2`
nuisance parameters.
"""
struct Equilibrium <: AbstractDemography end

"""
    free_rj(d::AbstractDemography) -> Bool

Whether the per-bin `r_j` are estimated.
"""
function free_rj end

free_rj(::NuisanceRj) = true
free_rj(::Equilibrium) = false

# -----------------------------------------------------------------------------
# Ancestral Ne
# -----------------------------------------------------------------------------

"""
    AbstractNeModel

Ratio `r = Ne_ancestral / Ne_current` of the lineage along which divergence
accumulated. The DFE fitted to the SFS is rescaled by `r` (S-scale parameters times
`r`, see `dfe_ne_multipliers`) before it predicts divergence, as GRAPES'
`-anc_to_rec_Ne_ratio` does.

A subtype implements `ne_ratio(m) -> Float64`, or `:free` to estimate the ratio
from the divergence counts.
"""
abstract type AbstractNeModel end

"""
    ConstantNe()

`r = 1`. Default.
"""
struct ConstantNe <: AbstractNeModel end

"""
    AncestralNeRatio(r)
    AncestralNeRatio(:free)

A fixed ratio `r > 0`, or `:free`: a two-stage estimate (Route A). Stage 1 fits the
DFE with `r = 1`; stage 2 finds the `r` in [1, 50] that maximises the likelihood of
(Dn, Ds) given that DFE (Brent search on log r); stage 3 refits at that `r`. The
result gains a column `Ne_ratio` and one parameter in `n_params`. The synonymous SFS
is not used to fit a demography.
"""
struct AncestralNeRatio <: AbstractNeModel
    ratio::Union{Float64,Symbol}
    function AncestralNeRatio(r::Union{Real,Symbol})
        if r isa Symbol
            r === :free || throw(ArgumentError("AncestralNeRatio: the only Symbol is :free, got :$r"))
            return new(r)
        end
        r > 0 || throw(ArgumentError("AncestralNeRatio: the ratio must be positive, got $r"))
        return new(Float64(r))
    end
end

"""
    ne_ratio(m::AbstractNeModel) -> Union{Float64,Symbol}

The ancestral-to-current Ne ratio, or `:free`.
"""
function ne_ratio end

ne_ratio(::ConstantNe) = 1.0
ne_ratio(m::AncestralNeRatio) = m.ratio

# -----------------------------------------------------------------------------
# DFE models
# -----------------------------------------------------------------------------

"""
    AbstractDFEModel

A distribution of fitness effects over S. A new DFE implements

    MKtest.dfe_parameters(m) -> NamedTuple    # name = (lower, upper), in parameter order

and EITHER

    MKtest.distribution(m, p...)              # a Distributions.jl law of |S|, S < 0

which gives the density, the tail masses and the CDF for free, OR

    MKtest.dfe_density!(out, m, p, s_grid)    # density at each S of `s_grid`, into `out`

(wrappers such as `WithLethal` use the second form). Optional overrides:
`dfe_log_scale`, `dfe_start_points`, `dfe_name`, `discrete_components`,
`dfe_mass_below`, `dfe_mass_above`, `dfe_cdf_continuous!`, `dfe_ne_multipliers`.
See the "How to add a new DFE" tutorial.
"""
abstract type AbstractDFEModel end

"""
    dfe_parameters(m::AbstractDFEModel) -> NamedTuple

The free parameters of `m`, as `name = (lower, upper)` in parameter order. The
vector `p` handed to `dfe_density!` and `distribution` follows this order.
"""
function dfe_parameters end

"""
    distribution(m::AbstractDFEModel, p...) -> Distributions.UnivariateDistribution

The law of |S| for deleterious mutations (S < 0), given the parameter values in
`dfe_parameters` order. Defining it is enough for `dfe_density!`, the tail masses
and the CDF (for `IntervalCDF`).
"""
function distribution end

"""
    dfe_density!(out, m, p, s_grid) -> out

Density of the DFE `m` with parameters `p` at every S of `s_grid`, written into `out`.
The generic method uses `distribution`: `pdf(d, -S)` for S < 0 and 0 otherwise.
"""
function dfe_density! end

# -----------------------------------------------------------------------------
# Likelihood kernel
# -----------------------------------------------------------------------------

"""
    AbstractLikelihoodKernel

The law each observed count (every SFS bin, synonymous and nonsynonymous, and each
divergence total) is scored with, given its expectation.

A subtype implements

    log_kernel(k, count, mu, logfact_count, disp) -> Float64

where `logfact_count = loggamma(count + 1)` and `disp` is the kernel's own fitted
parameter (`NaN` when it has none). A kernel with a parameter also defines
`kernel_parameters(k) -> (name = (lower, upper),)` (one parameter at most) and
`kernel_start(k) -> Vector{Float64}`.
"""
abstract type AbstractLikelihoodKernel end

"""
    PoissonKernel()

`k log(mu) - mu - log(k!)`. GRAPES, DFE-alpha, fitdadi and fastDFE. Default.
"""
struct PoissonKernel <: AbstractLikelihoodKernel end

"""
    NegativeBinomialKernel()

Negative binomial with mean `mu` and one shared dispersion `a` (`Var = mu + mu^2/a`),
polyDFE's kernel (Tataru et al. 2017). `a` is estimated on a log scale in
`[1e-4, 1e6]`; the upper bound is wider than polyDFE's 1000 because on human data
that bound binds. `a -> Inf` is the Poisson limit.
"""
struct NegativeBinomialKernel <: AbstractLikelihoodKernel end

"""
    kernel_parameters(k) -> NamedTuple

The kernel's own free parameters as `name = (lower, upper)`. Empty by default.
"""
kernel_parameters(::AbstractLikelihoodKernel) = NamedTuple()
kernel_parameters(::NegativeBinomialKernel) = (a = (QUASI_ZERO, 1.0e6),)

"""
    kernel_start(k) -> Vector{Float64}

Starting values of the kernel parameters. Empty by default.
"""
kernel_start(::AbstractLikelihoodKernel) = Float64[]
kernel_start(::NegativeBinomialKernel) = [10.0]

"""
    log_kernel(kernel, count, mu, logfact_count, disp) -> Float64

Log-probability of `count` under `kernel` with mean `mu`.
"""
log_kernel(::PoissonKernel, k::Real, mu::Float64, logfact::Float64, ::Float64) =
    k * log(mu) - mu - logfact

# polyDFE 2.0 `src/likelihood.c:get_log_neg_bin`, including its k == 0 branch.
function log_kernel(::NegativeBinomialKernel, k::Real, mu::Float64, logfact::Float64, a::Float64)
    k == 0 && return a * (log(a) - log(a + mu))
    return loggamma(k + a) - loggamma(a) - logfact +
           a * log(a) + k * log(mu) - (a + k) * log(a + mu)
end

# -----------------------------------------------------------------------------
# Role of divergence
# -----------------------------------------------------------------------------

"""
    AbstractDivergenceModel

What (Dn, Ds) contribute to the likelihood. A subtype implements
`uses_divergence(d) -> Bool` and `divergence_free(d) -> Bool`.
"""
abstract type AbstractDivergenceModel end

"""
    SaturatedDivergence()

Dn and Ds are fitted exactly (one free rate each); the DFE is fitted to the SFS and
alpha comes from the observed Dn. GRAPES' default. Default.
"""
struct SaturatedDivergence <: AbstractDivergenceModel end

"""
    DFEConstrainedDivergence()

Dn is predicted by the DFE and scored (GRAPES `-no_div_param`, the "[-A]" model of
Galtier 2016).
"""
struct DFEConstrainedDivergence <: AbstractDivergenceModel end

"""
    NoDivergence()

Divergence is not used (GRAPES `-no_div_data`); alpha is the DFE's prediction.
"""
struct NoDivergence <: AbstractDivergenceModel end

"""
    uses_divergence(d) -> Bool

Whether divergence enters the likelihood at all.
"""
function uses_divergence end
uses_divergence(::AbstractDivergenceModel) = true
uses_divergence(::NoDivergence) = false

"""
    divergence_free(d) -> Bool

Whether Dn is a free (saturated) quantity rather than a DFE prediction.
"""
function divergence_free end
divergence_free(::AbstractDivergenceModel) = true
divergence_free(::DFEConstrainedDivergence) = false

# -----------------------------------------------------------------------------
# Integrators
# -----------------------------------------------------------------------------

"""
    AbstractIntegrator

How the two sojourn integrals the engine needs are computed: the expected SFS bin
`j` at S, and GRAPES' dN correction term. A subtype implements the generic rule

    integrate(i, f, a::Float64, b::Float64) -> Float64

and may take over either integral with the hooks `expected_count(i, S, j, n, logc)`
and `dn_correction_term(i, S, n)`. Its value is part of the `DFECache` key.
"""
abstract type AbstractIntegrator end

"""
    integrate(i::AbstractIntegrator, f, a, b) -> Float64

Integral of `f` over `[a, b]`.
"""
function integrate end

# -----------------------------------------------------------------------------
# DFE quadrature
# -----------------------------------------------------------------------------

"""
    AbstractDFEQuadrature

How the DFE weight of each S cell is formed. Every integral has the form
`sum_i dfe_vals[i] * s_widths[i] * f(S_i)`, and a quadrature fills `dfe_vals` with
the cell's average density. A subtype implements

    fill_dfe_vals!(out, q, m::AbstractDFEModel, p, pre::PrecomputedData) -> out

(`pre.s_grid`, `pre.s_widths`, `pre.s_lo`, `pre.s_hi` are the cells).
"""
abstract type AbstractDFEQuadrature end

"""
    MidpointDensity()

`dfe_vals[i] = density(S_i)`, GRAPES' rule. Default.
"""
struct MidpointDensity <: AbstractDFEQuadrature end

"""
    IntervalCDF()

`dfe_vals[i] = (CDF(hi_i) - CDF(lo_i)) / width_i`: the exact mass of each cell
(fastDFE's rule). Better than the midpoint where the density is singular (a gamma
with shape < 1 at S -> 0-), at the cost of one CDF evaluation per cell edge.
"""
struct IntervalCDF <: AbstractDFEQuadrature end

# -----------------------------------------------------------------------------
# S grid
# -----------------------------------------------------------------------------

"""
    AbstractSGrid

The grid the DFE integral is discretised on. A subtype implements
`s_grid_points(g) -> Vector{Float64}` (ascending) and
`s_grid_widths(g, grid) -> Vector{Float64}`, and may override
`s_cell_bounds(g, grid, widths) -> (lo, hi)` (default: `S_i -+ w_i/2`). Its value is
part of the `DFECache` key.
"""
abstract type AbstractSGrid end

"""
    s_grid_points(g) -> Vector{Float64}

The S values the DFE is evaluated at, ascending (most deleterious first).
"""
function s_grid_points end

"""
    s_grid_widths(g, grid) -> Vector{Float64}

Quadrature weight (cell width) of each grid point.
"""
function s_grid_widths end

"""
    GrapesGrid()

GRAPES' grid: 4000 points spaced linearly in four blocks, [-500, -0.5], [-0.5, 0.5],
[0.5, 500], [500, 10000], with midpoint cells. Default.
"""
struct GrapesGrid <: AbstractSGrid end

"""
    LogGammaGrid(; int_bounds = (1e-4, 1000.0), Npts = 1000)

fitdadi's convention: log-spaced in |gamma| over `int_bounds` (S = 2 gamma), with a
log-spaced beneficial side up to S = 10000.
"""
struct LogGammaGrid <: AbstractSGrid
    int_bounds::Tuple{Float64,Float64}
    Npts::Int
end
function LogGammaGrid(; int_bounds = (1e-4, 1000.0), Npts::Integer = 1000)
    ib = (Float64(int_bounds[1]), Float64(int_bounds[2]))
    0 < ib[1] < ib[2] || throw(ArgumentError("int_bounds must be 0 < lo < hi in |gamma|, got $ib"))
    return LogGammaGrid(ib, Int(Npts))
end

# -----------------------------------------------------------------------------
# Optimizers
# -----------------------------------------------------------------------------

"""
    AbstractOptimizer

The minimiser of the negative log-likelihood. A subtype implements

    run_optimizer(o, f, bounds::Vector{Tuple{Float64,Float64}}, x0::Vector{Float64})
        -> (minimizer, minimum, converged)

`bounds` and `x0` are in optimizer space (log scale where a parameter is log-scaled).
"""
abstract type AbstractOptimizer end

"""
    LBFGSOptimizer(; max_iter = 10_000, f_reltol = 1e-8)

Box-constrained L-BFGS (`Optim.Fminbox(LBFGS(BackTracking()))`) with finite-difference
gradients. Default.
"""
struct LBFGSOptimizer <: AbstractOptimizer
    max_iter::Int
    f_reltol::Float64
end
LBFGSOptimizer(; max_iter::Int = 10_000, f_reltol::Float64 = 1e-8) = LBFGSOptimizer(max_iter, f_reltol)

"""
    NelderMeadOptimizer(; iterations = 3000)

Derivative-free simplex; the bounds are honoured by a logit reparametrisation.
"""
struct NelderMeadOptimizer <: AbstractOptimizer
    iterations::Int
end
NelderMeadOptimizer(; iterations::Int = 3000) = NelderMeadOptimizer(iterations)

"""
    run_optimizer(o::AbstractOptimizer, f, bounds, x0) -> (minimizer, minimum, converged)

Minimise `f` from `x0` inside `bounds`.
"""
function run_optimizer end

# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------

"""
    AbstractAlphaSummary

A set of output columns computed from a finished fit. A subtype implements

    summarize(s, ctx::FitContext) -> NamedTuple

whose fields become columns of the `dfe_ml` DataFrame. `FitContext` documents what
a fit exposes.
"""
abstract type AbstractAlphaSummary end

"""
    DFEMAlpha()

GRAPES' DFE-model alpha: `alpha = (Dn - E[Dn | S < nearly_neutral]) / Dn`, with its
likelihood CI (when `alpha_ci_threshold > 0`), omega_a, omega_na, the observed
omega, and polyDFE's two readings `alpha_div` (observed Dn as denominator) and
`alpha_dfe` (DFE-predicted Dn as denominator).
"""
struct DFEMAlpha <: AbstractAlphaSummary end

"""
    FWWAlpha()

The Fay-Wyckoff-Wu alpha, from SFS bins at frequency >= `FWW_threshold`, as GRAPES
reports it. NaN when undefined.
"""
struct FWWAlpha <: AbstractAlphaSummary end

"""
    summarize(s::AbstractAlphaSummary, ctx::FitContext) -> NamedTuple

The columns `s` adds to a `dfe_ml` row.
"""
function summarize end

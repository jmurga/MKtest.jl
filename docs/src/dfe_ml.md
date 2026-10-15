# DFE and alpha by maximum likelihood (`dfe_ml`)

`MKtest.dfe_ml` fits a distribution of fitness effects (DFE) and the rate of adaptive
evolution to the site frequency spectrum (SFS) and divergence of synonymous and
nonsynonymous sites. It is the GRAPES model (Galtier 2016) in pure Julia: the
nonsynonymous SFS is the DFE-weighted Poisson random field sojourn, demographic distortion
is absorbed by one nuisance factor `r_j` per frequency bin shared by both classes, and
alpha is the share of nonsynonymous substitutions not predicted by mutations with
`S < nearly_neutral`. Selection coefficients are `S = 4 Ne s`.

```julia
using MKtest
sfs, div, _ = MKtest.parse_sfs("data/MKdata_tgp_may2023.txt"; sample_size = 661)
sfs20 = MKtest.project(sfs, 20)                       # 20 individuals, 40 chromosomes

fit = MKtest.dfe_ml(sfs20, div, MKtest.GammaExpo())    # a DataFrame, one row
fit.alpha, fit.negGshape, fit.negGmean
```

`sfs` is the matrix `[DAC Pn Ps]` with every bin `1..n-1` (as `parse_sfs` returns it,
not cumulated) and `div` is `[Dn, Ds, Ln, Ls]`; a Vector of each fits several datasets,
one row each.

## Output

One row per dataset with:

| columns | content |
|---|---|
| `dataset`, `model`, `n` | dataset index, model name, chromosomes |
| `converged`, `loglik`, `aic`, `n_params`, `free_params`, `theta` | the fit; `free_params` lists the optimised parameters, comma separated |
| one column per parameter | the DFE's (`dfe_parameters` order), then `syno_orientation_error`, then the kernel's (`a`) |
| `Ne_ratio` | only under `AncestralNeRatio(:free)` |
| `alpha`, `alpha_ci_low`, `alpha_ci_high`, `omega_a`, `omega_a_ci_*`, `omega_na`, `omega_na_ci_*`, `omega_obs`, `alpha_div`, `alpha_dfe` | `DFEMAlpha()`; the CI is GRAPES' likelihood interval, `alpha_ci_threshold = 1.92` by default (0 disables it) |
| `expectation`, `integrator`, `quadrature`, `s_grid`, `dominance`, `demography`, `ne_model`, `div_model`, `kernel`, `optimizer` | the axis settings of the fit |

`alpha_div` uses the observed Dn as denominator (it equals `alpha` under the default
`SaturatedDivergence`); `alpha_dfe` uses the Dn the fitted DFE predicts. They agree when
the DFE explains the observed divergence. The data-only MK alphas (standard, FWW,
asymptotic) are separate methods: `MKtest.standard_mk`, `MKtest.fww_mk`,
`MKtest.asymptotic_mk`.

## Model comparison

```julia
cmp = MKtest.dfe_model_comparison(sfs20, div; models = [MKtest.GammaZero(), MKtest.GammaExpo()])
cmp.aic        # model, n_params, loglik, aic, delta_aic, weight, converged, alpha, omega_a
cmp.lrt        # model_null, model_alt, loglik_null, loglik_alt, df, D, p_value
cmp.fits       # the stacked dfe_ml rows

MKtest.dfe_model_comparison(vcat(fit1, fit2; cols = :union))   # on existing fits
```

Both tables carry `dataset` and are sorted by log-likelihood within a dataset; choose the
AIC-best model with `argmin(aic)`. A model is nested in another when its `free_params` are
a strict subset (GammaZero in GammaExpo, ExpoZero in GammaZero, GammaZero in
`WithLethal(GammaZero())`). Non-converged fits stay in the AIC table but not in the LRTs.

## The axes

Each keyword of `dfe_ml` is one independent axis: an abstract type whose interface the
engine calls, so a new subtype of one axis needs no change to any other.

| keyword | abstract type | default | others | a new subtype implements |
|---|---|---|---|---|
| `model` (argument) | `AbstractDFEModel` | -- | `GammaZero`, `GammaExpo`, `GammaGamma`, `ExpoZero`, `ExpoExpo`, `ExpoGamma`, `DisplGamma`, `ScaledBeta`, `LogNormalZero`, `NeuGamma`, `FGMBesselK`, `FixedSpikes`, `WithLethal(m)` | `dfe_parameters` + `distribution`, or `dfe_parameters` + `dfe_density!` ([tutorial](new_dfe.md)) |
| `expectation` | `AbstractSFSExpectation` | `PRFSojourn()` | -- | `expectation_counts(e, s_eff, n, integ)` |
| `dominance` | `AbstractDominance` | `Additive()` | `Dominance(h)`, `DominanceFunction(f)` | `effective_s(d, S)` |
| `demography` | `AbstractDemography` | `NuisanceRj()` | `Equilibrium()` | `free_rj(d)` |
| `ne_model` | `AbstractNeModel` | `ConstantNe()` | `AncestralNeRatio(r)`, `AncestralNeRatio(:free)` | `ne_ratio(m)` |
| `div_model` | `AbstractDivergenceModel` | `SaturatedDivergence()` | `DFEConstrainedDivergence()`, `NoDivergence()` | `uses_divergence(d)`, `divergence_free(d)` |
| `kernel` | `AbstractLikelihoodKernel` | `PoissonKernel()` | `NegativeBinomialKernel()` | `log_kernel(k, count, mu, logfact, disp)` |
| `integrator` | `AbstractIntegrator` | `ExactSojourn()` | `RombergIntegrator()`, `GaussLegendreSojourn()` | `integrate(i, f, a, b)` |
| `quadrature` | `AbstractDFEQuadrature` | `MidpointDensity()` | `IntervalCDF()` | `fill_dfe_vals!(out, q, m, p, pre)` |
| `s_grid` | `AbstractSGrid` | `GrapesGrid()` | `LogGammaGrid()` | `s_grid_points(g)`, `s_grid_widths(g, grid)` |
| `optimizer` | `AbstractOptimizer` | `LBFGSOptimizer()` | `NelderMeadOptimizer()` | `run_optimizer(o, f, bounds, x0)` |
| `summaries` | `AbstractAlphaSummary` | `(DFEMAlpha(),)` | -- | `summarize(s, ctx::FitContext)` |

Methods added from outside the package are written with the module prefix, e.g.
`MKtest.log_kernel(::MyKernel, k, mu, lf, d) = ...`.

GRAPES' options map to: `-no_div_param` = `div_model = DFEConstrainedDivergence()`,
`-no_div_data` = `NoDivergence()`, `-anc_to_rec_Ne_ratio r` = `ne_model = AncestralNeRatio(r)`,
`-fixed_param` = `fixed = (name = value,)`, `-no_syn_orient_error` = `no_syn_orient_error = true`,
`-fold` = `fold = true`, `-nb_rand_start` = `nb_rand_start`.

## Integrators and GRAPES parity

The default `ExactSojourn()` evaluates the expected SFS in closed form (Kummer's
confluent hypergeometric function), to about 1e-16 against a BigFloat quadrature.
GRAPES integrates with Romberg over `[1e-4, 1 - 1e-4]`, which loses mass from the rare
bins: bin 1 is low by 0.19% at n = 20, 0.99% at n = 100 and 6.4% at n = 661. To compare
with the GRAPES binary, pass `integrator = RombergIntegrator()`; the parity tests do. The
GRAPES binary itself is wrapped by `MKtest.grapes`.

## Fixing and bounding parameters

```julia
MKtest.dfe_ml(sfs20, div, MKtest.GammaExpo(); fixed = (pos_prop = 0.01,),
              bounds = (pos_prop = (1e-6, 0.5),))
```

Names are those of `dfe_parameters`, plus `syno_orientation_error` and the kernel's
parameter; an unknown name is an error. A free parameter that ends on a bound is
reported with a warning: it is not an estimate. GammaExpo's `pos_prop` upper bound of
0.1 is GRAPES'.

## Reusing tables

The expected-SFS tables depend on the sample size and the integrator, grid, dominance
and expectation, not on the data. `dfe_ml` builds them in a fresh `DFECache()` unless one
is passed:

```julia
cache = MKtest.DFECache()
fits = [MKtest.dfe_ml(s, d, MKtest.GammaZero(); cache = cache) for (s, d) in zip(sfs_list, div_list)]
```

A `DFECache` may be shared across threads. There is no module-level cache.

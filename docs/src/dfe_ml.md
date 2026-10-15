# DFE and alpha by maximum likelihood (`dfe_ml`)

`MKtest.dfe_ml` fits a distribution of fitness effects (DFE) to the SFS of synonymous
and nonsynonymous sites and estimates alpha from the substitutions it predicts. It is
the GRAPES model (Galtier 2016) in pure Julia:

- the nonsynonymous SFS is the DFE-weighted Poisson random field sojourn of mutations
  with selection coefficient `S = 4 Ne s`;
- demography and other distortions of the SFS are absorbed by one nuisance factor `r_j`
  per frequency bin, shared by both classes, and profiled out with theta;
- ancestral misorientation is a free rate, `misorient`;
- alpha is the share of observed Dn not explained by mutations with `S <= nearly_neutral`
  (default 5).

```julia
using MKtest
sfs, div, _ = MKtest.parse_sfs("data/MKdata_tgp_may2023.txt"; sample_size = 661)
s20 = MKtest.project(sfs, 20)                 # 20 individuals = 40 chromosomes

r = MKtest.dfe_ml(s20, div, MKtest.GammaZero())
```

```
DFEResult: 1 fit(s), 1 dataset(s), models GammaZero, 1 config(s)
  .inference  1 x 18
  .params     3 x 10
  .sfs        39 x 9
  .divergence 1 x 12
  .config     1 x 21

 Row | fit  dataset  model      converged  loglik    alpha    alpha_ci
   1 |   1        1  GammaZero       true  -321.377  0.11446  (0.103515, 0.125515)
```

`sfs` must hold every bin `1..n-1` (no frequency cutoff). Pass a `Vector` of SFS and of
divergence vectors to fit several datasets at once, with `labels` to name them.

## The result: `DFEResult`

`dfe_ml` returns a `DFEResult`: a set of tables (`DataFrame`s) with the same columns
whatever the model. Every row carries `fit` (one id per fitted model x dataset x
settings), `dataset` and `model`, so results of different models or calls stack with
`vcat` and join on `fit`.

| table | one row per | columns |
|---|---|---|
| `inference` | fit | `model`, `config`, `n` (chromosomes), `converged`, `loglik`, `n_params`, `theta`, `alpha`, `alpha_ci`, `omega_a`, `omega_a_ci`, `omega_na`, `omega_na_ci`, `omega_obs`, `alpha_div`, `alpha_dfe` |
| `params` | fitted parameter | `parameter`, `kind` (`dfe` / `nuisance`), `value`, `lower`, `upper`, `free`, `at_bound` |
| `sfs` | SFS bin | `j`, `pn_obs`, `ps_obs`, `pn_exp`, `ps_exp` (the fitted expectation), `r_j` |
| `divergence` | fit | `Dn_obs`, `Ds_obs`, `Dn_total`, `Dn_pred`, `lambda`, `Dn_sdel`, `Dn_wdel`, `Dn_wadv`, `Dn_sadv` |
| `config` | distinct setting | the model axes, `nearly_neutral`, `alpha_ci_threshold`, `fixed`, `bounds`, `mktest_version`, ... |
| `extra` | | engine-specific tables (a `Dict`) |

```julia
r.params
#  parameter  kind      value       lower  upper  free  at_bound
#  del_shape  dfe       0.0757      0.05   100.0  true  false
#  del_mean   dfe       1071.41     1.0    1.0e7  true  false
#  misorient  nuisance  1.38e-15    0.0    0.1    true  true

first(r.sfs, 3)
#  j  pn_obs   ps_obs   pn_exp    ps_exp    r_j
#  1  14838.0  12561.0  14760.0   12639.0   1.0
#  2   5227.0   5094.0   5273.5    5047.5   0.799
#  3   2968.0   3070.0   2994.0    3044.0   0.723

r.divergence
#  Dn_obs   Ds_obs   Dn_total  Dn_pred  lambda   Dn_sdel  Dn_wdel  Dn_wadv  Dn_sadv
#  35903.0  45088.0  35903.0   31793.5  0.00573  37.2     31756.3  0.0      4109.5
```

- **Intervals** are `(lo, hi)` tuples: `r.inference.alpha_ci[1]`. The CI is GRAPES'
  likelihood interval (`alpha_ci_threshold = 1.92`, half the chi-square(1) 95%
  quantile); `alpha_ci_threshold = 0` turns it off (NaN bounds).
- **`at_bound`** marks a free parameter that ended on its bound: its value is the
  bound, not an estimate (a warning is also printed).
- **`divergence`** splits the expected Dn by S class, as GRAPES' `prop_subst_*`:
  strongly deleterious (`S <= -5`), weakly deleterious, weakly advantageous (`0 < S <= 5`)
  and strongly advantageous. `alpha = Dn_sadv / Dn_total`, where `Dn_total` is the
  observed Dn (the default `SaturatedDivergence`).
- **Two readings of alpha.** `alpha_div` uses the observed Dn as denominator (it equals
  `alpha` by default); `alpha_dfe` uses the Dn the fitted DFE predicts (`Dn_pred`).
  They agree when the DFE explains the observed divergence.

Accessors:

```julia
MKtest.params_wide(r)              # one row per fit, one column per parameter
MKtest.param_values(r, :del_mean)  # one value per fit (missing where absent)
MKtest.expected_sfs(r, 1)          # [j pn_obs ps_obs pn_exp ps_exp] of fit 1
MKtest.rj(r, 1)                    # its r_j
both = vcat(r1, r2)                # stack results; fit ids are renumbered
```

Saving: `MKtest.write_dfe_result("out/tgp", r)` writes one CSV per table
(`out/tgp_inference.csv`, `_params.csv`, ...; intervals as `_low` / `_high`) and
`MKtest.read_dfe_result("out/tgp")` reads them back. `JLD2.save("r.jld2", "result", r)`
keeps everything in one file.

## Parameter names

Every parameter is named `<component>_<quantity>`, so the same quantity has the same
name in every model:

- components: `del` (S < 0), `ben` (S > 0), `neu` (a point mass at S = 0), `lethal`, or
  a family prefix (`beta_`, `fgm_`, `spike_`);
- quantities: `prop` (mixture weight), `mean` (mean |S| of the component, a positive
  number), `shape`, `sigma`, `shift`.

| model | parameters | GRAPES / earlier names |
|---|---|---|
| `GammaZero()` | `del_shape`, `del_mean` | `negGshape`, `negGmean` |
| `GammaExpo()` | + `ben_prop`, `ben_mean` | + `pos_prop`, `posGmean` |
| `GammaGamma()` | + `ben_prop`, `ben_shape`, `ben_mean` | + `pos_prop`, `posGshape`, `posGmean` |
| `ExpoZero()`, `ExpoExpo()`, `ExpoGamma()` | `del_mean` (+ the `ben_` set) | `negGmean` (+ the `pos` set) |
| `DisplGamma()` | `del_shape`, `del_mean`, `del_shift` | `negGshape`, `negGmean`, `s0` |
| `LogNormalZero()` | `del_mean`, `del_sigma` (sd of log abs(S)) | `negLNmean`, `negLNsigma` |
| `NeuGamma()` | `neu_prop`, `del_shape`, `del_mean` | `p_neu`, ... |
| `ScaledBeta()` | `beta_a`, `beta_b`, `beta_mass` | `Ba`, `Bb`, `med_prop` |
| `FGMBesselK()` | `fgm_m`, `fgm_nsz2`, `fgm_sigma`, `fgm_scale` | `BKm`, `BKnsz2`, `BKsigma`, `BKscale` |
| `FixedSpikes()` | `spike_w1` ... (stick-breaking fractions; `spike_weights` gives the weights) | `w1` ... |
| `WithLethal(m)` | + `lethal_prop` | + `p_lethal` |
| nuisance | `misorient`, `nb_dispersion` (negative binomial), `ne_ratio` (Route A) | `syno_orientation_error`, `a`, `Ne_ratio` |

`MKtest.from_legacy_name("negGmean")` gives `:del_mean` (for results written with the old
names, before 2026-10-15); `MKtest.grapes_param_name(:ben_prop)` gives `"pos_prop"`
(GRAPES' output columns and `-fixed_param` files).

## Fixing and bounding parameters

```julia
MKtest.dfe_ml(s20, div, MKtest.GammaExpo(); fixed = (ben_prop = 0.01,),
              bounds = (ben_mean = (1e-3, 100.0),))
```

Names are those of the model plus `misorient` and the kernel's parameter; an unknown
name is an error.

**`ben_prop` is bounded at 0.1 by default**, GRAPES' cap (its `MAX_ORIENTATION_ERROR`
reused as a mixture weight), kept so results compare with the GRAPES binary. When the
cap binds (`at_bound`), widen it for a scientific fit:
`bounds = (ben_prop = (1e-6, 0.99),)`. On human data the uncapped beneficial class has a
mean `S` well below 1: an effectively neutral class rather than adaptive mutations,
which is why GammaExpo is best reported next to NeuGamma.

## Model comparison

```julia
cmp = MKtest.dfe_model_comparison(s20, div;
          models = [MKtest.GammaZero(), MKtest.GammaExpo(), MKtest.NeuGamma()],
          bounds = (ben_prop = (1e-6, 0.99),))
cmp.aic
#  model      n_params  loglik    aic      delta_aic  weight  alpha
#  GammaExpo         5  -318.273  646.546  1.860      0.259   0.0591
#  NeuGamma          4  -318.343  644.685  0.0        0.656   0.0704
#  GammaZero         3  -321.377  648.753  4.068      0.086   0.1145
cmp.lrt
#  model_null  model_alt  df  D      p_value
#  GammaZero   GammaExpo   2  6.207  0.0449
#  GammaZero   NeuGamma    1  6.068  0.0138
cmp.result                # the DFEResult of the three fits
```

- `bounds` and `fixed` reach each model restricted to the parameters it has.
- A model is nested in another when its free parameters are a strict subset of the
  other's (GammaZero in GammaExpo and in NeuGamma; ExpoZero in GammaZero).
- Fits are compared only on the same data and the same likelihood function: models fitted
  with different integrators, kernels or projections are in different `group`s.
- Both tables are sorted by log-likelihood; pick the best model by `argmin(aic)`, not row
  1.
- `MKtest.dfe_model_comparison(r)` compares fits already in a `DFEResult`
  (e.g. `vcat(r1, r2)`).

**GammaExpo / GammaGamma alpha is weakly identified.** Small numbers of strong
beneficial mutations and large numbers of weak ones predict nearly the same SFS, so
alpha can move along a ridge at almost constant likelihood. Report it with its
log-likelihood and the comparison against GammaZero, not alone.

## The axes

Each keyword of `dfe_ml` is one independent axis: an abstract type whose interface the
engine calls, so a new subtype of one axis needs no change to another.

| keyword | abstract type | default | alternatives | a new subtype implements |
|---|---|---|---|---|
| model (argument) | `AbstractDFEModel` | | the models above; your own ([tutorial](new_dfe.md)) | `dfe_parameters` + `distribution` (or `dfe_density!`) |
| `expectation` | `AbstractSFSExpectation` | `PRFSojourn()` | | `expectation_counts(e, s_eff, n, integ)` |
| `dominance` | `AbstractDominance` | `Additive()` | `Dominance(h)`, `DominanceFunction(f)` | `effective_s(d, S)` |
| `demography` | `AbstractDemography` | `NuisanceRj()` | `Equilibrium()` (r_j = 1) | `free_rj(d)` |
| `ne_model` | `AbstractNeModel` | `ConstantNe()` | `AncestralNeRatio(r)`, `AncestralNeRatio(:free)` | `ne_ratio(m)` |
| `div_model` | `AbstractDivergenceModel` | `SaturatedDivergence()` | `DFEConstrainedDivergence()`, `NoDivergence()` | `uses_divergence(d)`, `divergence_free(d)` |
| `kernel` | `AbstractLikelihoodKernel` | `PoissonKernel()` | `NegativeBinomialKernel()` | `log_kernel(k, count, mu, logfact, disp)` |
| `integrator` | `AbstractIntegrator` | `ExactSojourn()` | `RombergIntegrator()`, `GaussLegendreSojourn()` | `integrate(i, f, a, b)` |
| `quadrature` | `AbstractDFEQuadrature` | `MidpointDensity()` | `IntervalCDF()` | `fill_dfe_vals!(out, q, m, p, pre)` |
| `s_grid` | `AbstractSGrid` | `GrapesGrid()` | `LogGammaGrid()` | `s_grid_points(g)`, `s_grid_widths(g, grid)` |
| `optimizer` | `AbstractOptimizer` | `LBFGSOptimizer()` | `NelderMeadOptimizer()` | `run_optimizer(o, f, bounds, x0)` |
| `summaries` | `AbstractAlphaSummary` | `(DFEMAlpha(),)` | your own | `summarize(s, ctx::FitContext)` |

Methods added from outside the package need the module prefix
(`MKtest.log_kernel(::MyKernel, k, mu, lf, d) = ...`).

GRAPES' options map to: `-no_div_param` = `div_model = DFEConstrainedDivergence()`,
`-no_div_data` = `NoDivergence()`, `-anc_to_rec_Ne_ratio r` = `ne_model = AncestralNeRatio(r)`,
`-fixed_param` = `fixed`, `-no_syn_orient_error` = `no_syn_orient_error = true`,
`-fold` = `fold = true`, `-nb_rand_start` = `nb_rand_start`.

## Integrators and the GRAPES binary

The default `ExactSojourn()` computes the expected SFS in closed form, to about 1e-16.
GRAPES integrates with Romberg over `[1e-4, 1 - 1e-4]`, which loses mass from the rare
bins (bin 1 is low by 0.19% at 20 chromosomes and 6.4% at 1322). To compare with the
binary, use `integrator = MKtest.RombergIntegrator()`:

```julia
MKtest.dfe_ml(s20, div, MKtest.GammaZero(); integrator = MKtest.RombergIntegrator())
# loglik -321.809, alpha 0.1162 (0.1052, 0.1272); the GRAPES binary: -321.809, 0.1163
```

`MKtest.grapes(sfs, div, "GammaZero"; n = 20)` runs the binary itself (installed through
CondaPkg) and returns its output table, with GRAPES' column names.

On the 1000 Genomes SFS at 10 and 20 individuals, `dfe_ml` under Romberg equals the
binary's log-likelihood at the binary's own parameters for every model tested. The fits
match for GammaZero. For GammaExpo and GammaGamma, `dfe_ml` reaches an equal or higher
likelihood, at a different point of the ridge, so alpha differs.

## Reusing tables

The expected-SFS tables depend on the sample size and on the integrator, grid, dominance
and expectation, not on the data. Pass one cache to reuse them across calls (safe
across threads):

```julia
cache = MKtest.DFECache()
res = vcat((MKtest.dfe_ml(s, d, MKtest.GammaZero(); cache = cache) for (s, d) in zip(sfs_list, div_list))...)
```

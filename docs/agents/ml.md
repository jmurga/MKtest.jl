# dfe_ml engine internals (`src/ml/`)

For agents working on the engine. The user-facing pages are `docs/src/dfe_ml.md` and
`docs/src/new_dfe.md`. Moved from dev (9d752e40) in Task 1, 2026-10-15; the Phase 0
design is `docs/agents/briefs/task1_phase0_report.md` on dev.

## Files

| file | content |
|---|---|
| `interface.jl` | constants, every abstract type with its interface stubs and the small default subtypes |
| `models.jl` | DFE traits (`dfe_parameters`, `dfe_log_scale`, `dfe_start_points`, `dfe_name`), the generic `distribution` path, the 12 shapes, `WithLethal`, `FixedSpikes`, tail masses, CDFs, `dfe_ne_multipliers` |
| `integrators.jl` | Romberg (GRAPES `qromb`), `GaussLegendreSojourn`, `ExactSojourn`, `expected_count` / `dn_correction_term` hooks, numerical tail masses (`_log_quad`, `_mass_below/_mass_above`) |
| `grids.jl` | `GrapesGrid`, `LogGammaGrid`, cell bounds, exterior edges, `fill_dfe_vals!` |
| `cache.jl` | `PRFSojourn` counts, `PrecomputedData`, `DnCorrectionTable`, `DFECache`, `get_precomputed`, `dn_correction`, `discrete_row` |
| `likelihood.jl` | `DataInput`, `ParamLayout`, `LikelihoodObjective` and `_loglik!`, CI, expected divergence |
| `fit.jl` | `build_layout`, starting points, `run_optimizer`, `_fit`, `FitContext`, summaries, Route A, `fit_model` |
| `api.jl` | `dfe_ml`, `dfe_model_comparison`, `_aic_lrt_tables`, internal `_objective_for` / `_loglik_at` |

## Algorithm (GRAPES)

1. S grid (`GrapesGrid`: 4000 points over [-500,-0.5], [-0.5,0.5], [0.5,500], [500,10000]).
2. `expected_counts[i, j]`: expected SFS bin j at S_i, theta = 1, from
   `expectation_counts(expectation, s_eff, n, integrator)`; `s_eff = effective_s(dominance, S)`.
3. Per evaluation: DFE weights `fill_dfe_vals!` x widths, `exp_sel = Ln * sum_i w_i counts[i, :]`
   (+ point masses + exterior terms if the expectation uses them), folded if asked;
   `exp_neut = Ls / j`, mixed with its complement by `syno_orientation_error`.
4. theta from bin 1 of Ps + Pn, `r_j` profiled from the other bins (if `free_rj`).
5. Log-likelihood: `log_kernel` over every bin; divergence saturated (constant) or
   DFE-predicted (`DFEConstrainedDivergence`).
6. Multi-start optimisation: model starts + box centre + `nb_rand_start` (MersenneTwister(42)).
7. At the optimum the same objective is re-evaluated, and `FitContext` carries the fitted
   SFS and the divergence decomposition to the summaries (DFEMAlpha, FWWAlpha).

## Traps

- **Name the integrator when quoting GRAPES parity.** Dev's human_test numbers (lnL
  -328.3279, alpha 0.1797) were under `RombergIntegrator()`; `ExactSojourn` gives -327.568
  / 0.178 on the same data. The difference is the quadrature, not a regression.
- **GRAPES' PseudoNewton stalls at its start point on flat surfaces** (DisplGamma,
  ScaledBeta on human data); a higher Julia lnL there is not a discrepancy. The pos_prop x
  posGmean ridge makes GammaExpo/GammaGamma land on different local optima (Delta lnL 1-2).
- **Ambiguity trap.** Generic fallbacks (`dfe_density!`, masses, CDFs,
  `discrete_components`, `fill_dfe_vals!`) are untyped beyond the model. Typing their
  trailing arguments would make a user's untyped `dfe_density!(out, m::MyType, p, s)`
  ambiguous.
- **Merging NamedTuples drops duplicates.** `dfe_parameters` of a wrapper is built by
  splatting; `WithLethal` therefore refuses a base that already has `p_lethal`.
- **Log scales of the built-ins are GRAPES', not the `lower > 0` default**
  (`negGshape` is linear). A user model gets the default; its optimizer path differs from
  a built-in with the same density.
- **One-ulp changes move ridge fits.** dev measured GammaExpo's alpha moving 0.8% under a
  1-ulp perturbation of the neutral term. The new engine sums the SFS in a plain loop where
  dev used BLAS `mul!`. `test/ml/dev_reference.jl` asserts the objective to 1e-10
  relative and does not assert ridge-model parameters.
- **Ne ratio**: the divergence side uses the DFE rescaled by `dfe_ne_multipliers` (S-scale
  parameters x r, `med_prop` / r as GRAPES). Unknown parameter names error rather than
  being skipped (dev skipped `negLNmean`, `BKscale` and FixedSpikes silently).
  Route A (`AncestralNeRatio(:free)`) needs `SaturatedDivergence`.
- **Optim 2.** main's LsqFit 0.16 needs NLSolversBase 8, so the branch runs Optim 2.3.2
  while dev ran 1.13.3.

## Divergence roles

| `div_model` | divergence in the LL | Dn | GRAPES |
|---|---|---|---|
| `SaturatedDivergence()` | yes, constant | free | default |
| `DFEConstrainedDivergence()` | yes | DFE-predicted | `-no_div_param` |
| `NoDivergence()` | no | -- | `-no_div_data` |

`r_j` corrects the SFS shape (recent demography); the Ne ratio corrects divergence
scaling (ancestral Ne). They are orthogonal and can both be free.

## Not moved (stay on dev)

pde/ and wf/ expectations, fitted demography, `dfe_bayes`, multi-population
`shared_params`, `TwoLineageNe` / `ThreeLineageDivergence` / `MultiLineageDivergence`,
`remove_poly`, String model names and `fixed_param` strings, the parallel-over-datasets
branch. Their interfaces (`AbstractSFSExpectation`, `default_demography`, `uses_exterior`)
are the plug-in points.

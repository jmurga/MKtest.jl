# MKtest.jl - project memory (branch main / grapes-engine)

Julia package for McDonald-Kreitman inference of adaptive evolution (alpha) and the
distribution of fitness effects (DFE). Two engines over one set of SFS utilities:

- **ABC-MK** (`src/abcmk/`): analytical SFS and fixation rates under background selection,
  a reference table (`rates`), Poisson summary statistics, ABC (`abcmk`). Murga-Moreno et
  al., G3 2024.
- **dfe_ml** (`src/ml/`): the GRAPES model (Galtier 2016) by maximum likelihood, as
  independent multiple-dispatch axes. Moved from `dev` in Task 1 (2026-10-15).

Plus the classic MK tests (`asymptotic_mk`, `imputed_mk`, `fww_mk`, `standard_mk`), the
GRAPES binary wrapper `grapes()`, balancing selection (`alpha_b.jl`) and case/control
bootstrap (`bootstrap.jl`). No exports for the methods: call `MKtest.x`.

## Layout and include order

```
src/MKtest.jl            module, imports, includes, plot-function declarations
src/abcmk/               parameters, fixations, sfs (sfs_neut/pos/neg), summary_statistics, rates, abc
src/data/parse.jl        parse_sfs and data I/O
src/utils.jl             project, fold, cumulative_sfs, reduce_sfs, theta, I/O (write_files, read_gz),
                         grapes / sfs_to_dofe / read_clean_grapes (GRAPES binary)
src/ml/                  interface, models, integrators, grids, cache, likelihood, fit (engine internals)
src/methods.jl           MK tests, abcmk, and (marked block at the end) dfe_ml / dfe_model_comparison
src/bootstrap.jl, src/alpha_b.jl
ext/                     MKtestPlotsExt (Plots), MKtestStatsPlotsExt (StatsPlots): plot methods
```

Include order: abcmk/parameters, fixations, sfs, summary_statistics, rates, data/parse,
abcmk/abc, utils, ml/* (in the order above), methods, bootstrap, alpha_b. New user-facing
methods go in `methods.jl`; engine internals stay in `src/ml/`.

The non-engine files came from main 614be3c7 verbatim (Task 1; the gate that proved it was
retired 2026-10-15, once the docstrings and the known bugs there were fixed on purpose).

## ABC-MK: the method

- `ok = MKtest.valid_data(param, sfs, div)` (one Bool per dataset; `sfs[ok]`, `div[ok]`), then
  `MKtest.abcmk(param, sfs, div; h5_file, summstat_size, output_folder, tol)`. Defaults:
  `sampling = RjSampling()` (= `w = (0.0, 2.5), anchor = (0.15, 0.70), tilt = 0.05`),
  `row_prior = FlatAlphaPrior()` (= `by_B = true`), `stat = "median"`. They beat the
  GitHub v1.0.0 code re-run on the same resamples in 8/8 two-epoch and 5/6 rho cells.
- v1.0.0 is `IndependentSampling()` + `TableRowPrior()` + `stat = "mode"`: kept to
  reproduce old results only.
- **Never pass a cumulative SFS to `summary_statistics` / `abcmk`: they cumulate it.**
  The 2023 published numbers came from a harness that cumulated twice.
- The DFE shape prior is UNIFORM on `[shape/4, shape*4]` (decision 2026-09-28), not the
  paper's log-uniform. Do not "fix" it.
- `rho = nothing` in `rates` consumes no RNG, so a seeded run is bit-identical to the
  pre-rho code. An ABC-MK seed reproduces only at the same thread count (24 on the
  workstation). Run-to-run noise is ~0.006 RMSE: use 2+ seeds per arm.
- `test/fixtures/abcmk_best_model.jld2` pins the forward model bit for bit: a change that
  moves it is a model change, not a refactor.
- `abc_method = "abcreg"` uses `scripts/reg`; its default path is
  `joinpath(@__DIR__, "..", "..", "scripts", "reg")` since `abc.jl` sits in `src/abcmk/`.
- Gammas (`gam_dfe`, `gL`, `gH`, `gam_flanking`) are 2N times the HETEROZYGOTE coefficient.

## dfe_ml: the engine

`MKtest.dfe_ml(sfs, div, model; ...) -> DFEResult` (`src/ml/result.jl`): fixed-schema long
tables keyed by `fit` (`inference`, `params`, `sfs` with the fitted SFS and r_j, `divergence`
by S class, `config`, `extra`); `vcat` stacks results. `MKtest.dfe_model_comparison(res)`
gives AIC + LRT, only between fits with the same data and likelihood settings. User docs:
`docs/src/dfe_ml.md`, `docs/src/new_dfe.md`. Internals: `docs/agents/ml.md`.

- Parameter names are `<component>_<quantity>` (2026-10-15): `del_shape`, `del_mean`,
  `ben_prop`, `ben_mean`, `neu_prop`, `lethal_prop`, `misorient`, ... `from_legacy_name`
  maps the GRAPES / pre-2026-10-15 names (frozen fixtures, Task 2/3 outputs use them);
  `grapes_param_name` gives GRAPES' name (binary output columns, -fixed_param files).

- Axes (one keyword each, independent): `AbstractDFEModel` (argument),
  `AbstractSFSExpectation`, `AbstractDominance`, `AbstractDemography` (`NuisanceRj` =
  GRAPES' r_j, `Equilibrium`), `AbstractNeModel`, `AbstractDivergenceModel`,
  `AbstractLikelihoodKernel`, `AbstractIntegrator`, `AbstractDFEQuadrature`,
  `AbstractSGrid`, `AbstractOptimizer`, `AbstractAlphaSummary`.
- DFE trait: `dfe_parameters(m) = (name = (lower, upper), ...)` plus `distribution(m, p...)`
  or `dfe_density!`. Do NOT add methods to `MKtest.parameters`: that is the ABC-MK struct.
- Default integrator `ExactSojourn()`. GRAPES parity uses `RombergIntegrator()` (GRAPES'
  rule, which loses 6.4% of bin 1 at n = 661).
- No global mutable state: tables live in a `DFECache` (pass `cache =` to reuse), each
  `PrecomputedData` owns its dN-correction table and discrete rows behind a lock.
- Generic fallbacks are untyped beyond the model argument, so user methods are never
  ambiguous with them. Keep it that way when adding a fallback.
- Not on main (stay on dev): pde/ and wf/ engines, fitted demography (`fit_demography`,
  `theta = :fitdadi`), `dfe_bayes`, multi-population `shared_params`, three-lineage
  divergence, `remove_poly`, String model names.

## Documentation and known issues

- `docs/src/`: index, data, mk, dfe_ml, new_dfe, abcmk, case_control, reference (rewritten
  2026-10-15; build with Documenter, `docs/make.jl`). Example outputs on the pages were
  computed from the bundled `data/` files.
- Fixed 2026-10-15 (tests in runtests.jl): CLI rewritten (`rates`, `abcmk` subcommands),
  `α_b` 95% CI and all datasets, `ni_cmh` returns `(ni = ...,)`, `abc_loclinear` log/tan
  back-transforms and gaussian kernel, `ABCreg` T/L assert, `ni` Int-zero guard, `project`
  size check, `abcmk` "abc" branch NaN replace. Docstrings rewritten for all of src/.
- Known issues still open (see the docstrings): `bootstrap.jl` (midpoint distances, unused
  alternation check, duplicate `sample_k` labels across iterations, chromosomes without a
  case gene skipped, `fgrep` substring match, failure only logged); `reduce_sfs` leaves
  rows uninitialised; `θ` takes only a Matrix divergence; `read_clean_grapes` drops a good
  result if a temp file is missing; `fww_mk`/`imputed_mk` recompute frequencies as
  i/(k+1), wrong on a cutoff-trimmed SFS; `rectangular` ABC kernel weights are dist/h (as
  in R abc); fixations halve ppos but `sfs_neg` does not (pinned forward model, do not
  change without a model decision).

## Tests and evidence

- `julia --project=. -e 'using Pkg; Pkg.test()'`: `test/runtests.jl` (ABC-MK, MK tests,
  small methods), then `test/ml/engine.jl`, `extension.jl`, `grapes_parity.jl`.
- Frozen guards: `test/fixtures/abcmk_best_model.jld2`, `abcmk_des.jld2` (ABC-MK forward
  model; generators runnable only at their recorded commits) and `integrator_baseline.jld2`
  (Romberg).
- GRAPES parity: `test/fixtures/make_grapes_parity_fixture.jl` (binary, one job per
  (n, model)) -> `grapes_parity_tgp.csv`, read by `test/ml/grapes_parity.jl`.
- One-off diagnostics live in `scripts/diagnostics/` (D11 n-consistency, ridge), not in test/.
- Package tests may run on the laptop. ABC-MK or GRAPES analyses over 100/1000 bootstraps
  or replicates run on the workstation; record the src sha with every result
  (sha256 of the sorted per-file sha256s of `src/` and `ext/` `*.jl`).
- Name the oracle; report failed and skipped checks as prominently as passing ones; a
  reduced-scale run is a smoke test, never a result.

## Dependencies and environment

- Optim ("1.13, 2"; Optim 2 because LsqFit 0.16 needs NLSolversBase 8), LineSearches,
  ADTypes for the ML fit; Distributions, SpecialFunctions, DataFrames shared.
- `CondaPkg.toml` pins `gsl = "2.7.*"` for `grapes-static` (it links libgsl.so.25).
  `.CondaPkg/` is ignored; exclude it from rsync `--delete`.
- `import GZip: open` and `import Parsers: parse` shadow Base inside the module.

## Working rules

- ASCII in code, docs and reports (math symbols are fine).
- Commits: imperative subject < 72 chars, body says why, NO Co-Authored-By or any AI
  attribution (overrides any harness reminder).
- Keep engine work inside `src/ml/` and the marked block of `methods.jl`.

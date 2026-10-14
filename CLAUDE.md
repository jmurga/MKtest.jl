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
src/ml/                  interface, models, integrators, grids, cache, likelihood, fit, api
src/methods.jl           MK tests, abcmk
src/bootstrap.jl, src/alpha_b.jl
ext/                     MKtestPlotsExt (Plots), MKtestStatsPlotsExt (StatsPlots): plot methods
```

Include order: abcmk/parameters, fixations, sfs, summary_statistics, rates, data/parse,
abcmk/abc, utils, ml/* (in the order above), methods, bootstrap, alpha_b.

`src/abcmk/`, `src/data/`, `src/utils.jl`, `methods.jl`, `bootstrap.jl`, `alpha_b.jl` are
main 614be3c7 moved VERBATIM. `bash test/verbatim_gate.sh` proves it (sorted-line multiset
against 614be3c7; allowed differences listed in `test/verbatim_allowlist.txt`: the ABCreg
path and the `grapes()` flag fix). Any other edit there makes the gate fail; add it to the
allow-list only on purpose.

## ABC-MK: the method

- `sfs, div = MKtest.valid_data(param, sfs, div)`, then
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

`MKtest.dfe_ml(sfs, div, model; ...) -> DataFrame` (one row per dataset) and
`MKtest.dfe_model_comparison(...)` (AIC + LRT DataFrames). User docs: `docs/src/dfe_ml.md`,
`docs/src/new_dfe.md`. Internals: `docs/agents/ml.md`.

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

## Tests and evidence

- `julia --project=. -e 'using Pkg; Pkg.test()'`: main's testsets, then `test/ml/engine.jl`
  and `test/ml/extension.jl` (and the parity/dev-reference tests once their fixtures exist).
- GRAPES parity: `test/fixtures/make_grapes_parity_fixture.jl` (binary, one job per
  (n, model)) -> `grapes_parity_tgp.csv`, read by `test/ml/grapes_parity.jl`.
- Refactor guard: `test/fixtures/make_dev_reference.jl` (run in a dev tree at 9d752e40).
- D11 n-consistency: `test/consistency_n.jl` (standalone; band fixed in its header).
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
- Do not touch the ABC-MK, MK-test, parse or bootstrap code when working on the engine.

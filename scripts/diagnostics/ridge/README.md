# Ridge diagnostics: why GammaExpo / GammaGamma alpha differs from the GRAPES binary

Run on the laptop, 2026-10-15, src sha 622e6e0b, Julia 1.13.1, fastDFE 1.4.1. Each
script names its inputs; none needs the GRAPES binary (they read grapes_parity_tgp.csv).

| script | question | run |
|---|---|---|
| `alpha_at_grapes.jl` | is alpha defined as in GRAPES? Julia's alpha and lnL at GRAPES' fitted parameters, against Julia's optimum | `julia --project=. scripts/diagnostics/ridge/alpha_at_grapes.jl` |
| `drift_synthetic.jl` | is alpha's rise with n a property of the method? GammaZero data from the engine's own model, projected | `JULIA_NUM_THREADS=4 julia --project=. scripts/diagnostics/ridge/drift_synthetic.jl` |
| `fastdfe_export.jl` + `fit_fastdfe.py` | where does an independent implementation (fastDFE GammaExp) put the optimum? | `julia --project=. .../fastdfe_export.jl <dir>`; `<fastdfe venv>/bin/python -I .../fit_fastdfe.py <dir>` |

`params.tsv` and `fastdfe_results.json` are the recorded inputs and output of the fastDFE
comparison. Results and interpretation: `docs/agents/briefs/task1_phase1_handoff.md` on dev.

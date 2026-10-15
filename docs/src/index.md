# MKtest.jl

MKtest.jl estimates the rate of adaptive evolution from polymorphism and divergence:
alpha, the proportion of nonsynonymous substitutions fixed by positive selection, and
the distribution of fitness effects (DFE) of new mutations. It brings three families of
methods behind one data format:

| family | what it is | entry points |
|---|---|---|
| **MK tests** | the McDonald-Kreitman test and its extensions, computed from the counts alone | `standard_mk`, `fww_mk`, `imputed_mk`, `asymptotic_mk`, `dos`, `ni` |
| **dfe_ml** | the GRAPES model (Galtier 2016) by maximum likelihood: a DFE fitted to the site frequency spectrum (SFS), alpha from the divergence it predicts | `dfe_ml`, `dfe_model_comparison` |
| **ABC-MK** | approximate Bayesian computation over an analytical model with background selection and weakly and strongly beneficial mutations (Uricchio et al. 2019; Murga-Moreno et al. 2024) | `rates`, `abcmk` |

and the tools around them: parsing gene-level data into an SFS, projecting and folding
it, and testing a gene set against matched control sets.

## Installation

MKtest.jl needs Julia 1.10 or later.

```julia
using Pkg
Pkg.add(url = "https://github.com/jmurga/MKtest.jl")
```

Plotting is optional and loads with Plots.jl or StatsPlots.jl
(`using Plots` / `using StatsPlots`; see [MK tests](mk.md) and
[Case and control sets](case_control.md)). The GRAPES binary used by `MKtest.grapes` is
installed on first use through CondaPkg.

Most methods run on several datasets at once and use threads: start Julia with
`julia -t auto` (or set `JULIA_NUM_THREADS`).

## Conventions

- **Nothing is exported.** Call every function as `MKtest.f`.
- **Sample sizes are in individuals.** `parse_sfs(...; sample_size = 661)` reads 661
  diploids (1322 chromosomes, `ploidy = 2`); `project(sfs, 20)` projects to 20
  individuals (40 chromosomes, an SFS of 39 bins). Inbred or haploid lines use
  `ploidy = 1`.
- **The SFS is raw, not cumulative**: a matrix `[DAC Pn Ps]` with one row per derived
  allele count. Methods that need a cumulative spectrum compute it themselves.
- **Divergence is `[Dn, Ds, Ln, Ls]`**: nonsynonymous and synonymous substitutions, then
  the number of nonsynonymous and synonymous sites.
- **Selection coefficients are population scaled**: `S = 4 Ne s` in dfe_ml, `gamma = 2 N s`
  in the ABC-MK parameters.

## Quick start

```julia
using MKtest

# 1000 Genomes, 661 African individuals: one SFS for all genes
sfs, div, _ = MKtest.parse_sfs("data/MKdata_tgp_may2023.txt"; sample_size = 661)

MKtest.standard_mk(sfs, div)                          # the MK test
MKtest.asymptotic_mk(sfs, div; cutoff = [0.0, 0.7])   # asymptotic alpha

# the DFE and alpha by maximum likelihood, on the SFS projected to 20 individuals
r = MKtest.dfe_ml(MKtest.project(sfs, 20), div, MKtest.GammaZero())
r.inference          # alpha, its CI, omega_a, the log-likelihood
r.params             # the fitted DFE (del_shape, del_mean) and nuisance parameters
```

The data files used throughout these pages are in the package's `data/` folder
(`pkgdir(MKtest)/data`); [Input data](data.md) describes them.

## Pages

- [Input data](data.md): the file format, the bundled datasets, `parse_sfs`, projection.
- [MK tests](mk.md): standard, FWW, imputed and asymptotic MK, DoS, NI, alpha(x).
- [DFE and alpha by maximum likelihood](dfe_ml.md): `dfe_ml`, its result, model comparison.
- [How to add a new DFE](new_dfe.md): extending `dfe_ml` with your own distribution.
- [ABC-MK](abcmk.md): the analytical model, the reference table, the inference.
- [Case and control sets](case_control.md): a gene set against matched controls.
- [Reference](reference.md): every documented function.

## Citation

- ABC-MK: Murga-Moreno J, Casillas S, Barbadilla A, Uricchio L, Enard D. *An efficient
  and robust ABC approach to infer the rate and strength of adaptation.* G3 (2024).
- The analytical model: Uricchio LH, Petrov DA, Enard D. *Exploiting selection at linked
  sites to infer the rate and strength of adaptation.* Nat Ecol Evol (2019).
- dfe_ml follows GRAPES: Galtier N. *Adaptive protein evolution in animals and the
  effective population size hypothesis.* PLoS Genet (2016).

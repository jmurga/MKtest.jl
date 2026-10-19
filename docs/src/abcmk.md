# ABC-MK

ABC-MK estimates alpha, its weak and strong components (alpha_W, alpha_S) and the DFE
of deleterious mutations by approximate Bayesian computation (Murga-Moreno et al. 2024).
It compares the observed alpha(x) curve and divergence with those of an analytical model
(Uricchio et al. 2019) that includes:

- a gamma DFE of deleterious mutations (shape, mean `gam_dfe`);
- weakly (`gL`) and strongly (`gH`) beneficial mutations, the part of alpha due to the
  weak ones (`al_low`, alpha_W), and the total alpha (`al_tot`);
- background selection `B` from deleterious flanking sites (`gam_flanking`), with
  recombination `ρ`, and its reduction of the fixation probability of beneficial
  mutations (selective interference).

Three steps: a reference table of model predictions (`rates`, once per population and
sample size), summary statistics of the data sampled against it, and the ABC itself
(`abcmk` runs both).

## 1. The model parameters

```julia
using MKtest
param = MKtest.parameters(N = 1000, n = 661, dac = [2, 4, 5, 10, 20, 50, 200, 661, 925])
```

`MKtest.parameters` holds the model. The fields you set:

| field | default | meaning |
|---|---|---|
| `N` | 1000 | population size of the model |
| `n` | 500 | sample size (individuals; `nn = 2n` chromosomes) |
| `dac` | `[2,4,5,10,20,50,200,500,700]` | derived allele counts at which alpha(x) is summarised: ascending, `< 2n`, no singletons |
| `shape`, `gam_dfe` | 0.184, -457 | DFE of deleterious mutations (shape, mean `2 N s`) |
| `gL`, `gH` | 10, 500 | weakly and strongly beneficial `2 N s` |
| `al_low`, `al_tot` | 0.2, 0.4 | weak and total alpha |
| `B`, `ρ`, `gam_flanking`, `Lf` | 0.999, 0.001, -500, 2e5 | background selection, recombination, flanking selection and length |
| `cutoff` | `[0.0, 1.0]` | frequency range of the data |

Most of these are drawn from priors when the table is built; set `N`, `n`, `dac` and
`cutoff` to match the data.

## 2. The reference table: `rates`

```julia
MKtest.rates(param; iterations = 100_000, output_file = "rates_tgp.jld2")
```

Each iteration draws a model from the priors and solves its expected SFS (at the
`dac`) and fixation rates for 37 values of `B`. Priors (keywords, ranges):

| keyword | default | prior |
|---|---|---|
| `alpha` | `[0.1, 0.9]` | `al_tot ~ U` |
| `gam_dfe` | `[-1000, -200]` | mean deleterious `2 N s` |
| `gL`, `gH` | `[1, 10]`, `[100, 1000]` | weakly / strongly beneficial `2 N s` |
| `gam_flanking` | `[-1000, -500]` | flanking `2 N s` |
| `rho` | `nothing` (= `param.ρ`) | a value, or `[lo, hi]` for a uniform prior |
| `B_bins` | 0.1 to 0.975 by 0.025, and 0.999 | background selection values |

The DFE shape is drawn uniformly in `[shape/4, 4 shape]`, and the weak share of alpha
from `0:0.05:1` (`weak = false` for none). The table is written to `output_file`
(JLD2), under the group `"N/n/cutoff"`, with the drawn models, the expected
synonymous and nonsynonymous tails at each DAC, the fixation rates and the options used.

This is the expensive step: build it once per population, sample size and cutoff, on
many threads (`julia -t auto`), and reuse it for every gene set.

## 3. The data

ABC-MK takes the raw SFS (never a cumulative one: it computes the alpha(x) tails
itself) and the divergence of one or several datasets:

```julia
sfs, div, _ = MKtest.parse_sfs("data/MKdata_tgp_may2023.txt"; sample_size = 661,
                               gene_list = "data/DNAonly_nonVIPs.txt.gz")   # 1000 sets

ok = MKtest.valid_data(param, sfs, div)     # one Bool per dataset
sfs, div = sfs[ok], div[ok]
```

`valid_data` drops the datasets the sampler cannot use: those with `Dn` or `Ds` equal
to 0, or with no polymorphism above some DAC of `param.dac`.

## 4. The inference: `abcmk`

```julia
res = MKtest.abcmk(param, sfs, div; h5_file = "rates_tgp.jld2", summstat_size = 1e5,
                   output_folder = "abc_vips", tol = 1000)
res[:inference]     # one row per dataset
```

`abcmk`:

1. samples `summstat_size` rows of the table and turns each into the summary
   statistics the dataset would show under that model (`summary_statistics`), written to
   `output_folder`;
2. runs the ABC with the `tol` closest rows (`abc_method = "abcreg"`, the ABCreg binary
   shipped in `scripts/reg`, or `"abc"`, a pure Julia rejection and local-linear ABC);
3. summarises each posterior (`summary_abc`).

It returns an `OrderedDict`:

- `:inference`: one row per dataset, with the posterior statistic of `α_weak`,
  `α_strong`, `α`, `ωₐ_weak`, `ωₐ_strong`, `ωₐ`, `ωₙₐ`, the DFE (`γ₋`, `β`), `γ₊`,
  `γ₊₊` and `B` (x1000);
- `:quantiles`: the same values as strings with their 95% interval;
- `:stats`: every statistic (minimum, quantiles, median, mean, mode, maximum).

### Settings

The defaults are the measured best configuration:

| keyword | default | alternative |
|---|---|---|
| `sampling` | `RjSampling()`: the data's synonymous SFS distortion is carried onto the model's spectra | `IndependentSampling()`: independent Poisson draws (v1.0.0) |
| `row_prior` | `FlatAlphaPrior()`: table rows drawn flat in alpha within each `B` | `TableRowPrior()`: rows as drawn (v1.0.0) |
| `stat` | `"median"` | `"mode"` (v1.0.0), `"mean"` |
| `alpha`, `B_bins` | the whole table | restrict the prior to an alpha range or `B` values |

`sampling = MKtest.IndependentSampling(), row_prior = MKtest.TableRowPrior(), stat = "mode"`
reproduces ABC-MK v1.0.0.

### Without files: `in_memory = true`

By default `abcmk` writes one `summstat_i.txt.gz` per dataset (about 5 MB at
`summstat_size = 1e5`) and reads it back for the ABC. With `in_memory = true` each
dataset's summaries are built, used and dropped in memory, one dataset per thread, and no
`output_folder` is needed:

```julia
res = MKtest.abcmk(param, sfs, div; h5_file = "rates_tgp.jld2", summstat_size = 1e5,
                   tol = 1000, in_memory = true)
```

The summaries are the same random draws as on disk at the same seed and thread count, so
`abc_method = "abc"` gives exactly the file route's posteriors, and `abc_method = "abcreg"`
runs `abcreg_loclinear`, a Julia port of the ABCreg algorithm (mean/sd scaling, GSL
quantile radius, Epanechnikov weights, GSL's balanced weighted SVD, no heteroscedastic
correction). It equals the binary to its printed precision (6 significant digits) with the
same accepted draws; cite [ABCreg](https://doi.org/10.1186/1471-2156-10-35) as for the
binary.

This is the route for many datasets, such as the controls of a case/control analysis. For
1000 human control sets (n = 661, 8 DACs, a 3.7M-row table, `summstat_size = 1e5`,
`tol = 1000`, 8 threads on a 4-core laptop, one run of the same code on the development
branch) the in-memory ABCreg took 71 s against 685 s for `summary_statistics` plus the
binary (470 s writing the files, 215 s in ABCreg), with the same posteriors to 6 digits
for all 1000 sets.

### Many calls on one table: `load_rates`

Each `abcmk` call reads the rates table and computes the row-prior weights again. When the
same table serves several calls in one Julia session (cases and controls, several gene
sets), load it once and pass it as `rates`:

```julia
t = MKtest.load_rates("rates_tgp.jld2", param)
cases    = MKtest.abcmk(param, sfs_cases, div_cases; rates = t, summstat_size = 1e5,
                        tol = 1000, in_memory = true)
controls = MKtest.abcmk(param, sfs_ctl, div_ctl; rates = t, summstat_size = 1e5,
                        tol = 1000, in_memory = true)
```

The results are the same as with `h5_file`, which keeps working; `rates` also accepts a
path. Give one of the two. The table must match `param` (`N`, `n`, `cutoff` and every DAC
of `param.dac`), and stays in memory, with its cached weights, as long as `t` does.
`summary_statistics` takes `rates` too.

On a 3.7M-row human table with 1000 control sets (8 threads, laptop, five warm calls each,
measured on the development branch) the setup before the Poisson draws took a median
1.20 s per call before this change, 0.76 s through `h5_file` now and 0.24 s with
`rates = t`, plus 7.7 s once for `load_rates`, most of it compiling the table's types on
first use.

### Things to know

- **The SFS must be raw.** Passing a cumulative SFS cumulates it twice.
- **Run-to-run noise.** Two runs with different seeds differ by about 0.006 in alpha
  RMSE; compare configurations over two or more seeds. A seed reproduces a run only at
  the same thread count.
- **Empirical data** are usually summarised below 70% frequency
  (`cutoff = [0.0, 0.7]` in the parameters and the table) to limit the effect of
  ancestral misorientation; simulations need no cutoff.
- `summary_statistics` clears earlier `summstat_*`, `alphas_*` and `out_*` files in
  `output_folder`; use one folder per analysis.

The model functions behind the table (`analytical_alpha`, `fix_pos_sim`, `sfs_pos`, ...)
are listed in the [Reference](reference.md).

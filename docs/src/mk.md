# MK tests

The McDonald-Kreitman (MK) family compares polymorphism and divergence at nonsynonymous
(N) and synonymous (S) sites. Under neutrality `Dn/Ds = Pn/Ps`; adaptive substitutions
inflate `Dn`. All functions take a raw SFS and `[Dn, Ds, Ln, Ls]` (see
[Input data](data.md)), and also accept a `Vector` of each (one result row per dataset,
returned as a `DataFrame`).

The examples use the whole 1000 Genomes table:

```julia
sfs, div, alpha_x = MKtest.parse_sfs("data/MKdata_tgp_may2023.txt"; sample_size = 661)
```

## Standard MK test

`alpha = 1 - (Pn / Ps) (Ds / Dn)` (Smith & Eyre-Walker 2002), with a Fisher exact test
of the 2x2 table.

```julia
MKtest.standard_mk(sfs, div)
# (alpha = -0.598, p_value = 0.0, omega = 0.284, omega_a = -0.170, omega_d = 0.453)
```

`omega = (Dn/Ln) / (Ds/Ls)`, split into its adaptive (`omega_a = alpha * omega`) and
non-adaptive (`omega_d`) parts. Slightly deleterious polymorphisms make the standard
alpha negative, as here.

## Fay, Wyckoff and Wu (FWW)

Polymorphisms below a frequency `cutoff` are removed before computing alpha (Fay et al.
2002), which discards most slightly deleterious variants; 15% is the usual choice
(Charlesworth & Eyre-Walker 2008).

```julia
MKtest.fww_mk(sfs, div; cutoff = 0.15)
# (alpha = 0.008, p_value = 0.580, omega = 0.284, omega_a = 0.0023, omega_d = 0.281)
```

## Imputed MK test

Instead of removing low-frequency polymorphisms, the imputed MK test (Murga-Moreno et al.
2022) subtracts the slightly deleterious ones it estimates from the shape of the SFS
below `cutoff`, which keeps more information for gene-level tests.

```julia
MKtest.imputed_mk(sfs, div; cutoff = 0.15)
# (alpha = 0.008, p_value = 0.331, b = 0.172, f = 0.281, d = 0.547, omega = 0.284, ...)
```

`b`, `f` and `d` are the fractions of new nonsynonymous mutations estimated to be
slightly deleterious, neutral and strongly deleterious.

## Asymptotic MK test

alpha is computed as a function of frequency, `alpha(x)` from the polymorphisms at
frequency `>= x`, and the curve is extrapolated to `x = 1` (Messer & Petrov 2013),
which removes the bias of slightly deleterious polymorphisms. `cutoff` restricts the
frequencies used in the fit; `[0.0, 0.7]` is common for empirical data.

```julia
amk = MKtest.asymptotic_mk(sfs, div; cutoff = [0.0, 0.7])
# amk.asymptotic_mk = 0.0118, CI [0.0081, 0.0154]; fit a + b exp(-c x) with
# a = 0.0118, b = -0.318, c = 21.7

MKtest.asymptotic_mk(sfs, div; cutoff = [0.0, 0.7], model = MKtest.LinearFit())
# asymptotic_mk = 0.0823, CI [0.0679, 0.0966]
```

The fitted model is `ExponentialFit()` (`a + b exp(-c x)`, the default) or `LinearFit()`
(`a + b x`). Only the chosen model is fitted: if the fit fails the result is NaN, with
no silent fallback to the other. The CI is the 95% interval of the asymptote. The
result also holds the `f` and `alpha_x` points used.

`asymptotic_mk` also takes the `alpha_x` table from `parse_sfs`, or a bare alpha(x)
vector: `MKtest.asymptotic_mk(alpha_x; cutoff = [0.0, 0.7])`.

`MKtest.α_x(sfs, div)` returns the alpha(x) table `[DAC frequency alpha(x)]` itself.

## Direction of selection and neutrality index

```julia
MKtest.dos(sfs, div)   # (dos = -0.117, p_value = 0.0)   Stoletzki & Eyre-Walker 2011
MKtest.ni(sfs, div)    # (ni = 1.598, p_value = 0.0)     Rand & Kann 1996
```

`MKtest.ni_cmh(sfs_vector, div_vector)` combines several datasets into one weighted NI
(Stoletzki & Eyre-Walker 2011); it returns the value alone, without a test.

## Many datasets and multiple testing

With a `Vector` of SFS and divergence (for instance from a gene list with one gene per
line and `mode = :gene`), every test returns a `DataFrame`, and `correction = :fdr` or
`:bonferroni` adds a `pvalue_adjusted` column:

```julia
gs, gd, _ = MKtest.parse_sfs("data/MKdata_tgp_may2023.txt"; sample_size = 661,
                             gene_list = "data/ensembl_list.txt", mode = :gene)
MKtest.standard_mk(gs, gd; correction = :fdr)
```

## Plots

With Plots.jl loaded:

```julia
using Plots
MKtest.plot_asymptotic_mk(amk)                    # alpha(x), the fit and its asymptote
MKtest.plot_α_x(alpha_x; cutoff = [0.0, 0.7])     # alpha(x) curves
MKtest.plot_α_x(sfs, div; sample_size = 50)       # from the SFS, projected first
```

`plot_α_x` takes several curves at once and `control = ...` for a band of control
curves; `output_path = "fig.pdf"` saves the figure (png, svg, pdf, ps or html).

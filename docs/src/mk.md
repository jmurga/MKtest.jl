# MK approaches
We included other MK approaches in our module. All the functions use the formated SFS and divergence data described at [Summary statistics section](summstat.md).

## Standard MKT
The standard McDonald & Kreitman test ([McDonald and Kreitman, 1991]) was developed to be applied to protein coding sequences, combining both divergence ($D$) and polymorphism ($P$) sites, and categorizing mutations as synonymous ($P_S$, $D_S$) and non-synonymous ($P_N$, $D_N$). 

If all mutations are either strongly deleterious or neutral, then $Di/D0$ is expected to roughly equal $Pi/P0$. In contrast, if positive selection is operating in the region, adaptive mutations rapidly reach fixation and contribute more to divergence than polymorphism compared to neutral mutations, and then $Di/D0 > Pi/P0$. Assuming that adaptive mutations contribute little to polymorphism but substantially to divergence, the proportion of non-synonymous substitutions is inferred following Smith and Eyre-Walker (2002).

$\alpha = 1 - (\frac{P_N}{P_S}\cdot\frac{D_S}{D_N})$

Please check [`MKtest.standardMK`](@ref) to obtain more info.

```julia
adap = MKtest.parameters(n=661,cutoff=[0.0,1.0])
alpha, sfs, divergence, m= MKtest.parse_sfs(adap, data = "analysis/tgp.txt")

mk = MKtest.standardMK(sfs[1],divergence[1])
```

## Fay, Waycoff and Wu MK extension
[Fay et al. (2002)]([fwwMK](https://doi.org/10.1038/4151024a)) proposed an approach that removes all polymorphisms segregating at a frequency ($j$) below a given threshold (usually $j > 5\%–15\%$). Although there is no consensus about what this threshold should be used, [J. Charlesworth & Eyre-Walker (2008)](https://doi.org/10.1093/molbev/msn005) demonstrated that  estimates are robust using a frequency threshold of 15%, below which most slightly deleterious polymorphisms are found and removed. The estimates are reasonably accurate only when the rate of adaptive evolution is high and the Distribution of Fitness Effects (DFE) of deleterious mutations is leptokurtic.

$\alpha = 1 - (\frac{P_{N (j>15\%)}}{P_{S (j>15\%)}\cdot\frac{D_S}{D_N})$

Please check [`MKtest.fwwMK`](@ref) to obtain more info.

```julia
fww = MKtest.fwwMK(adap,sfs[1],divergence[1],cutoff=0.15)
```

## imputed MKT (in preparation)
The imputed MKT (impMKT) is a modification of the Fay, Waycoff, and Wu MK extension (fwwMK) ([Fay et al. (2002)]([fwwMK](https://doi.org/10.1038/4151024a)) to improve gene-by-gene analyses. The method propose the imputation of slightly deleterious mutations at the SFS rather than removing all variants below a frequency threshold. The imputedMK aims to maximize the information to test the excess of divergence ratio relative to polymorphism at the gene level.

$\alpha$ is estimated as

$\alpha_{imputed} = 1 - \left( \frac{P_{N} - P_{wd}}{P_{S}} \cdot \frac{D_{N}}{D_{S}} \right)$

where $P_{wd}$ is

$P_{wd} \approx P_{wd (j < 15\%)} = P_{N (j<15\%)} - \frac{P_{N (j>15\%) } \cdot P_{S (j<15\%)}}{P_{S (j>15\%)}}$

Please check [`MKtest.imputedMK`](@ref) to obtain more info.

```julia
imputed = MKtest.imputedMK(adap,sfs[1],divergence[1],cutoff=0.15)
```

## Asymptotic MKT
Proposed by Messer and Petrov (2013). This extension is robust to the presence of selective sweeps (genetic hitchhiking) and the segregation of slightly deleterious polymorphisms substitutions (BGS). In this approach, the authors defined $\alpha$ as a function that depends on the SFS of alleles. Therefore, $\alpha$ is estimated in different frequency intervals ($x$), and these values are then adjusted to an exponential function. An exponential fit is suitable as the non-synonymous allele frequency is expected to decay exponentially over the respective levels of synonymous polymorphisms (Messer & Petrov, 2013).

$\alpha$ is estimated as

$\alpha_{fit(x)} = a+b \cdot e^{-cx}$

Please check [`MKtest.aMK`](@ref) to obtain more info.

```julia
amk, ci, model = MKtest.aMK(adap,alpha[1])
```

## Grapes
Grapes is a ML method that can estimate the expected proportion of adaptive fixations given the inferred DFE from the MK data. Grapes assumes can model the DFE in the form of two different versions of the Fisher's geometric model, and a model assuming a Beta-shaped distribution of weak effect mutations, a Gamma distribution as well as a exponential distibution. Please check (Grapes repository)[https://github.com/BioPP/grapes] a cite (Galtier 2016)[https://doi.org/10.1371/journal.pgen.1005774] if you use `MKtest.grapes` function.

```julia
grapes_df = grapes(sfs,divergence,m,"GammaExpo","analysis/",20)
```
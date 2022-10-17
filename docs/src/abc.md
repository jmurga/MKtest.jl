# ABC inference from empirical data

At this point, you have a folder containing summary statistics and observed data to perform ABC inference. As explained in our [home page](index.md), we performed the ABC inference using [ABCreg](https://github.com/molpopgen/ABCreg). However, you can used other ABC software to perform the inference.

We link [ABCreg](https://github.com/molpopgen/ABCreg) with Julia to perform ABC inference. If you are going to use ABCreg to make inferences from our software directly, please [cite the publication](https://doi.org/10.1186/1471-2156-10-35). Remember you need to install ABCreg before continue. Please check [home page](index.md) to install ABCreg.

It is possible to perform the inference through Julia. We set the tolerance value such that 2500 acceptances were recorded for the inference

```julia
MKtest.ABCreg(analysis_folder="analysis/",S=size(adap.dac,1),tol=0.025,abcreg="/home/jmurga/ABCreg/src/reg");
```

The function will output one file per bootstrapped replicas containing the posteriors distributions. The posterior distributions contains five columns corresponding to :
 - α weak: Contribution of weak selecction to $\alpha$
 - α strong: Contribution of strong selecction to $\alpha$
 - α: Adaptation rate
 - γ: Negative selection coefficient
 - β: DFE shape parameter

We used R to estimate the Maximum-A-Posteriori (MAP) from posterior distributions following ABCreg examples. We linked Julia and R internally. The module contains functions to perform the estimations without quit the Julia session.

If you will perform MAP estimates and plot using our module, be sure you have installed R and the following packages: ggplot2 and data.table, locfit. 

```julia
MKtest.source_plot_map_r("analysis/script.jl")
posterior, tgp_map = MKtest.plot_map(analysis_folder="analysis/");
DataFrames.describe(tgp_map)
```

```
 Row │ variable  mean          min           median        max          nmissing  eltype   
     │ Symbol    Float64       Float64       Float64       Float64      Int64     DataType 
─────┼─────────────────────────────────────────────────────────────────────────────────────
   1 │ aw           0.108927     0.010857       0.108796      0.19754          0  Float64
   2 │ as           0.0506607   -0.00750128     0.0514826     0.134143         0  Float64
   3 │ a            0.152842     0.0962341      0.149083      0.233131         0  Float64
   4 │ gam_neg    1184.81       512.458       1277.47       1903.12             0  Float64
   5 │ shape        0.142934     0.128369       0.14189       0.167394         0  Float64
```

![image](https://raw.githubusercontent.com/jmurga/MKtest.jl/master/docs/src/figure2.svg)

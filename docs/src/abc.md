# ABC inference from empirical data

At this point, you have a folder containing summary statistics and observed data to perform ABC inference. As explained in our [home page](index.md), we performed the ABC inference using [ABCreg](https://github.com/molpopgen/ABCreg). However, you can used other ABC software to perform the inference.

We link [ABCreg](https://github.com/molpopgen/ABCreg) with Julia to perform ABC inference. If you are going to use ABCreg to make inferences from our software directly, please [cite the publication](https://doi.org/10.1186/1471-2156-10-35). Remember you need to install ABCreg before continue. Please check [home page](index.md) to install ABCreg.

It is possible to perform the inference through Julia. We set the tolerance value such that 2500 acceptances were recorded for the inference

```julia
posteriors = MKtest.ABCreg(analysis_folder="mktest/",S=size(adap.dac,1),tol=0.025,abcreg="/home/jmurga/ABCreg/src/reg");
```
The function will output one file per dataset containing the posteriors distributions. The posterior distributions contains five columns corresponding to:
 - α weak: Contribution of weak selecction to $\alpha$
 - α strong: Contribution of strong selecction to $\alpha$
 - α: Adaptation rate
 - γ: Negative selection coefficient
 - β: DFE shape parameter

You can check multiple statistics from posteriors distribution using `MKtest.summary_abc` function. Please check [```MKtest.summary_abc```](@ref), you can approximate the inference using different statistics such as the mode, the mean or the median from posterior.

```julia
df_summary, parameter_inference = MKtest.summary_abc(posteriors)
```

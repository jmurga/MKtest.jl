# Command-Line Interface

We develop a Command-Line Interface (CLI) in case you want to avoid Julia interpreter. You can easly download [abcmk_cli.jl](https://raw.githubusercontent.com/jmurga/MKtest.jl/master/scripts/abcmk_cli.jl). The CLI have different functions to perform the whole pipeline as explained in [Infering the rate and strength of adaptation](empirical.md) section

```bash
julia abcmk_cli.jl  
```

```bash
See --help of each command for usages  
  rates  
  summaries  
  abc
```

To reproduce the examples you can follow the steps described at [Empirical estimation](https://jmurga.github.io/MKtest.jl/dev/empirical/#Computational-pipeline-1) section

## Estimate rates
To perform the rate estimations you can use the function ```rates``` at the CLI. The function works simirarly to the function [```MKtest.rates```](@ref). You can check the argument at [Rates](https://jmurga.github.io/MKtest.jl/dev/rates/#Estimating-fixation-and-polymorphic-rates-considering-generalized-model-of-selection-and-linkage-1) sectino or using the macro ```@doc MKtest.rates``` in the Julia interpreter.

```julia
julia -t8 abcmk_cli.jl rates --help
Function to solve analytical fixation rates and the expected SFS. The function will create N random models from prior values. Use the arguments to defined the input range for each parameter.

If rho and/or theta are not set, default values will be used (0.001).

To parallelize the estimations please be sure you start up Julia using --threads/-t option and set the number of cores.

The function returns a HDF5 file containing models solved and rates. The rates will be used to compute summary statistics required at ABC.

Please check the documentation to get more info about models parameters or detailed arguments description https://jmurga.github.io/MKtest.jl/dev/cli/ to check model

Positional Arguments:
  N: Int64
  samples: Int64
  iterations: Int64

Optional Arguments:
  --alpha: String (default: 0.1,0.9)
  --gam-dfe: String (default: -1000,-200)
  --gam-flanking: String (default: -1000,-500)
  --gL: String (default: 5,10)
  --gH: String (default: 400,1000)
  --dac: String (default: 1,2,4,5,10,20,50,100,200,400,500,661,925,1000)
  --shape: Float64 (default: 0.184)
  --rho: Float64 (default: 0.001)
  --theta: Float64 (default: 0.001)
  --cutoff: String (default: 0.0,1.0)
  --output: String (default: rates.jld2)
```

```
julia -t8 abcmk_cli.jl rates 10000 661 100000 --alpha 0.1,0.9 --gam-dfe -1000,-200 --gam-flanking -1000,-500 --gL 1,10 --gH 200,2000 --output mktest/rates_cli.jld2
```
## Estimate summary statistics

To estimate summary statistics, we used the estimated analytical rates and empirical data following the description at section [Empirical estimation](empirical.md).


```bash
julia abcmk_cli summaries --help
Estimate summary statistics using observed data and analytical rates.

Positional Arguments:
  N: Int64
  samples: Int64
  data: String

Optional Arguments:
  --genes: String (default: <genes.txt>)
  --dac: String (default: 2,4,5,10,20,50,200,661,925)
  --cutoff: String (default: 0.0,1.0)
  --rates: String (default: rates.jld2)
  --folder: String (default: <folder>)
  --summsize: Int64 (default: 100000)
```

```bash
julia abcmk_cli.jl summaries --folder mktest/ --rates mktest/rates.jld2 --samples 661 --dac 2,4,5,10,20,50,200,661,925 --summsize 1000000
```

The function will output observed data (*alphas_N.txt*) and summary statistics (*summaries_N.txt*) in the selected folder. These file will be used at ABC inference to generate posterior distributions.

```bash
julia abcmk_cli.jl summaries --analysis_folder mktest/ --rates  mktest/rates.jld2 --samples 661 --replicas 100 --summstatSize 100000 --dac 2,4,5,10,20,50,200,661,925
```

## Perform ABC inference
At this point, you have a folder containing summary statistics and observed data to perform ABC inference. As explained in our [home page](index.md), we performed the ABC inference using [ABCreg](https://github.com/molpopgen/ABCreg). However, you can used other ABC software to perform the inference.

We link [ABCreg](https://github.com/molpopgen/ABCreg) with Julia to perform ABC inference. If you are going to use ABCreg to make inferences from our software directly, please [cite the publication](https://doi.org/10.1186/1471-2156-10-35). Remember you need to install ABCreg before continue. Please check [home page](index.md) to install ABCreg.

It is possible to perform the inference through Julia. The function will output one file per bootstrapped replicas containing the posteriors distributions. We set the tolerance value to record 1000 values for the regression.  The posterior distributions contains five columns corresponding to:

 - α weak: Contribution of weak selecction to $\alpha$
 - α strong: Contribution of strong selecction to $\alpha$
 - α weak: Adaptation rate
 - γ: Negative selection coefficient
 - β: Negative selection coefficient


```bash
julia abcmk_cli.jl abc --help

Positional Arguments:
  nsumm: Int64
  tol: Float64

Optional Arguments:
  --folder: String (default: <folder>)
  --abcreg: String (default: reg)

```


```bash
julia -t8 abc 9 0.025 --folder mktest/ --S 9 --abcreg /home/jmurga/ABCreg/src/reg
```

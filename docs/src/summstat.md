# Summary statistics
To estimate summary statistics, we used the estimated analytical rates and empirical data following the description at section [Empirical estimation](empirical.md)

Before starting the summary statistics, consider parallelizing the process using Julia Distributed computing. If you are following the tutorial step by step, do not input the following commands.

```julia
using MKtest
```

Before to compute the summary statistic, declare a model specifying the samples size and a DAC.

```julia
adap = MKtest.parameters(n=661,dac=[2,4,5,10,20,50,200,661,925])
```

Note you can only input DAC already estimated, nonetheles you can perform any subset from the estimated DAC. To check the estimated DAC you can follow the hierarchy of the ```h5``` variable.

```julia 
using JLD2
# Check hierarchy
h5   = jldopen("analysis/rates.jld2")
h5

JLDFile /home/jmurga/analysis/rates.jld2 (read-only)
 └─📂 1000
    └─📂 661
       ├─🔢 models
       ├─🔢 neut
       ├─🔢 sel
       ├─🔢 dsdn
       └─🔢 dac
```

```julia
# Checking estimated dac, string pattern inside the HDF5 variable
h5["1000/661/dac"]

14-element Vector{Int64}:
    1
    2
    4
    5
   10
   20
   50
  100
  200
  400
  500
  661
  925
 1000
```

The function ```MKtest.summary_statistics``` use the SFS and divergence Matrix output from ```MKtest.parse_sfs```. We include the argument ```bootstrap``` to perform bootstrap analysis following [polyDFE](https://github.com/paula-tataru/polyDFE) manual. In the following example we boostrap the SFS and divegence file 100 times subseting 10^5 summary statistic for each dataset:

```julia
@time summstat = MKtest.summary_statistics(param=adap,h5_file="analysis/rates.jld2",sfs=sfs,divergence=divergence,analysis_folder="analysis/",summstat_size=10^5,bootstrap=100);
```

Nonetheless, you can read your SFS and divergence files using the packages CSV and DataFrames to create the input Matrix. In such case, please be sure that sfs and divergence arguments are of type ```Vector``` using square braces (```[]```) comprehesion

```
using CSV, DataFrames
sfs = CSV.read("/path/to/sfs_file.txt", header=false, DataFrame) |> Matrix
divergence = CSV.read("/path/to/divergence_file.txt", header=false, DataFrame) |> Matrix

# Note we modified both variables into a Vector using square braces ([]) comprehesion
@time summstat = MKtest.summary_statistics(param=adap,h5_file="analysis/rates.jld2",sfs=[sfs],divergence=[divergence],analysis_folder="analysis/",summstat_size=10^5,replicas=100,bootstrap=true);

```

The function will create summary statistic files and the observed data files (*summaries_N.txt* and *alpha_N.txt* respectively). Both files will be used to perform the ABC inference. Each file will be used to infer $\alpha$ using the bootstrapped SFS.

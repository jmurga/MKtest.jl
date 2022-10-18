# Summary statistics
To estimate summary statistics, we used the estimated analytical rates and empirical data following the description at section [Empirical estimation](empirical.md)

Before executing the rate estimation, start-up Julia using `-t` option to add the desired number of threads to parallelize the estimation

```bash
julia -t8
```

Before to compute the summary statistic, declare a model specifying the sample size and cutoff used to estimate the rates.

```julia
adap = MKtest.parameters(N=10000,n=661,dac=[2,4,5,10,20,50,200,661,925],cutoff=[0.0,1.0])
```

Note you can only input DAC already estimated, nonetheles you can perform any subset from the estimated DAC. To check the estimated DAC you can follow the hierarchy of the ```h5``` variable.

```julia 
using JLD2
# Check hierarchy
h5   = jldopen("analysis/rates.jld2")
h5

JLDFile /home/jmurga/analysis/rates.jld2 (read-only)
 └─📂 10000
    └─📂 661
       └─📂 cutoff=[0.0,1.0]
          ├─🔢 models
          ├─🔢 neut
          ├─🔢 sel
          ├─🔢 dsdn
          └─🔢 dac
```

```julia
# Checking estimated dac, string pattern inside the HDF5 variable
h5["10000/661/dac"]

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

The function ```MKtest.summary_statistics``` use the SFS and divergence Matrix output from ```MKtest.parse_sfs```.

```julia
@time summstat = MKtest.summary_statistics(param=adap,h5_file="analysis/rates.jld2",sfs=sfs,divergence=divergence,analysis_folder="analysis/",summstat_size=10^5,bootstrap=100);
```

Nonetheless, you can read your SFS and divergence files using the packages `CSV` and `DataFrames` to create the input Vectors. In such case, please be sure that sfs and divergence arguments are of type ```Vector``` using square braces (```[]```) comprehesion

```
using CSV, DataFrames
sfs = CSV.read("/path/to/sfs_file.txt", header=false, DataFrame) |> Matrix
divergence = CSV.read("/path/to/divergence_file.txt", header=false, DataFrame) |> Matrix

# Note we modified both variables into a Vector using square braces ([]) comprehesion
@time summstat = MKtest.summary_statistics(param=adap,h5_file="analysis/rates.jld2",sfs=[sfs],divergence=[divergence],analysis_folder="analysis/",summstat_size=10^5);

```

The function will create $N$ summary statistic files and the observed data files depending on the length of variables `sfs` and `divergence` parsed with `MKtest.parse_sfs` function (*summaries_N.txt* and *alpha_N.txt* respectively, see [Data](data.md)). Both files will be used to perform the ABC inference.

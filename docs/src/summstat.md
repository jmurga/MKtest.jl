# Summary statistics
To estimate summary statistics, we need to used the estimated expected fixation rates and frequency spectra as well as empirical data.

The module includes functions to parse TGP from [Uricchio et al. (2019)](https://doi.org/10.1038/s41559-019-0890-6) and DGN from [Murga-Moreno et al. (2019)](https://doi.org/10.1093/nar/gkz372
). Please to parse TGP or DGN raw data into SFS and divergence counts, first download raw files deposited in our repository:  

 - [TGP](https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/tgp.txt)
 - [DGN Zambia population](https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/dgn_ral.txt)  
 - [DGN Raleigh population](https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/dgn_zi.txt)  

If you are going to parse DGN, you need to change the value of the argument *isoline* to *true*. Following the [Murga-Moreno et al. (2019)](https://doi.org/10.1093/nar/gkz372) sample size for each population is:

 - Zambia population: 154
 - RAL population: 160

## Parsing genomic data

Once you have downloaded the files, you can use the function ```MKtest.parse_sfs``` to convert the data into SFS and divergence counts. Please check [`MKtest.parse_sfs`](@ref) to get more info or execute:

```julia

download("https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/tgp.txt","mktest/tgp.txt")

adap = MKtest.parameters(n=661,cutoff=[0.0,1.0])
alpha, sfs, divergence, m= MKtest.parse_sfs(adap, data = "mktest/tgp.txt")
```

It is possible to directly subset genes IDs using the same ids deposited at you data. You can use a (column list file)[https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/ensembl_list.txt] or (CSV-like file)[https://raw.githubusercontent.com/jmurga/MKtest.jl/main/data/test_nonVIPs.txt] to subset the file list. If you use CSV-like each row will be parsed independtly.

```julia
download("https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/ensembl_list.txt","mktest/ensembl_list.txt")

alpha, sfs, divergence = MKtest.parse_sfs(adap, data = "mktest/tgp.txt",gene_list = "mktest/ensembl_list.txt")
```

```julia
download("https://raw.githubusercontent.com/jmurga/MKtest.jl/main/data/example_bootstrap.txt","mktest/example_bootstrap.txt")

# In our case, eachrow is a bootstrapped set
alpha, sfs, divergence, m = MKtest.parse_sfs(adap, data = "mktest/tgp.txt",gene_list = "mktest/example_bootstrap.txt")
```

## Estimating summary statistics
Once you have the SFS and divergence data variable you can compute the summary statistic. To do it declare a model specifying the sample size corresponding to your data as well as the cutoff used to estimate the rates.

```julia
adap = MKtest.parameters(N=1000,n=661,dac=[2,4,5,10,20,50,200,661,925],cutoff=[0.0,1.0])
```

Note you can only input DAC already estimated, nonetheles you can perform any subset from the estimated DAC. To check the estimated DAC you can follow the hierarchy of the hdf5 file.

```julia 
using JLD2
# Check hierarchy
h5   = jldopen("mktest/rates.jld2")
h5

JLDFile /home/jmurga/mktest/rates.jld2 (read-only)
 └─📂 1000
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

The function ```MKtest.summary_statistics``` will produce the simulated (summary statistics) and obversed data to run the ABC inference

```julia
@time df = MKtest.summary_statistics(adap,h5_file="mktest/rates.jld2",sfs=sfs,divergence=divergence,analysis_folder="mktest/",summstat_size=10^5);
```
The function will create $N$ summary statistic files and the observed data files depending on the length of variables `sfs` and `divergence` parsed with `MKtest.parse_sfs` function (*summaries_N.txt* and *alpha_N.txt* respectively, see [Data](data.md)). Both files will be used to perform the ABC inference.

You can also read any SFS and divergence files using the packages `CSV` and `DataFrames` to create the input Vectors. In such case, please be sure that sfs and divergence arguments are of type ```Vector``` using square braces (```[]```) comprehesion

```
using CSV, DataFrames
sfs = CSV.read("/path/to/sfs_file.txt", header=false, DataFrame) |> Matrix
divergence = CSV.read("/path/to/divergence_file.txt", header=false, DataFrame) |> Matrix

# Note we modified both variables into a Vector using square braces ([]) comprehesion
@time df = MKtest.summary_statistics(adap,h5_file="mktest/rates.jld2",sfs=[sfs],divergence=[divergence],analysis_folder="mktest/",summstat_size=10^5);
```


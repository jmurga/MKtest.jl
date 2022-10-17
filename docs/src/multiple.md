# Infering the adapation rate from multiple datasets
As we explained at section [Input Data](input.data) it is possible to subset the TGP dataset given a matrix of Ensembl/Flybase ids. In this section we will use whole-genome protein-coding information, Viral Interacting Proteins and Non-Viral Interacting Proteins to infer the empirical adaptation.

The next step assume that you execute the previous sections. Please before to execute this example, be sure you already performed the [rates estimates](rates.md) and your *analysis/* folder contain the TGP data.

We are going to create three separate directories into the *analysis/* directory we are working with to perform the estimations. You can create the folder manually, through bash, or directly in the Julia interpreter

```julia
using Distributed
addprocs(7)
@everywhere using MKtest
using CSV, DataFrames, JLD2

adap = MKtest.parameters(n=661,dac=[2,4,5,10,20,50,200,661,925])

mkpath("analysis/wg/")
mkpath("analysis/vips/")
mkpath("analysis/nonvips/")
```

Once you have the folder please download from our repository the files containing the VIPs and Non-VIPs ids.

```julia
download("https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/vip_list.txt","analysis/vips/vip_list.txt")
download("https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/nonvip_list.txt","analysis/nonvips/nonvip_list.txt")
```

Now we are going to parse our three dataset to output the SFS and divergence files into each analysis folder. We we follow the example provided at [Input Data](input.data) section. Once you have the data parsed, you can estimate the summary statistic following [Summary statistic](summstat.md) section. Load the rates and select a DAC before to continue the analysis. And finally you can perform the ABC inference using ABCreg or another ABC software using the files *alphas.txt* and *summaries.txt* deposited in each folder. We will perform the inference using ABCreg. Please compile ABCreg before to perform the execution as described in the [Installation](index.md)

 - Whole-Genome dataset
 
```julia
alpha, sfs, divergence = MKtest.parse_sfs(sample_size = 661, data = "analysis/tgp.txt")

@time summstat = MKtest.summary_statistics(param=adap,sfs=sfs,divergence=divergence,h5_file="analysis/rates.jld2",analysis_folder="analysis/wg/",summstat_size=10^5,bootstrap=100);

MKtest.ABCreg(analysis_folder="analysis/wg/",P=5,S=size(adap.dac,1),tol=0.025,abcreg="/home/jmurga/ABCreg/src/reg");
```

 - VIPs dataset
```julia
vips = String.(CSV.read("analysis/vips/vip_list.txt",header=false,DataFrame) |> Array)

alpha, sfs, divergence = MKtest.parse_sfs(sample_size = 661, data = "analysis/tgp.txt",gene_list=vips)

@time summstat = MKtest.summary_statistics(param=adap,sfs=sfs,divergence=divergence,h5_file="analysis/rates.jld2",analysis_folder="analysis/vips/",summstat_size=10^5,bootstrap=100);

MKtest.ABCreg(analysis_folder="analysis/vips/",P=5,S=size(adap.dac,1),tol=0.025,abcreg="/home/jmurga/ABCreg/src/reg");
```

 - Non-VIPs dataset
```julia
nonvips_list = String.(CSV.read("analysis/nonvips/nonvip_list.txt",DataFrame) |> Array)

alpha, sfs, divergence = MKtest.parse_sfs(sample_size = 661, data = "analysis/tgp.txt",gene_list=nonvips_list)

@time summstat = MKtest.summary_statistics(param=adap,sfs=sfs,divergence=divergence,h5_file="analysis/rates.jld2",analysis_folder="analysis/nonvips/",summstat_size=10^5,bootstrap=100);

MKtest.ABCreg(analysis_folder="analysis/nonvips/",P=5,S=size(adap.dac,1),tol=0.025,abcreg="/home/jmurga/ABCreg/src/reg");

```

Once the posterior distributions are estimated you can perform the Maximum-A-Posterior (MAP) estimates. We performed the MAP estimates following ABCreg examples. Please be sure you have installed R and the packages(ggplot2, locfit and data.table). This packages are already installed in Docker and Singularity images. Load the R functions to estimate and plot MAP executing the following command

```julia
MKtest.source_plot_map_r("analysis/script.jl")
```

 - Whole-Genome dataset
```julia
tgpmap = MKtest.plot_map(analysis_folder="analysis/wg/");
DataFrames.describe(tgpmap[2])

5×7 DataFrame
 Row │ variable  mean          min          median        max          nmissing  eltype   
     │ Symbol    Float64       Float64      Float64       Float64      Int64     DataType 
─────┼────────────────────────────────────────────────────────────────────────────────────
   1 │ aw           0.104377     0.0350986     0.110202      0.184566         0  Float64
   2 │ as           0.0481464   -0.004137      0.0441584     0.114982         0  Float64
   3 │ a            0.151626     0.110054      0.148903      0.203543         0  Float64
   4 │ gamNeg    1302.71       612.844      1474.7        1742.59             0  Float64
   5 │ shape        0.141999     0.130137      0.141801      0.15542          0  Float64


```

 - VIPs dataset
```julia
vipsmap = MKtest.plot_map(analysis_folder="analysis/vips/");
DataFrames.describe(vipsmap[2])

5×7 DataFrame
 Row │ variable  mean        min          median      max          nmissing  eltype   
     │ Symbol    Float64     Float64      Float64     Float64      Int64     DataType 
─────┼────────────────────────────────────────────────────────────────────────────────
   1 │ aw          0.131834   -0.016945     0.14184      0.284644         0  Float64
   2 │ as          0.160363    0.0780512    0.157498     0.27026          0  Float64
   3 │ a           0.284903    0.229213     0.283095     0.369928         0  Float64
   4 │ gamNeg    807.123     583.229      657.89      1594.79             0  Float64
   5 │ shape       0.202337    0.17996      0.201092     0.235535         0  Float64
```

 - Non-VIPs dataset
```julia
nonvipsmap = MKtest.plot_map(analysis_folder="analysis/nonvips/");
DataFrames.describe(nonvipsmap[2])

5×7 DataFrame
 Row │ variable  mean          min          median        max          nmissing  eltype   
     │ Symbol    Float64       Float64      Float64       Float64      Int64     DataType 
─────┼────────────────────────────────────────────────────────────────────────────────────
   1 │ aw           0.0982744    0.0320596     0.098271      0.196821         0  Float64
   2 │ as           0.0438451   -0.001272      0.0450831     0.099763         0  Float64
   3 │ a            0.140408     0.0964904     0.141531      0.196149         0  Float64
   4 │ gamNeg    1441.17       681.867      1547.24       1773.32             0  Float64
   5 │ shape        0.134675     0.12122       0.135073      0.148782         0  Float64

```

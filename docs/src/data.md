# Parsing genomic data
The module includes functions to parse TGP from Uricchio et al. (2019) and DGN from Murga-Moreno et al. (2019). In addition, the module have a function to parse SFS and divergence from multi-FASTA data following Murga-Moreno et al. (2019)

Please to parse raw data into SFS and divergence counts, first download raw files deposited in our repository:  

 - [TGP](https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/tgp.txt)
 - [DGN Zambia population](https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/dgn_ral.txt)  
 - [DGN Raleigh population](https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/dgnZi.txt)  

```julia
mkdir("analysis/")
download("https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/tgp.txt","analysis/tgp.txt")
download("https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/dgn_ral.txt","analysis/dgn_ral.txt")
```
## Parsing TGP and DGN data manually
Once you have downloaded the files, you can use the function ```MKtest.parse_sfs``` to convert the data into SFS and divergence counts. Please check [`MKtest.parse_sfs`](@ref) to get more info o execute:

```julia
alpha, sfs, divergence = MKtest.parse_sfs(sample_size = 661, data = "analysis/tgp.txt")
```

To save the data, you can use CSV and DataFrames packages

```julia
using CSV, DataFrames
CSV.write("analysis/tpg_sfs.tsv",DataFrame(sfs,:auto),delim='\t',header=false)
CSV.write("analysis/tpg_div.tsv",DataFrame(permutedims(divergence),:auto),delim='\t',header=false)
```

It is possible to directly subset genes IDs using Ensembl or Flybase id. Use a variable of type ```Matrix{String} or Vector{String}``` into the argument *gene_list*. The input ```Matrix{String}``` will be subset by columns creating *ncol* SFS and divergence Matrix.

```julia
download("https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/ensembl_list.txt","analysis/ensembl_list.txt")
ensembl_list = CSV.read("analysis/ensembl_list.txt",header=false,DataFrame) |> Array

alpha, sfs, divergence = MKtest.parse_sfs(sample_size = 661, data = "analysis/tgp.txt",gene_list = ensembl_list)
```

If you are going to perform a bootstrap analysis you can input a ```Matrix{String}``` to subset the bootstraped gene ids and create the SFS and divergence files.

```julia
download("https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/example_bootstrap","analysis/example_bootstrap.txt")

bootstrap_list = String.(Array(CSV.read("analysis/example_bootstrap.txt",DataFrame))[:,2:end])

# In our case, eachrow is a bootstrapped set
alpha, sfs, divergence = MKtest.parse_sfs(sample_size = 661, data = "analysis/tgp.txt",gene_list = bootstrap_list)
```


If you are going to parse DGN, you need to change the value of the argument *isoline* to *true*. Following the Murga-Moreno et al. (2019) sample size for each population is:

 - Zambia population: 154
 - RAL population: 160

```julia
alpha, sfs, divergence = MKtest.parse_sfs(sample_size = 160, data = "analysis/dgn_ral.txt",isolines=true)
```
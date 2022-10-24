# Parsing genomic data
The module includes functions to parse TGP from [Uricchio et al. (2019)](https://doi.org/10.1038/s41559-019-0890-6) and DGN from [Murga-Moreno et al. (2019)](https://doi.org/10.1093/nar/gkz372).

Please to parse raw data into SFS and divergence counts, first download raw files deposited in our repository:  

 - [TGP](https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/tgp.txt)
 - [DGN Zambia population](https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/dgn_ral.txt)  
 - [DGN Raleigh population](https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/dgnZi.txt)  

```julia
mkdir("analysis/")
download("https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/tgp.txt","analysis/tgp.txt")
```
## Parsing TGP and DGN data manually
Once you have downloaded the files, you can use the function ```MKtest.parse_sfs``` to convert the data into SFS and divergence counts. Please check [`MKtest.parse_sfs`](@ref) to get more info or execute:

```julia
adap = MKtest.parameters(n=661,cutoff=[0.0,1.0])
alpha, sfs, divergence, m= MKtest.parse_sfs(adap, data = "analysis/tgp.txt")
```

It is possible to directly subset genes IDs using the same ids deposited at you data. You can use a (column list file)[https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/ensembl_list.txt] or (CSV-like file)[https://raw.githubusercontent.com/jmurga/MKtest.jl/main/data/test_nonVIPs.txt] to subset the file list. If you use CSV-like each row will be parsed independtly.

```julia
download("https://raw.githubusercontent.com/jmurga/MKtest.jl/master/data/ensembl_list.txt","analysis/ensembl_list.txt")

alpha, sfs, divergence = MKtest.parse_sfs(adap, data = "analysis/tgp.txt",gene_list = "analysis/ensembl_list.txt")
```

```julia
download("https://raw.githubusercontent.com/jmurga/MKtest.jl/main/data/example_bootstrap.txt","analysis/example_bootstrap.txt")

# In our case, eachrow is a bootstrapped set
alpha, sfs, divergence, m = MKtest.parse_sfs(adap, data = "analysis/tgp.txt",gene_list = "analysis/example_bootstrap.txt")
```

If you are going to parse DGN, you need to change the value of the argument *isoline* to *true*. Following the [Murga-Moreno et al. (2019)](https://doi.org/10.1093/nar/gkz372) sample size for each population is:

 - Zambia population: 154
 - RAL population: 160

```julia
download("https://raw.githubusercontent.com/jmurga/MKtest.jl/main/data/dgn_ral.txt","analysis/dgn_ral.txt")
adap = MKtest.parameters(n=160,cutoff=[0.0,1.0])
alpha, sfs, divergence, m = MKtest.parse_sfs(adap, data = "analysis/dgn_ral.txt",isolines=true)
```
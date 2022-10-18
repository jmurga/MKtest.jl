# Input data

To estimate summary statistics, you need to provide empirical SFS and divergence files. As explained in section [data](data.md), you can directly parse TGP or DGN data using our module. Nonetheless, you can input any other SFS and divergence file.

You can easily use Julia to download the files using Julia or Bash. You also can download the files manually.

Once you have estimated the analytical rates and parsed the SFS and divergence information into variables, you can estimate the summary statistics.

```julia
adap = MKtest.parameters(n=661,cutoff=[0.0,1.0])
alpha, sfs, divergence, m= MKtest.parse_sfs(adap, data = "analysis/tgp.txt")
```
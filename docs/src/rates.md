# Estimating fixation and polymorphic rates considering generalized model of selection and linkage

Before executing the rate estimation, you need to load *Distributed* module and add some threads

```julia
using Distributed
addprocs(7)
```

Then you need to declare the ```MKtest``` module in all the threads using ```@everywhere``` macro. Otherwise, the ```MKtest``` module will perform the estimation just using the main core

```julia
@everywhere using MKtest
```

Declare a variable containing some basic information about your model. We used a sample size of 661 to perform later analysis over TGP data. The selected Derived Alleles Counts (DAC) will be used to compute summary statistics and perform ABC inference. It is possible to subset any of the selected DAC values when computing summary statistics.

```julia
adap = MKtest.parameters(n=661,dac=[1,2,4,5,10,20,50,100,200,400,500,661,925,1000])
```

Note that ```adap``` contains about mutation rate, recombination rates, DFE, BGS and probabilities of fixations. To check all the arguments you can access to the function documentation using ```@doc MKtest.parameter```

```julia
@doc MKtest.parameter

  Mutable structure containing the variables required to solve the analytical approach. All the functions are solve using the internal values of the structure. For this reason, adap is the
  only exported variable. Adap should be change before the perform the analytical approach, in other case, $\alpha_{(x)}$ will be solve with the default values.

  Parameters
  ≡≡≡≡≡≡≡≡≡≡≡≡

    •  gam_neg::Int64: Selection coefficient for deleterious alleles

    •  gL::Int64: Selection coefficient for weakly benefical alleles

    •  gH::Int64: Selection coefficient for strongly benefical alleles

    •  al_low::Float64: Proportion of α due to weak selection

    •  al_tot::Float64: α

    •  θ_noncoding::Float64: Mutation rate defining BGS strength

    •  θ_coding::Float64: Mutation rate on coding region

    •  al::Float64: DFE shape parameter

    •  be::Float64: DFE scale parameter

    •  B::Float64: BGS strength

    •  B_bins::Array{Float64,1}: BGS values to simulate

    •  ppos_l::Float64: Fixation probabily of weakly beneficial alleles

    •  ppos_h::Float64: Fixation probabily of strongly beneficial alleles

    •  N::Int64: Population size

    •  n::Int64: Sample size

    •  Lf::Int64: Flanking region length

    •  rho::Float64: Recombination rate

    •  TE::Float64
```

We will estimate the empirical adaptation rate using TGP data using the estimated DFE parameters at Boyko et al (2008). Nonetheless, shape and scale DFE parameter are flexible in our model.

```julia
adap.al = 0.184
adap.be = abs(0.184/-457)
```

Now the variable ```adap``` contains sample size, DAC and DFE information. The function ```MKtest.rates``` will perform the analytical estimation of *N* independent models regarding DFE, BGS, mutation rate, and recombination rate. In the following example, we used the function ```MKtest.rates``` to input the prior distributions. The function will randomize the input values to solve *N* independent estimation regarding our model. 

```julia
@doc MKtest.rates
  Function to solve randomly N scenarios. The function will create N models, defined by
  MKtest.parameters(), to estimate analytically fixation and polymorphic rates for each
  model. The rates will be used to compute summary statistics required at ABC. The function
  output a HDF5 file containing the solved models, the selected DAC and the analytical
  rates.

  If ρ and/or theta are set to nothing, the function will input random values given the
  range 0.0005:0.0005:0.01. Otherwise you can fix the values.

  If gL is set to nothing, the function will not account the role of the weakly selected
  alleles in the estimation.

  Arguments
  ≡≡≡≡≡≡≡≡≡≡≡

    •  param::parameters: mutable structure containing the model

    •  gH::Array{Int64,1} : Range of strong selection coefficients

    •  gL::Union{Array{Int64,1},Nothing}: Range of weak selection coefficients

    •  gam_neg::Array{Int64,1} : Range of deleterious selection coefficients

    •  theta::Union{Float64,Nothing} : Population-scaled mutation rate on coding region

    •  rho::Union{Float64,Nothing} : Population-scaled recombination rate

    •  iterations::Int64 : Number of solutions

    •  output::String : File to output HDF5 file

  Returns
  ≡≡≡≡≡≡≡≡≡

    •  Array: summary statistics

    •  Output: HDF5 file containing models solved and rates.

```

```julia
@time df = MKtest.rates(param = adap,gH=200:2000,gL=1:10,gam_neg=-2000:-200,iterations = 10^5,output="analysis/rates.jld2");
```

The function will create a HDF5 file containing the solved models, fixation rates, polymorphic rates, and the selected DAC. This information will be used later to estimate summary statistics.


Note that [```MKtest.rates```](@ref) is the most resource and time-consuming function. In our case, the function will estimate 10^5 independent models. Each model solves the estimation for all possible BGS values. We used BGS values from 0.1 to 0.999 in 5% increments. In total, the example will produce 3.7 million estimates. We have used a hierarchical data structure (HDF5) to facilitate model parameters and rates storage.

The following example took about 1.5 hours to execute on the hardware described at section [Infering the rate and strength of adaptation](empirical.md)

If you have a system with few resources, it is possible [to download pre-computed TGP and DGN data rates](https://imkt.uab.cat/files/inputs/rates.jld2). Please go to the next section to continue the inference using the pre-computed rates.
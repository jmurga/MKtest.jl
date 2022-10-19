# Estimating expected fixation rates and frequency spectra considering generalized model of selection and linkage

The first step is to solve *N* genetic models to obtain the analytical estimations. The process is automatize to create random combination model parameters given priors distributions. The expected fixation rates and frequency spectra, polymorphic rates, and model information will be used to estimate informative summary statistics needed to perform ABC inference. The following section show how to define a genetic a model, input prior distributions and solve the models.

Before executing the rates estimation, start-up Julia using `-t` option to add the desired number of threads to parallelize the process

```bash
julia -t8
```

```julia
using MKtest
mkpath("analysis/")
```

You must to declare a variable containing some basic information about your model using the function ```MKtest.paramterers```. Note that ```MKtest.parameters``` contains information about mutation rate, recombination rate, DFE, BGS and probabilities of fixations. To check all the arguments you can access to the function documentation using ```@doc MKtest.parameter```

```julia
@doc MKtest.parameters
  Mutable structure containing the variables required to solve the analytical
  approach. All the functions are solve using the internal values of the structure.
  You should declare a mutable structure to the perform the analytical estimations.

  Parameters
  ≡≡≡≡≡≡≡≡≡≡≡≡

    •  gam_flanking::Int64: Selection coefficient for deleterious alleles

    •  gL::Int64: Selection coefficient for weakly benefical alleles

    •  gH::Int64: Selection coefficient for strongly benefical alleles

    •  al_low::Float64: Proportion of α due to weak selection

    •  al_tot::Float64: α

    •  θ_flanking::Float64: Mutation rate defining BGS strength

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

    •  ρ::Float64: Recombination rate

    •  TE::Float64

```


We will estimate the empirical adaptation rate using TGP data using the estimated DFE parameters at [Boyko et al (2008)](https://doi.org/10.1371/journal.pgen.1000083). Hence, we used a sample size of 661 to perform the analysis and the infered deleterious DFE. Note that shape DFE parameter is flexible and modified by a factor of 4 and scale DFE parameter is modified from a prior distribution when executing ```MKtest.rates``` function. In addition we will input the selected Derived Alleles Counts (DAC) to later compute summary statistics and perform ABC inference. It is possible to subset any of the selected DAC values when computing summary statistics. If you want to exclude any variant bellow or above a frequency threshold you can use the argument `cutoff`.


```julia
adap = MKtest.parameters(N=10000,n=661,dac=[1,2,4,5,10,20,50,100,200,400,500,661,925,1000],gam_dfe=-457,shape=0.184,cutoff=[0.0,1.0])
```

Now the variable ```adap``` contains sample size, DAC and deleterious DFE information. The function ```MKtest.rates``` will perform the analytical estimation of *N* independent models regarding DFE, BGS, mutation rate, and recombination rate. Please note you can edit such values when using ```MKtest.parameters```.

Now you can used the function ```MKtest.rates``` to input the prior distributions. The function will randomize the input values to solve *N* independent estimation regarding our model.

```julia
  rates(param,gH,gL,gam_flanking,gam_dfe,alpha,iterations,output)

  Function to solve randomly N scenarios. The function will create N models,
  defined by MKtest.parameters(), to solve analytically fixation rates and the
  expected SFS for each model. The rates will be used to compute summary statistics
  required at ABC inference. The function output a HDF5 file containing the solved
  models, the selected DAC and the analytical solutions.

  Arguments
  ≡≡≡≡≡≡≡≡≡≡≡

    •  param::parameters: mutable structure containing the model.

    •  gH::Array{Int64,1}: Range of strong selection coefficients.

    •  gL::Union{Array{Int64,1},Nothing}: Range of weak selection coefficients.

    •  gam_flanking::Array{Int64,1}: Range of deleterious selection coefficients at the flanking region.

    •  gam_dfe::Array{Int64,1}: Range of deleterious selection coefficients at the coding region.

    •  alpha::Vector{Float64}: Range of α value to solve

    •  iterations::Int64: Number of solutions.

    •  output::String: File to output HDF5 file.

  Returns
  ≡≡≡≡≡≡≡≡≡

    •  DataFrame: models solved.

    •  Output: HDF5 file containing models solved and rates.


```

```julia
@time df = MKtest.rates(adap,gH=[200,2000],gL=[1,10],gam_dfe=[-2000,-200],gam_flanking=[-1000,-500],iterations = 10,output="analysis/rates.jld2");
```

The function will create a HDF5 file containing the solved models, the expected fixation rates and frequency spectra, and the selected DAC. This information will be used later to estimate summary statistics.

Note that [```MKtest.rates```](@ref) is the most resource and time-consuming function. In our case, the function will estimate 10^5 independent models. Each model solves the estimation for all possible BGS values. We used BGS values from 0.1 to 0.999 in 5% increments (defined at `adap.B_range`). In total, the example will produce 3.7 million estimates. We have used a hierarchical data structure (HDF5) to facilitate model parameters and rates storage.

The following example took about 1.5 hours to execute on the hardware described at section [Infering the rate and strength of adaptation](empirical.md)

# Model parameters

*adap* is the only variable exported from *MKtest* module. It is a Mutable structure contaning the variables required to solve the analytical approach. Any value can be easly changed. Remember *adap* should be change before the execution, in other case, $\alpha_{(x)}$ will be solve with the default values. To change all the values at once, you can use [`MKtest.parameters`](@ref) in order to set specific models in a new variable.

```@docs
MKtest.parameters
MKtest.analytical_alpha
MKtest.Br
MKtest.set_θ
MKtest.binom_op
MKtest.Φ
```

# MKtest estimation
## Fixations
```@docs
MKtest.fix_neut
MKtest.fix_neg
MKtest.p_fix
MKtest.fix_pos_sim
```

## Polymorphism
```@docs
MKtest.sfs_neut
MKtest.sfs_pos
MKtest.sfs_neg
MKtest.cumulative_sfs
MKtest.reduce_sfs
```

## Rates
```@docs
MKtest.rates
```
## Summary statistics
```@docs
MKtest.poisson_fixation
MKtest.poisson_polymorphism
MKtest.sampling_summaries
MKtest.summary_statistics
```

## Inference tools
```@docs
MKtest.parse_sfs
MKtest.ABCreg
MKtest.summary_abc
```

## MKtest
```@docs
MKtest.aMK
MKtest.imputedMK
MKtest.fwwMK
MKtest.standardMK
MKtest.grapes
```

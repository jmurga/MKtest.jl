# Model parameters

*adap* is the only variable exported from *MKtest* module. It is a Mutable structure contaning the variables required to solve the analytical approach. Any value can be easly changed. Remember *adap* should be change before the execution, in other case, $\alpha_{(x)}$ will be solve with the default values. To change all the values at once, you can use [`MKtest.parameters`](@ref) in order to set specific models in a new variable.

```@docs
MKtest.parameters
MKtest.analytical_alpha
MKtest.assertion_params
MKtest.Br
MKtest.set_θ
MKtest.alpha_exp_sim_low
MKtest.alpha_exp_sim_tot
MKtest.solv_eqns
MKtest.set_ppos
MKtest.binom_op
MKtest.Φ
MKtest.analytical_alpha
```

# MKtest estimation
## Fixations
```@docs
MKtest.fix_neut
MKtest.fix_neg_b
MKtest.p_fix
MKtest.fix_pos_sim
```

## Polymorphism
```@docs
MKtest.sfs_neut
MKtest.sfs_pos
MKtest.sfs_pos_float
MKtest.sfs_neg
MKtest.cumulative_sfs
MKtest.reduce_sfs
```

## Rates
```@docs
MKtest.simulate_models
MKtest.solve_model
MKtest.solve_rates
MKtest.rates
```
## Summary statistics
```@docs
MKtest.poisson_fixation
MKtest.poisson_polymorphism
MKtest.sampling_summaries
MKtest.summary_statistics
MKtest.filter_expected
MKtest.pol_correction
```

## Inference tools
```@docs
MKtest.parse_sfs
MKtest.get_pol_div
MKtest.data_to_poisson
MKtest.write_files
MKtest.ABCreg
MKtest.get_mode
MKtest.estimates
MKtest.summary_abc
MKtest.bootstrap_data
```

## MKtest
```@docs
MKtest.aMK
MKtest.imputedMK
MKtest.fwwMK
MKtest.standardMK
MKtest.grapes
```

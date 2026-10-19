# Reference

Nothing is exported: call every function as `MKtest.f`.

## Input data

```@docs
MKtest.parse_sfs
MKtest.project
MKtest.fold
MKtest.cumulative_sfs
MKtest.reduce_sfs
MKtest.bootstrap_data
```

## MK tests

```@docs
MKtest.standard_mk
MKtest.fww_mk
MKtest.imputed_mk
MKtest.asymptotic_mk
MKtest.ExponentialFit
MKtest.LinearFit
MKtest.α_x
MKtest.dos
MKtest.ni
MKtest.ni_cmh
MKtest.grapes
```

## dfe_ml

```@docs
MKtest.dfe_ml
MKtest.dfe_model_comparison
MKtest.DFEResult
MKtest.params_wide
MKtest.param_values
MKtest.expected_sfs
MKtest.rj
MKtest.write_dfe_result
MKtest.read_dfe_result
MKtest.from_legacy_name
MKtest.grapes_param_name
MKtest.DFECache
```

### DFE models

```@docs
MKtest.AbstractDFEModel
MKtest.GammaZero
MKtest.GammaExpo
MKtest.GammaGamma
MKtest.ExpoZero
MKtest.ExpoExpo
MKtest.ExpoGamma
MKtest.DisplGamma
MKtest.ScaledBeta
MKtest.LogNormalZero
MKtest.NeuGamma
MKtest.FGMBesselK
MKtest.FixedSpikes
MKtest.spike_weights
MKtest.WithLethal
MKtest.dfe_parameters
MKtest.distribution
MKtest.dfe_density!
MKtest.dfe_log_scale
MKtest.dfe_start_points
MKtest.dfe_name
MKtest.discrete_components
MKtest.dfe_mass_below
MKtest.dfe_mass_above
MKtest.dfe_cdf_continuous
MKtest.dfe_cdf_continuous!
MKtest.dfe_ne_multipliers
```

### Axes

```@docs
MKtest.AbstractSFSExpectation
MKtest.PRFSojourn
MKtest.expectation_counts
MKtest.neutral_shape
MKtest.uses_exterior
MKtest.default_s_grid
MKtest.default_demography
MKtest.AbstractDominance
MKtest.Additive
MKtest.Dominance
MKtest.DominanceFunction
MKtest.effective_s
MKtest.AbstractDemography
MKtest.NuisanceRj
MKtest.Equilibrium
MKtest.free_rj
MKtest.AbstractNeModel
MKtest.ConstantNe
MKtest.AncestralNeRatio
MKtest.ne_ratio
MKtest.AbstractDivergenceModel
MKtest.SaturatedDivergence
MKtest.DFEConstrainedDivergence
MKtest.NoDivergence
MKtest.uses_divergence
MKtest.divergence_free
MKtest.AbstractLikelihoodKernel
MKtest.PoissonKernel
MKtest.NegativeBinomialKernel
MKtest.log_kernel
MKtest.kernel_parameters
MKtest.kernel_start
MKtest.AbstractIntegrator
MKtest.RombergIntegrator
MKtest.GaussLegendreSojourn
MKtest.ExactSojourn
MKtest.integrate
MKtest.expected_count
MKtest.dn_correction_term
MKtest.AbstractDFEQuadrature
MKtest.MidpointDensity
MKtest.IntervalCDF
MKtest.fill_dfe_vals!
MKtest.AbstractSGrid
MKtest.GrapesGrid
MKtest.LogGammaGrid
MKtest.s_grid_points
MKtest.s_grid_widths
MKtest.s_cell_bounds
MKtest.AbstractOptimizer
MKtest.LBFGSOptimizer
MKtest.NelderMeadOptimizer
MKtest.run_optimizer
MKtest.AbstractAlphaSummary
MKtest.DFEMAlpha
MKtest.summarize
MKtest.FitContext
```

## ABC-MK

```@docs
MKtest.parameters
MKtest.rates
MKtest.load_rates
MKtest.RatesTable
MKtest.valid_data
MKtest.summary_statistics
MKtest.abcmk
MKtest.RjSampling
MKtest.IndependentSampling
MKtest.FlatAlphaPrior
MKtest.TableRowPrior
MKtest.ABCreg
MKtest.abc
MKtest.abc_loclinear
MKtest.abcreg_loclinear
MKtest.summary_abc
MKtest.analytical_alpha
MKtest.set_ppos!
MKtest.set_θ!
MKtest.Φ
MKtest.binom_op
MKtest.fix_neut
MKtest.fix_neg
MKtest.p_fix
MKtest.fix_pos_sim
MKtest.sfs_neut
MKtest.sfs_pos
MKtest.sfs_neg
```

## Case and control sets

```@docs
MKtest.bootstrap_parameters
MKtest.bootstrap
MKtest.get_distance
MKtest.get_bootstrap
MKtest.case_control_table
MKtest.empirical_test
MKtest.empirical_pvalues
```

## Balancing selection

```@docs
MKtest.α_b
MKtest.parse_private_share
```

## Plots

These load with `using Plots` (`plot_α_x`, `plot_asymptotic_mk`) or `using StatsPlots`
(`plot_violins_empirical`).

```@docs
MKtest.plot_α_x
MKtest.plot_asymptotic_mk
MKtest.plot_violins_empirical
```

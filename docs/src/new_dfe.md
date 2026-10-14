# How to add a new DFE

A DFE in `dfe_ml` is a subtype of `MKtest.AbstractDFEModel`. It needs two methods in the
simplest case, and the engine derives everything else from them.

## 1. A DFE from a distribution

```julia
using MKtest, Distributions
import MKtest: AbstractDFEModel

struct NormalDFE <: AbstractDFEModel end

MKtest.dfe_parameters(::NormalDFE) = (μ = (-1e4, 1e5), σ = (1, 1e5))   # name = (lower, upper)
MKtest.distribution(::NormalDFE, μ, σ) = truncated(Normal(μ, σ), 0, Inf)

fit = MKtest.dfe_ml(sfs, div, NormalDFE())
fit.μ, fit.σ, fit.alpha
```

- `dfe_parameters` returns a NamedTuple: the keys are the parameter names (they become
  output columns and the names `fixed` and `bounds` accept), the values are the bounds,
  and the order is the parameter order.
- `distribution` receives the parameter values in that order and returns the law of |S|
  of deleterious mutations (S < 0), as any `Distributions.UnivariateDistribution` with
  support on [0, Inf).

From these two methods the engine builds:

- the density on the S grid: `pdf(d, -S)` for S < 0, 0 otherwise;
- the tail masses (`cdf`, `ccdf`) and the CDF, so `quadrature = IntervalCDF()` works;
- the optimizer scale: a parameter is optimised on a log scale when its lower bound is
  positive (`σ` here, not `μ`);
- three starting points, at 1/4, 1/2 and 3/4 of every range on that scale. The engine
  always adds the centre of the box, plus `nb_rand_start` seeded random points.

### Why the `MKtest.` prefix

`dfe_parameters` and `distribution` are functions of `MKtest`, and you are adding methods
to them. Julia requires the qualified name (`MKtest.dfe_parameters(...) = ...`) or an
explicit `import MKtest: dfe_parameters`. A bare `dfe_parameters(::NormalDFE) = ...`
would define a new, unrelated function in your module, and `dfe_ml` would never see it.

## 2. A wrapper that changes a density

A wrapper writes the density itself with `dfe_density!` instead of `distribution`. This
one adds a lethal class to any DFE:

```julia
struct WithLethal{M} <: AbstractDFEModel
    dfe::M
end

MKtest.dfe_parameters(m::WithLethal) = (MKtest.dfe_parameters(m.dfe)..., p_lethal = (0.0, 0.99))
MKtest.dfe_density!(out, m::WithLethal, p, s) =
    (MKtest.dfe_density!(out, m.dfe, p[1:end-1], s); out .*= 1 - p[end])

MKtest.dfe_ml(sfs, div, WithLethal(MKtest.GammaZero()))
```

`dfe_density!(out, m, p, s)` writes the density of `m` with parameters `p` at every S of
`s` into `out`. Without `distribution`, the tail masses and the CDF fall back to
numerical quadrature of this density.

`MKtest.WithLethal` is the built-in version of the same wrapper. It is written without
allocations (a view instead of `p[1:end-1]`), keeps the base model's log scales and
starting points, and refuses to wrap a model that already has `p_lethal`.

## 3. Optional methods

| method | default | override when |
|---|---|---|
| `MKtest.dfe_log_scale(m)` | lower bound > 0 | a positive parameter is better optimised linearly (shapes, proportions) |
| `MKtest.dfe_start_points(m)` | 3 points across the box | you know good starts (the built-ins use GRAPES') |
| `MKtest.dfe_name(m)` | the type name | the output table should say something else |
| `MKtest.discrete_components(m, p)` | none | the DFE has point masses, as `(S, weight)` pairs (`NeuGamma`, `FixedSpikes`) |
| `MKtest.dfe_mass_below(m, p, b)`, `dfe_mass_above(m, p, a)`, `dfe_cdf_continuous!(out, m, p, edges)` | from `distribution`, else numerical | a closed form is faster or more accurate |
| `MKtest.dfe_ne_multipliers(m, r)` | rescale parameters named `negGmean`, `posGmean`, `s0`, `negLNmean`, `BKscale` by r | the model is used with an Ne ratio other than 1 and has other S-scale parameters |

Write the trailing arguments of an overriding method untyped (as in `dfe_density!(out,
m::WithLethal, p, s)`): the engine's generic fallbacks are untyped beyond the model, so
the method that names your type is always the more specific one.

## 4. Check it

```julia
m = NormalDFE()
MKtest.dfe_param_names(m)                                    # [:μ, :σ]
lay = MKtest.build_layout(m, MKtest.PoissonKernel(), true)   # bounds, scales, nuisance slots
s = collect(range(-1000.0, 0.0; length = 101)); out = similar(s)
MKtest.dfe_density!(out, m, [100.0, 300.0], s)               # the density you meant
MKtest.dfe_ml(sfs, div, m; quadrature = MKtest.IntervalCDF())
```

`test/ml/extension.jl` runs these examples, including a check that the generic density
equals the hand-written `pdf` and that the `IntervalCDF` cell masses equal differences of
`cdf`.

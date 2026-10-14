# =============================================================================
# Public entry points: dfe_ml and dfe_model_comparison
# =============================================================================

_sfs_list(sfs::AbstractMatrix) = [sfs]
function _sfs_list(sfs::AbstractVector)
    all(x -> x isa AbstractMatrix, sfs) ||
        throw(ArgumentError("sfs must be a Matrix [DAC Pn Ps] or a Vector of such matrices"))
    return collect(sfs)
end

_div_list(d::AbstractVector{<:Real}) = [Float64.(d)]
_div_list(d::AbstractMatrix{<:Real}) = [Float64.(vec(d))]
function _div_list(d::AbstractVector)
    return [x isa AbstractArray{<:Real} ? Float64.(vec(x)) :
            throw(ArgumentError("divergence must be [Dn, Ds, Ln, Ls] or a Vector of those")) for x in d]
end

_axis_string(x) = repr(x)

"""
    dfe_ml(sfs, divergence, model::AbstractDFEModel; kwargs...) -> DataFrame

Maximum-likelihood fit of the DFE `model` and alpha to the SFS and divergence of one
or several datasets, with GRAPES' model (Galtier 2016). One row per dataset.

# Arguments
 - `sfs`: a matrix `[DAC Pn Ps]` with every bin `1..n-1` (unfolded; `n` chromosomes),
   as `parse_sfs` returns it, or a Vector of such matrices.
 - `divergence`: `[Dn, Ds, Ln, Ls]`, or a Vector of those (one per SFS).
 - `model`: an `AbstractDFEModel` instance (`GammaZero()`, `GammaExpo()`,
   `WithLethal(GammaZero())`, a user model, ...).

# Keywords (data)
 - `project_sfs::Union{Nothing,Int} = nothing`: project the SFS to this many diploid
   individuals (`2 project_sfs` chromosomes) with `project` first.
 - `fold::Bool = false`: fold the SFS (`syno_orientation_error` is then not fitted).
 - `Ln_poly`, `Ls_poly`: site counts of the polymorphism data, when they differ from
   the divergence ones (default: `Ln`, `Ls`).

# Keywords (model axes; each is independent of the others)
 - `ne_model::AbstractNeModel = ConstantNe()`: `AncestralNeRatio(r)` or
   `AncestralNeRatio(:free)` (adds the column `Ne_ratio`).
 - `div_model::AbstractDivergenceModel = SaturatedDivergence()`:
   `DFEConstrainedDivergence()` (GRAPES `-no_div_param`) or `NoDivergence()`
   (`-no_div_data`).
 - `demography::Union{Nothing,AbstractDemography} = nothing`: `NuisanceRj()` (GRAPES'
   per-bin `r_j`) or `Equilibrium()`; `nothing` takes `default_demography(expectation)`.
 - `dominance::AbstractDominance = Additive()`: `Dominance(h)`, `DominanceFunction(f)`.
 - `expectation::AbstractSFSExpectation = PRFSojourn()`.
 - `integrator::AbstractIntegrator = ExactSojourn()`; `RombergIntegrator()` is
   GRAPES' own rule (use it to compare with the binary).
 - `quadrature::AbstractDFEQuadrature = MidpointDensity()`; or `IntervalCDF()`.
 - `s_grid::Union{Nothing,AbstractSGrid} = nothing`: `nothing` takes
   `default_s_grid(expectation)` (`GrapesGrid()`).
 - `kernel::AbstractLikelihoodKernel = PoissonKernel()`; or `NegativeBinomialKernel()`.
 - `optimizer::AbstractOptimizer = LBFGSOptimizer()`; or `NelderMeadOptimizer()`.
 - `summaries = (DFEMAlpha(), FWWAlpha())`: the `AbstractAlphaSummary` column sets.

# Keywords (fit)
 - `fixed::NamedTuple = (;)`: pin parameters by name, e.g. `(pos_prop = 0.0,)`.
 - `bounds::NamedTuple = (;)`: replace bounds by name, e.g. `(pos_prop = (1e-6, 0.5),)`.
 - `no_syn_orient_error::Bool = false`: do not fit `syno_orientation_error`.
 - `nb_rand_start::Int = 0`: extra random starting points (seeded, reproducible).
 - `nearly_neutral::Real = 5`: S threshold of adaptive substitutions in alpha.
 - `FWW_threshold::Real = 0.15`: minimum frequency of the FWW alpha.
 - `alpha_ci_threshold::Real = 0.0`: log-likelihood drop of the alpha CI (GRAPES uses
   1.92); 0 disables it and the CI columns are NaN.
 - `cache::DFECache = DFECache()`: precomputed tables; pass one to reuse them across
   calls (safe across threads).

# Returns
A `DataFrame`, one row per dataset, with the columns `dataset`, `model`, `n`
(chromosomes), `converged`, `loglik`, `aic`, `n_params`, `free_params` (the
optimised names, comma separated), `theta`, one column per parameter (DFE, then
`syno_orientation_error`, then the kernel's), `Ne_ratio` under `AncestralNeRatio(:free)`,
the summary columns (`alpha`, `alpha_ci_low`, `alpha_ci_high`, `omega_a`, ...,
`omega_obs`, `alpha_div`, `alpha_dfe`, `alpha_fww`), and one column per axis setting
(`expectation`, `integrator`, `quadrature`, `s_grid`, `dominance`, `demography`,
`ne_model`, `div_model`, `kernel`, `optimizer`).
"""
function dfe_ml(sfs, divergence, model::AbstractDFEModel;
                project_sfs::Union{Nothing,Int} = nothing,
                fold::Bool = false,
                Ln_poly::Real = NaN,
                Ls_poly::Real = NaN,
                ne_model::AbstractNeModel = ConstantNe(),
                div_model::AbstractDivergenceModel = SaturatedDivergence(),
                demography::Union{Nothing,AbstractDemography} = nothing,
                dominance::AbstractDominance = Additive(),
                expectation::AbstractSFSExpectation = PRFSojourn(),
                integrator::AbstractIntegrator = ExactSojourn(),
                quadrature::AbstractDFEQuadrature = MidpointDensity(),
                s_grid::Union{Nothing,AbstractSGrid} = nothing,
                kernel::AbstractLikelihoodKernel = PoissonKernel(),
                optimizer::AbstractOptimizer = LBFGSOptimizer(),
                summaries = (DFEMAlpha(), FWWAlpha()),
                fixed::NamedTuple = NamedTuple(),
                bounds::NamedTuple = NamedTuple(),
                no_syn_orient_error::Bool = false,
                nb_rand_start::Int = 0,
                nearly_neutral::Real = 5,
                FWW_threshold::Real = 0.15,
                alpha_ci_threshold::Real = 0.0,
                cache::DFECache = DFECache())
    sfs_l = _sfs_list(sfs)
    div_l = _div_list(divergence)
    length(sfs_l) == length(div_l) ||
        throw(ArgumentError("got $(length(sfs_l)) SFS and $(length(div_l)) divergence vectors"))
    if project_sfs !== nothing
        sfs_l = [project(Matrix{Float64}(x), project_sfs) for x in sfs_l]
    end
    dem = demography === nothing ? default_demography(expectation) : demography
    grid = s_grid === nothing ? default_s_grid(expectation) : s_grid

    rows = Vector{OrderedDict{Symbol,Any}}(undef, length(sfs_l))
    for i in eachindex(sfs_l)
        n_i = Int(round(maximum(sfs_l[i][:, 1]))) + 1
        pre = get_precomputed(cache, n_i; integrator = integrator, grid = grid,
                              dominance = dominance, expectation = expectation)
        data = build_data_input(sfs_l[i], div_l[i], fold, uses_divergence(div_model);
                                Ln_poly = Float64(Ln_poly), Ls_poly = Float64(Ls_poly),
                                neutral_shape = pre.neutral_shape)
        row = fit_model(data, pre, model; kernel = kernel, quad = quadrature, optimizer = optimizer,
                        ne_model = ne_model, div_model = div_model, demography = dem,
                        no_syn_orient_error = no_syn_orient_error, fixed = fixed, bounds = bounds,
                        nb_rand_start = nb_rand_start, nearly_neutral = Float64(nearly_neutral),
                        FWW_threshold = Float64(FWW_threshold),
                        alpha_ci_threshold = Float64(alpha_ci_threshold), summaries = summaries)
        out = OrderedDict{Symbol,Any}(:dataset => i)
        merge!(out, row)
        for (k, v) in (:expectation => expectation, :integrator => integrator,
                       :quadrature => quadrature, :s_grid => grid, :dominance => dominance,
                       :demography => dem, :ne_model => ne_model, :div_model => div_model,
                       :kernel => kernel, :optimizer => optimizer)
            _put!(out, k, _axis_string(v))
        end
        rows[i] = out
    end
    return DataFrame([NamedTuple{Tuple(keys(r))}(Tuple(values(r))) for r in rows])
end

# -----------------------------------------------------------------------------
# Model comparison
# -----------------------------------------------------------------------------

"""
    _aic_lrt_tables(names, n_params, logliks, is_nested; extra_cols = (), skip)
        -> (aic::DataFrame, lrt::DataFrame)

AIC with Akaike weights for every model, and a likelihood-ratio test for every pair
with `is_nested(i, j)` (model `i` a special case of model `j`). `skip[i]` keeps a
model out of the LRT pairs (a non-converged fit) but in the AIC table. Both tables
are sorted by log-likelihood, best first; select on `argmin(aic)`, not row 1.
"""
function _aic_lrt_tables(names::Vector{String}, n_params::Vector{Int},
                         logliks::Vector{Float64}, is_nested;
                         extra_cols = (),
                         skip::AbstractVector{Bool} = fill(false, length(names)))
    aics = @. 2 * n_params - 2 * logliks
    daic = aics .- minimum(aics)
    wraw = exp.(-0.5 .* daic)
    ord = sortperm(aics)

    aic_df = DataFrame(model = names[ord], n_params = n_params[ord], loglik = logliks[ord],
                       aic = aics[ord], delta_aic = daic[ord], weight = (wraw ./ sum(wraw))[ord])
    for (nm, vals) in extra_cols
        aic_df[!, nm] = vals[ord]
    end

    rows = NamedTuple{(:model_null, :model_alt, :loglik_null, :loglik_alt, :df, :D, :p_value),
                      Tuple{String,String,Float64,Float64,Int,Float64,Float64}}[]
    for i in eachindex(names), j in eachindex(names)
        i == j && continue
        (skip[i] || skip[j]) && continue
        is_nested(i, j) || continue
        df = n_params[j] - n_params[i]
        df < 1 && continue
        D = 2.0 * (logliks[j] - logliks[i])
        D < 0.0 && continue          # the richer model must not fit worse
        push!(rows, (model_null = names[i], model_alt = names[j], loglik_null = logliks[i],
                     loglik_alt = logliks[j], df = df, D = D, p_value = ccdf(Chisq(df), D)))
    end
    lrt_df = isempty(rows) ?
        DataFrame(model_null = String[], model_alt = String[], loglik_null = Float64[],
                  loglik_alt = Float64[], df = Int[], D = Float64[], p_value = Float64[]) :
        DataFrame(rows)

    return (aic = sort!(aic_df, :loglik, rev = true),
            lrt = sort!(lrt_df, [:loglik_null, :loglik_alt], rev = true))
end

_free_set(s::AbstractString) = Set(split(s, ","; keepempty = false))

"""
    dfe_model_comparison(fits::DataFrame) -> (aic = DataFrame, lrt = DataFrame)
    dfe_model_comparison(sfs, divergence; models, cache = DFECache(), kwargs...)
        -> (aic = DataFrame, lrt = DataFrame, fits = DataFrame)

AIC (with Akaike weights) and likelihood-ratio tests between DFE models fitted to
the same data, per dataset. `fits` is `dfe_ml` output for several models stacked
(`vcat(...; cols = :union)`). Model `i` is nested in `j` when its `free_params` are
a strict subset of `j`'s (GammaZero in GammaExpo, ExpoZero in GammaZero,
GammaZero in `WithLethal(GammaZero())`), and non-converged fits are left out of the
LRTs. Both tables carry a `dataset` column and are sorted by log-likelihood within
dataset; pick the AIC-best model with `argmin(aic)`.

The second form fits `models` (default: the eleven continuous shapes) with `dfe_ml`,
sharing one `DFECache`, and forwards `kwargs` to it.
"""
function dfe_model_comparison(fits::DataFrame)
    for c in (:dataset, :model, :loglik, :n_params, :free_params, :converged)
        hasproperty(fits, c) || throw(ArgumentError("dfe_model_comparison: column $c missing; pass dfe_ml output"))
    end
    aics, lrts = DataFrame[], DataFrame[]
    for ds in unique(fits.dataset)
        sub = fits[fits.dataset .== ds, :]
        names = String.(sub.model)
        allunique(names) || throw(ArgumentError("dataset $ds: model names repeat: $names"))
        free = _free_set.(sub.free_params)
        nested(i, j) = free[i] ⊊ free[j]
        extra = Pair{Symbol,Vector}[:converged => collect(Bool, sub.converged)]
        for c in (:alpha, :omega_a)
            hasproperty(sub, c) && push!(extra, c => collect(sub[!, c]))
        end
        t = _aic_lrt_tables(names, collect(Int, sub.n_params), collect(Float64, sub.loglik), nested;
                            extra_cols = extra, skip = .!collect(Bool, sub.converged))
        insertcols!(t.aic, 1, :dataset => fill(ds, nrow(t.aic)))
        insertcols!(t.lrt, 1, :dataset => fill(ds, nrow(t.lrt)))
        push!(aics, t.aic); push!(lrts, t.lrt)
    end
    return (aic = vcat(aics...), lrt = vcat(lrts...))
end

const _DEFAULT_COMPARISON_MODELS = (ExpoZero(), GammaZero(), ExpoExpo(), GammaExpo(), ExpoGamma(),
                                    GammaGamma(), DisplGamma(), ScaledBeta(), LogNormalZero(),
                                    NeuGamma(), FGMBesselK())

function dfe_model_comparison(sfs, divergence;
                              models = collect(AbstractDFEModel, _DEFAULT_COMPARISON_MODELS),
                              cache::DFECache = DFECache(), kwargs...)
    fits = DataFrame[]
    for m in models
        @info "Fitting $(dfe_name(m))"
        push!(fits, dfe_ml(sfs, divergence, m; cache = cache, kwargs...))
    end
    all_fits = vcat(fits...; cols = :union)
    t = dfe_model_comparison(all_fits)
    return (aic = t.aic, lrt = t.lrt, fits = all_fits)
end

# -----------------------------------------------------------------------------
# Internal: the objective of one dataset, and the log-likelihood at given
# parameters (used by the parity tests and fixtures)
# -----------------------------------------------------------------------------

function _objective_for(model::AbstractDFEModel, sfs::AbstractMatrix, div;
                        integrator::AbstractIntegrator = ExactSojourn(),
                        grid::AbstractSGrid = GrapesGrid(),
                        dominance::AbstractDominance = Additive(),
                        expectation::AbstractSFSExpectation = PRFSojourn(),
                        kernel::AbstractLikelihoodKernel = PoissonKernel(),
                        quadrature::AbstractDFEQuadrature = MidpointDensity(),
                        div_model::AbstractDivergenceModel = SaturatedDivergence(),
                        demography::AbstractDemography = NuisanceRj(),
                        ne_ratio::Float64 = 1.0, fold::Bool = false,
                        no_syn_orient_error::Bool = false,
                        fixed::NamedTuple = NamedTuple(), cache::DFECache = DFECache())
    n = Int(round(maximum(sfs[:, 1]))) + 1
    pre = get_precomputed(cache, n; integrator = integrator, grid = grid,
                          dominance = dominance, expectation = expectation)
    data = build_data_input(sfs, first(_div_list(div)), fold, uses_divergence(div_model);
                            neutral_shape = pre.neutral_shape)
    lay = build_layout(model, kernel, (!fold) && (!no_syn_orient_error); fixed = fixed)
    return LikelihoodObjective(data, pre, model, lay; constrained_div = !divergence_free(div_model),
                               no_syn_orient_error = no_syn_orient_error,
                               free_rj = free_rj(demography), ne_ratio = ne_ratio,
                               kernel = kernel, quad = quadrature)
end

# `theta`: natural-scale values of every parameter of the layout, in order (DFE
# parameters, then syno_orientation_error, then the kernel's).
function _loglik_at(model::AbstractDFEModel, sfs::AbstractMatrix, div, theta::AbstractVector{<:Real};
                    kwargs...)
    obj = _objective_for(model, sfs, div; kwargs...)
    lay = obj.layout
    length(theta) == length(lay.names) ||
        throw(ArgumentError("_loglik_at: $(length(theta)) values for the parameters $(lay.names)"))
    for (i, idx) in enumerate(lay.fixed_idx)
        obj.full_param[idx] = lay.fixed_vals[i]
    end
    x = [_to_opt(Float64(theta[i]), lay.log_scale[i]) for i in lay.opt_idx]
    return -obj(x)
end

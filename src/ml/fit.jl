# =============================================================================
# Fitting: parameter layout, starting points, optimizers, the fit, Route A, and
# the reported summaries
# =============================================================================

"""
    build_layout(model, kernel, include_oer; fixed = (;), bounds = (;)) -> ParamLayout

The parameter vector of a fit: the DFE parameters (`dfe_parameters` order), then
`syno_orientation_error` (when the SFS is unfolded and it is not switched off), then
the kernel's parameter. `bounds` replaces bounds and `fixed` pins values, both by
name; an unknown name is an error.
"""
function build_layout(model::AbstractDFEModel, kernel::AbstractLikelihoodKernel, include_oer::Bool;
                      fixed::NamedTuple = NamedTuple(), bounds::NamedTuple = NamedTuple())
    dp = dfe_parameters(model)
    names = collect(Symbol, keys(dp))
    bnds = Tuple{Float64,Float64}[(Float64(b[1]), Float64(b[2])) for b in values(dp)]
    ls = collect(Bool, dfe_log_scale(model))
    length(ls) == length(names) ||
        throw(ArgumentError("$(dfe_name(model)): dfe_log_scale has $(length(ls)) entries for $(length(names)) parameters"))
    n_spec = length(names)

    oer_idx = 0
    if include_oer
        push!(names, :syno_orientation_error); push!(bnds, (0.0, 0.1)); push!(ls, false)
        oer_idx = length(names)
    end

    kp = kernel_parameters(kernel)
    length(kp) <= 1 || throw(ArgumentError("a likelihood kernel may have at most one parameter"))
    kernel_idx = 0
    for (k, b) in pairs(kp)
        push!(names, k); push!(bnds, (Float64(b[1]), Float64(b[2]))); push!(ls, b[1] > 0)
        kernel_idx = length(names)
    end
    allunique(names) || throw(ArgumentError("duplicate parameter names: $names"))

    _check_names(nt, what) = for k in keys(nt)
        k in names || throw(ArgumentError("$what: $k is not a parameter of this fit. " *
                                          "Parameters: $(join(names, ", "))"))
    end
    _check_names(bounds, "bounds")
    _check_names(fixed, "fixed")
    for (k, b) in pairs(bounds)
        lo, hi = Float64(b[1]), Float64(b[2])
        lo < hi || throw(ArgumentError("bounds for $k: need lo < hi, got ($lo, $hi)"))
        bnds[findfirst(==(k), names)] = (lo, hi)
    end

    fixed_idx = Int[findfirst(==(k), names) for k in keys(fixed)]
    fixed_vals = Float64[Float64(v) for v in values(fixed)]
    opt_idx = [i for i in eachindex(names) if !(i in fixed_idx)]
    return ParamLayout(names, bnds, ls, opt_idx, fixed_idx, fixed_vals, n_spec, oer_idx, kernel_idx)
end

_to_opt(x::Float64, logscale::Bool) = logscale ? log(max(x, 1e-300)) : x
_from_opt(y::Float64, logscale::Bool) = logscale ? exp(y) : y

function _opt_bounds(lay::ParamLayout)
    return [(lay.log_scale[i] ? log(max(lay.bounds[i][1], 1e-300)) : lay.bounds[i][1],
             lay.log_scale[i] ? log(lay.bounds[i][2]) : lay.bounds[i][2]) for i in lay.opt_idx]
end

# The model's starts (with syno_orientation_error = 0.01 and the kernel's start),
# clamped to the bounds, then the centre of the box in optimizer space, then
# `nb_rand_start` uniform draws in optimizer space (MersenneTwister(42), so a fit is
# reproducible).
function starting_points(model::AbstractDFEModel, lay::ParamLayout, kernel::AbstractLikelihoodKernel,
                         nb_rand_start::Int)
    ks = kernel_start(kernel)
    starts = Vector{Vector{Float64}}()
    for s in dfe_start_points(model)
        length(s) == lay.n_spec ||
            throw(ArgumentError("$(dfe_name(model)): a start point has $(length(s)) values for $(lay.n_spec) parameters"))
        full = zeros(length(lay.names))
        full[1:lay.n_spec] .= s
        lay.oer_idx != 0 && (full[lay.oer_idx] = 0.01)
        lay.kernel_idx != 0 && (full[lay.kernel_idx] = ks[1])
        for i in eachindex(full)
            full[i] = clamp(full[i], lay.bounds[i][1], lay.bounds[i][2])
        end
        push!(starts, [_to_opt(full[i], lay.log_scale[i]) for i in lay.opt_idx])
    end
    ob = _opt_bounds(lay)
    push!(starts, [(lo + hi) / 2.0 for (lo, hi) in ob])
    if nb_rand_start > 0
        rng = MersenneTwister(42)
        for _ in 1:nb_rand_start
            push!(starts, [rand(rng) * (hi - lo) + lo for (lo, hi) in ob])
        end
    end
    return starts
end

function _nudge_to_interior!(x::Vector{Float64}, bounds::Vector{Tuple{Float64,Float64}})
    for i in eachindex(x)
        lo, hi = bounds[i]
        if hi <= lo
            x[i] = lo
            continue
        end
        eps = (hi - lo) * 1e-8
        if x[i] <= lo
            x[i] = lo + eps
        elseif x[i] >= hi
            x[i] = hi - eps
        end
    end
    return x
end

function run_optimizer(opt::LBFGSOptimizer, f, bounds::Vector{Tuple{Float64,Float64}}, x0::Vector{Float64})
    result = Optim.optimize(f, first.(bounds), last.(bounds), x0,
                            Optim.Fminbox(Optim.LBFGS(linesearch = LineSearches.BackTracking())),
                            Optim.Options(iterations = opt.max_iter, f_reltol = opt.f_reltol);
                            autodiff = ADTypes.AutoFiniteDiff())
    return Optim.minimizer(result), Optim.minimum(result), Optim.converged(result)
end

function run_optimizer(opt::NelderMeadOptimizer, f, bounds::Vector{Tuple{Float64,Float64}}, x0::Vector{Float64})
    lo = first.(bounds)
    hi = last.(bounds)
    to_q(p) = @. log((p - lo) / (hi - p))
    to_p(q) = @. lo + (hi - lo) / (1 + exp(-q))
    r = Optim.optimize(q -> f(to_p(q)), to_q(clamp.(x0, lo .+ 1e-7 .* (hi .- lo), hi .- 1e-7 .* (hi .- lo))),
                       Optim.NelderMead(), Optim.Options(iterations = opt.iterations))
    return to_p(Optim.minimizer(r)), Optim.minimum(r), Optim.converged(r)
end

# -----------------------------------------------------------------------------
# One fit
# -----------------------------------------------------------------------------

struct _FitOutcome{O<:LikelihoodObjective}
    obj::O                 # left at the best point: its buffers hold the fitted SFS
    x::Vector{Float64}
    loglik::Float64
    theta::Float64
    converged::Bool
    n_opt::Int
end

function _fit(data::DataInput, pre::PrecomputedData, model::AbstractDFEModel;
              kernel, quad, optimizer, ne_ratio::Float64, constrained_div::Bool,
              free_rj::Bool, no_syn_orient_error::Bool, fixed::NamedTuple, bounds::NamedTuple,
              nb_rand_start::Int)
    include_oer = (!data.folded) && (!no_syn_orient_error)
    lay = build_layout(model, kernel, include_oer; fixed = fixed, bounds = bounds)
    obj = LikelihoodObjective(data, pre, model, lay; constrained_div = constrained_div,
                              no_syn_orient_error = no_syn_orient_error, free_rj = free_rj,
                              ne_ratio = ne_ratio, kernel = kernel, quad = quad)
    if isempty(lay.opt_idx)
        x = Float64[]
        converged = true
    else
        ob = _opt_bounds(lay)
        best_val, x, converged = Inf, Float64[], false
        for x0 in starting_points(model, lay, kernel, nb_rand_start)
            _nudge_to_interior!(x0, ob)
            xm, val, conv = run_optimizer(optimizer, obj, ob, x0)
            if val < best_val
                best_val, x, converged = val, collect(Float64, xm), conv
            end
        end
    end
    _set_params!(obj, x)
    ll, theta = _loglik!(obj)
    return _FitOutcome(obj, x, ll, theta, converged, length(lay.opt_idx))
end

# -----------------------------------------------------------------------------
# What a finished fit exposes to the summaries
# -----------------------------------------------------------------------------

"""
    FitContext

What a finished fit hands to `summarize`:

- `data::DataInput`, `pre::PrecomputedData`, `model`, `params` (DFE parameters at the
  optimum, natural scale), `full_params` (with the nuisance parameters), `names`;
- `theta`, `loglik`, `exp_neut`, `exp_sel` (fitted expected SFS), `kernel`, `disp`;
- the divergence decomposition at the ancestral Ne: `totN` (the denominator: the
  observed Dn under `SaturatedDivergence`, the DFE prediction otherwise), `pred_totN`
  (always the DFE prediction), `exp_div_N1`, `exp_div_N2`, `exp_div_N3` (expected Dn
  from mutations with S <= -nearly_neutral, <= 0, <= nearly_neutral),
  `effective_Ln_div`, `lambda` (neutral substitution rate per site);
- the settings `nearly_neutral`, `alpha_ci_threshold`,
  `divergence_free`.
"""
struct FitContext{M<:AbstractDFEModel,P<:PrecomputedData,K<:AbstractLikelihoodKernel}
    data::DataInput
    pre::P
    model::M
    params::Vector{Float64}
    full_params::Vector{Float64}
    names::Vector{Symbol}
    theta::Float64
    loglik::Float64
    exp_neut::Vector{Float64}
    exp_sel::Vector{Float64}
    kernel::K
    disp::Float64
    totN::Float64
    pred_totN::Float64
    exp_div_N1::Float64
    exp_div_N2::Float64
    exp_div_N3::Float64
    effective_Ln_div::Float64
    lambda::Float64
    nearly_neutral::Float64
    alpha_ci_threshold::Float64
    divergence_free::Bool
end

function FitContext(f::_FitOutcome; nearly_neutral::Float64,
                    alpha_ci_threshold::Float64)
    obj = f.obj
    data, pre = obj.data, obj.pre
    theta = f.theta
    pdiv = obj.param_spec .* obj.ne_mult
    dvals = similar(pre.s_grid)
    fill_dfe_vals!(dvals, obj.quad, obj.model, pdiv, pre, obj.scratch)
    theta_div = theta * obj.ne_ratio
    # With no divergence data any length cancels in alpha; GRAPES' ARBITRARY_LENGTH.
    effL = data.div_available ? data.Ln_div : 100_000.0
    pred(lam, minS) = expected_divergence_from_dfe(dvals, obj.model, pdiv, pre, effL, lam,
                                                   theta_div, data.n; minS = minS)
    opt_lambda = 1.0
    if !data.div_available
        totN = pred(1.0, Inf)
    elseif !obj.constrained_div
        totN = data.fixN
        opt_lambda = (data.fixS / data.Ls_div) - theta / data.n
    else
        prov_lambda = (data.fixS / data.Ls_div) - theta / data.n
        prov_totN = pred(prov_lambda, Inf)
        coeff = (data.fixS + data.fixN) / (prov_lambda * data.Ls_div + prov_totN)
        opt_lambda = prov_lambda * coeff
        totN = prov_totN * coeff
    end
    nn = Float64(nearly_neutral)
    return FitContext(data, pre, obj.model, copy(obj.param_spec), copy(obj.full_param),
                      copy(obj.layout.names), theta, f.loglik, copy(obj.exp_neut), copy(obj.exp_sel),
                      obj.kernel, _kernel_disp(obj), totN, pred(opt_lambda, Inf),
                      pred(opt_lambda, -nn), pred(opt_lambda, 0.0), pred(opt_lambda, nn),
                      effL, opt_lambda, nn, alpha_ci_threshold, !obj.constrained_div)
end

function summarize(::DFEMAlpha, c::FitContext)
    d = c.data
    alpha = (c.totN - c.exp_div_N3) / c.totN
    omega_a = (alpha * c.totN * d.Ls_div) / (d.fixS * c.effective_Ln_div)
    omega_na = ((1.0 - alpha) * c.totN * d.Ls_div) / (d.fixS * c.effective_Ln_div)
    omega_obs = (d.div_available && d.Ls_div > 0 && d.Ln_div > 0) ?
                (d.fixN / d.Ln_div) / (d.fixS / d.Ls_div) : NaN
    alpha_div = (d.div_available && d.fixN > 0) ? (d.fixN - c.exp_div_N3) / d.fixN : NaN
    alpha_dfe = c.pred_totN > 0 ? (c.pred_totN - c.exp_div_N3) / c.pred_totN : NaN

    a_lo = a_hi = oa_lo = oa_hi = ona_lo = ona_hi = NaN
    if d.div_available && c.divergence_free && c.alpha_ci_threshold > 0 && d.fixN > 0 && d.fixS > 0
        opt_lnL = loglik_from_expectations(d, c.exp_neut, c.exp_sel, d.fixS, d.fixN, c.kernel, c.disp)
        a_lo, a_hi = alpha_ci_bounds(d, c.exp_neut, c.exp_sel, d.fixS, d.fixN, c.exp_div_N3,
                                     opt_lnL, c.alpha_ci_threshold, c.kernel, c.disp)
        if alpha != 0.0
            oa_lo, oa_hi = omega_a * (a_lo / alpha), omega_a * (a_hi / alpha)
        end
        ona_lo = (omega_a + omega_na) * (1.0 - a_hi)
        ona_hi = (omega_a + omega_na) * (1.0 - a_lo)
    end
    return (alpha = alpha, alpha_ci_low = a_lo, alpha_ci_high = a_hi,
            omega_a = omega_a, omega_a_ci_low = oa_lo, omega_a_ci_high = oa_hi,
            omega_na = omega_na, omega_na_ci_low = ona_lo, omega_na_ci_high = ona_hi,
            omega_obs = omega_obs, alpha_div = alpha_div, alpha_dfe = alpha_dfe)
end

# Warn for every free DFE parameter whose estimate sits on a bound: a clamped value
# is not an estimate. (Nuisance parameters are skipped: syno_orientation_error = 0 is
# a common, legitimate optimum.)
function _warn_pinned_bounds(lay::ParamLayout, full::Vector{Float64}, model::AbstractDFEModel)
    hits = String[]
    for i in lay.opt_idx
        i <= lay.n_spec || continue
        v = full[i]
        lo, hi = lay.bounds[i]
        if isapprox(v, hi; rtol = 1e-4)
            push!(hits, "$(lay.names[i]) = $(round(v, sigdigits = 6)) at its UPPER bound $hi")
        elseif isapprox(v, lo; rtol = 1e-4, atol = lo == 0 ? 1e-10 : 0.0)
            push!(hits, "$(lay.names[i]) = $(round(v, sigdigits = 6)) at its LOWER bound $lo")
        end
    end
    isempty(hits) && return
    @warn "$(dfe_name(model)): $(length(hits)) parameter(s) ON a bound; a clamped parameter is " *
          "not an estimate:\n  " * join(hits, "\n  ") *
          "\nWiden with `bounds = (name = (lo, hi),)` if the bound is not physical."
end

# -----------------------------------------------------------------------------
# Route A: the Ne ratio from the divergence counts, given the DFE from the SFS
# -----------------------------------------------------------------------------

# 1-D Brent search of log(r) in [0, log 50] on the likelihood of (Dn, Ds) given the
# stage-1 DFE. Returns 1 at the lower bound when that likelihood only decreases.
function _estimate_ne_ratio_1d(f::_FitOutcome)
    obj = f.obj
    data, pre, model = obj.data, obj.pre, obj.model
    theta = f.theta
    disp = _kernel_disp(obj)
    lambda = (data.fixS / data.Ls_div) - theta / data.n
    lambda < 1e-12 && (lambda = 1e-12)
    p0 = copy(obj.param_spec)
    pw = similar(p0)
    dvals = similar(pre.s_grid)
    kernel = obj.kernel

    function neg_div_ll(log_ratio::Float64)
        ratio = exp(log_ratio)
        pw .= p0 .* dfe_ne_multipliers(model, ratio)
        fill_dfe_vals!(dvals, obj.quad, model, pw, pre, obj.scratch)
        exp_divN = expected_divergence_from_dfe(dvals, model, pw, pre, data.Ln_div, lambda,
                                                theta * ratio, data.n)
        exp_divS = lambda * data.Ls_div
        tot = exp_divN + exp_divS
        tot <= 0.0 && return 1e10
        coeff = (data.fixN + data.fixS) / tot
        en = exp_divN * coeff
        es = exp_divS * coeff
        (en <= 0.0 || es <= 0.0) && return 1e10
        kernel isa PoissonKernel && return -(data.fixN * log(en) - en + data.fixS * log(es) - es)
        return -(log_kernel(kernel, data.fixN, en, 0.0, disp) + log_kernel(kernel, data.fixS, es, 0.0, disp))
    end

    res = Optim.optimize(neg_div_ll, 0.0, log(50.0), Optim.Brent())
    return exp(Optim.minimizer(res))
end

# -----------------------------------------------------------------------------
# The fit, as one output row
# -----------------------------------------------------------------------------

function fit_model(data::DataInput, pre::PrecomputedData, model::AbstractDFEModel;
                   kernel::AbstractLikelihoodKernel, quad::AbstractDFEQuadrature,
                   optimizer::AbstractOptimizer, ne_model::AbstractNeModel,
                   div_model::AbstractDivergenceModel, demography::AbstractDemography,
                   no_syn_orient_error::Bool, fixed::NamedTuple, bounds::NamedTuple,
                   nb_rand_start::Int, nearly_neutral::Float64,
                   alpha_ci_threshold::Float64, summaries)
    r = ne_ratio(ne_model)
    common = (kernel = kernel, quad = quad, optimizer = optimizer,
              free_rj = free_rj(demography), no_syn_orient_error = no_syn_orient_error,
              fixed = fixed, bounds = bounds, nb_rand_start = nb_rand_start)
    route_a = r === :free
    if route_a
        divergence_free(div_model) || throw(ArgumentError(
            "AncestralNeRatio(:free) estimates the ratio from Dn and Ds, so Dn cannot also be " *
            "predicted by the DFE; use div_model = SaturatedDivergence()"))
        s1 = _fit(data, pre, model; ne_ratio = 1.0, constrained_div = false, common...)
        r_ml = (data.div_available && data.fixN > 0 && data.fixS > 0) ? _estimate_ne_ratio_1d(s1) : 1.0
        fit = _fit(data, pre, model; ne_ratio = r_ml, constrained_div = false, common...)
    else
        r_ml = Float64(r)
        fit = _fit(data, pre, model; ne_ratio = r_ml, constrained_div = !divergence_free(div_model),
                   common...)
    end

    lay = fit.obj.layout
    _warn_pinned_bounds(lay, fit.obj.full_param, model)
    ctx = FitContext(fit; nearly_neutral = nearly_neutral,
                     alpha_ci_threshold = alpha_ci_threshold)

    n_params = fit.n_opt + (route_a ? 1 : 0)
    free = [string(lay.names[i]) for i in lay.opt_idx]
    route_a && push!(free, "Ne_ratio")

    row = OrderedDict{Symbol,Any}(
        :model => dfe_name(model), :n => data.n, :converged => fit.converged,
        :loglik => fit.loglik, :aic => 2.0 * n_params - 2.0 * fit.loglik,
        :n_params => n_params, :free_params => join(free, ","), :theta => fit.theta)
    for (i, nm) in enumerate(lay.names)
        _put!(row, nm, fit.obj.full_param[i])
    end
    route_a && _put!(row, :Ne_ratio, r_ml)
    if fit.obj.exterior
        # Mass beyond the outer cell edge, all assigned to that row; reported so a fit
        # resting on it is visible.
        _put!(row, :lethal_mass, dfe_mass_above(model, ctx.params, pre.s_ext_hi))
        _put!(row, :neutral_mass, pre.s_ext_lo > 0 ? dfe_mass_below(model, ctx.params, pre.s_ext_lo) : 0.0)
    end
    for s in summaries
        for (k, v) in pairs(summarize(s, ctx))
            _put!(row, k, v)
        end
    end
    return row
end

function _put!(row::OrderedDict{Symbol,Any}, k::Symbol, v)
    haskey(row, k) && throw(ArgumentError("output column $k is produced twice (a parameter or summary name clashes)"))
    row[k] = v
    return row
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

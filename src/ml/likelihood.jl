# =============================================================================
# Data, the likelihood objective, and the divergence integrals
# =============================================================================

"""
    DataInput

One dataset as the likelihood sees it: the (optionally folded) SFS counts,
divergence counts and site numbers, and the neutral expectation at theta = 1.
"""
struct DataInput
    n::Int
    n_bins::Int
    folded::Bool
    specN::Vector{Float64}
    specS::Vector{Float64}
    specN_int::Vector{Int}
    specS_int::Vector{Int}
    logfactN::Vector{Float64}
    logfactS::Vector{Float64}
    Ln_poly::Float64
    Ls_poly::Float64
    Ln_div::Float64
    Ls_div::Float64
    fixN::Float64
    fixS::Float64
    div_available::Bool
    neutral_base::Vector{Float64}
end

"""
    fold_sfs!(out, unfolded) -> out

Fold an unfolded SFS of `n - 1` bins into `n / 2` bins (`out[j] = x[j] + x[n-j]`,
the middle bin once).
"""
function fold_sfs!(out::AbstractVector{Float64}, unfolded::AbstractVector{Float64})
    n = length(unfolded) + 1
    folded_len = n ÷ 2
    length(out) == folded_len ||
        throw(DimensionMismatch("fold_sfs!: out has length $(length(out)), expected $folded_len"))
    @inbounds for j in 1:folded_len
        comp = n - j
        out[j] = j == comp ? unfolded[j] : unfolded[j] + unfolded[comp]
    end
    return out
end

"""
    compute_neutral_base(Ls, n, folded, shape) -> Vector{Float64}

Expected neutral SFS at theta = 1 and `r_j = 1`: `(Ls/j) * shape[j]`, folded after
the shape is applied.
"""
function compute_neutral_base(Ls::Float64, n::Int, folded::Bool, shape::AbstractVector{Float64})
    length(shape) == n - 1 ||
        error("neutral_shape has length $(length(shape)), expected $(n - 1) for n=$n")
    if folded
        n_bins = n ÷ 2
        base = Vector{Float64}(undef, n_bins)
        for j in 1:n_bins
            base[j] = Ls / j * shape[j]
            if (n % 2 != 0) || (j != n ÷ 2)
                base[j] += Ls / (n - j) * shape[n - j]
            end
        end
        return base
    else
        return [Ls / j * shape[j] for j in 1:(n - 1)]
    end
end

"""
    build_data_input(sfs, div, fold, use_div; Ln_poly, Ls_poly, neutral_shape) -> DataInput

`sfs` has columns `[DAC, Pn, Ps]` (unfolded, `n - 1` rows, so `n = max(DAC) + 1`);
`div` is `[Dn, Ds, Ln, Ls]`. `Ln_poly`/`Ls_poly` default to the divergence site
counts.
"""
function build_data_input(sfs_mat::AbstractMatrix, div::AbstractVector, fold::Bool, use_div::Bool;
                          Ln_poly::Float64 = NaN, Ls_poly::Float64 = NaN,
                          neutral_shape::Union{Nothing,AbstractVector{Float64}} = nothing)
    size(sfs_mat, 2) >= 3 ||
        throw(ArgumentError("the SFS matrix needs the columns [DAC, Pn, Ps]"))
    length(div) == 4 ||
        throw(ArgumentError("divergence must be [Dn, Ds, Ln, Ls], got length $(length(div))"))

    n = Int(round(maximum(sfs_mat[:, 1]))) + 1
    size(sfs_mat, 1) == n - 1 ||
        throw(ArgumentError("the SFS has $(size(sfs_mat, 1)) rows but its largest DAC is $(n - 1); " *
                            "dfe_ml needs every bin 1..n-1 (do not trim the SFS)"))
    specN = Vector{Float64}(sfs_mat[:, 2])
    specS = Vector{Float64}(sfs_mat[:, 3])

    fixN, fixS, Ln_div, Ls_div = Float64(div[1]), Float64(div[2]), Float64(div[3]), Float64(div[4])
    div_available = use_div && (fixN > 0 || fixS > 0)
    if use_div && fixN <= 0 && fixS <= 0
        @warn "Divergence counts are zero (Dn=$(fixN), Ds=$(fixS)). Running with polymorphism data only."
    end

    Ln_poly = isnan(Ln_poly) ? Ln_div : Ln_poly
    Ls_poly = isnan(Ls_poly) ? Ls_div : Ls_poly

    if fold
        specN_f = Vector{Float64}(undef, n ÷ 2)
        specS_f = Vector{Float64}(undef, n ÷ 2)
        fold_sfs!(specN_f, specN)
        fold_sfs!(specS_f, specS)
        specN, specS = specN_f, specS_f
    end

    specN_int = [Int(round(x)) for x in specN]
    specS_int = [Int(round(x)) for x in specS]
    logfactN = [loggamma(x + 1) for x in specN_int]
    logfactS = [loggamma(x + 1) for x in specS_int]

    sh = neutral_shape === nothing ? ones(n - 1) : neutral_shape
    neutral_base = compute_neutral_base(Ls_poly, n, fold, sh)

    return DataInput(n, length(specN), fold, specN, specS, specN_int, specS_int,
                     logfactN, logfactS, Ln_poly, Ls_poly, Ln_div, Ls_div, fixN, fixS,
                     div_available, neutral_base)
end

# -----------------------------------------------------------------------------
# Parameter layout: the DFE parameters, then misorient, then the
# kernel's parameter. Fixed parameters keep their slot and are not optimised.
# -----------------------------------------------------------------------------

struct ParamLayout
    names::Vector{Symbol}
    bounds::Vector{Tuple{Float64,Float64}}
    log_scale::Vector{Bool}
    opt_idx::Vector{Int}
    fixed_idx::Vector{Int}
    fixed_vals::Vector{Float64}
    n_spec::Int
    oer_idx::Int
    kernel_idx::Int
end

# -----------------------------------------------------------------------------
# The objective
# -----------------------------------------------------------------------------

"""
    LikelihoodObjective

The negative log-likelihood as a callable of the optimizer-space vector, with every
buffer preallocated so an evaluation allocates nothing. Type parameters are the
model, kernel, quadrature and `PrecomputedData` types, so the hot path dispatches
statically. Internal.
"""
struct LikelihoodObjective{M<:AbstractDFEModel,K<:AbstractLikelihoodKernel,
                           Q<:AbstractDFEQuadrature,P<:PrecomputedData}
    data::DataInput
    pre::P
    model::M
    layout::ParamLayout
    constrained_div::Bool        # Dn predicted by the DFE (DFEConstrainedDivergence)
    no_syn_orient_error::Bool
    free_rj::Bool
    ne_ratio::Float64
    ne_mult::Vector{Float64}     # per DFE parameter, see `dfe_ne_multipliers`
    exterior::Bool
    kernel::K
    quad::Q
    dfe_vals::Vector{Float64}
    param_spec::Vector{Float64}
    full_param::Vector{Float64}
    exp_sel_unfolded::Vector{Float64}
    exp_sel::Vector{Float64}
    exp_neut::Vector{Float64}
    r::Vector{Float64}
    w::Vector{Float64}
    dfe_vals_div::Vector{Float64}
    param_spec_div::Vector{Float64}
    neutral_row::Vector{Float64}  # S = 0 expected counts, unfolded
    del_row::Vector{Float64}      # most deleterious row, unfolded
    scratch::Vector{Float64}      # quadrature scratch, length(pre.s_edges)
end

function LikelihoodObjective(data::DataInput, pre::PrecomputedData, model::AbstractDFEModel,
                             layout::ParamLayout; constrained_div::Bool, no_syn_orient_error::Bool,
                             free_rj::Bool, ne_ratio::Float64,
                             kernel::AbstractLikelihoodKernel, quad::AbstractDFEQuadrature)
    data.n == pre.n || throw(ArgumentError("data has n = $(data.n) but the tables have n = $(pre.n)"))
    n_s = length(pre.s_grid)
    full_param = zeros(length(layout.names))
    @inbounds for (i, idx) in enumerate(layout.fixed_idx)
        full_param[idx] = layout.fixed_vals[i]
    end
    n_unf = data.n - 1
    return LikelihoodObjective(data, pre, model, layout, constrained_div, no_syn_orient_error,
        free_rj, ne_ratio, dfe_ne_multipliers(model, ne_ratio), uses_exterior(pre.expectation),
        kernel, quad,
        zeros(n_s), zeros(layout.n_spec), full_param, zeros(n_unf), zeros(data.n_bins),
        zeros(data.n_bins), ones(data.n_bins), zeros(n_s), zeros(n_s), zeros(layout.n_spec),
        [pre.neutral_shape[j] / j for j in 1:n_unf],
        [pre.expected_counts[1, j] for j in 1:n_unf],
        zeros(length(pre.s_edges)))
end

# Optimizer space -> natural parameters, into `full_param` and its DFE slice.
@inline function _set_params!(obj::LikelihoodObjective, x_opt::AbstractVector{Float64})
    lay = obj.layout
    @inbounds for k in eachindex(lay.opt_idx)
        idx = lay.opt_idx[k]
        obj.full_param[idx] = lay.log_scale[idx] ? exp(x_opt[k]) : x_opt[k]
    end
    @inbounds for i in 1:lay.n_spec
        obj.param_spec[i] = obj.full_param[i]
    end
    return obj
end

# Expected synonymous and nonsynonymous SFS at the current parameters, into
# `exp_neut` / `exp_sel`, with theta and r_j profiled from Ps + Pn (GRAPES'
# `expected_SFS_div`). Returns theta. `dfe_vals` holds the DFE weights afterwards.
function _expected_sfs!(obj::LikelihoodObjective)
    data, pre = obj.data, obj.pre
    fill_dfe_vals!(obj.dfe_vals, obj.quad, obj.model, obj.param_spec, pre, obj.scratch)

    @inbounds for i in eachindex(obj.dfe_vals)
        obj.w[i] = obj.dfe_vals[i] * pre.s_widths[i]
    end
    # sum_i w_i * counts[i, j], column by column (contiguous in i)
    C = pre.expected_counts
    @inbounds for j in eachindex(obj.exp_sel_unfolded)
        acc = 0.0
        @simd for i in eachindex(obj.w)
            acc += obj.w[i] * C[i, j]
        end
        obj.exp_sel_unfolded[j] = acc
    end

    # Point masses, on unfolded bins before the Ln scaling.
    for (S_k, w_k) in discrete_components(obj.model, obj.param_spec)
        row = discrete_row(pre, Float64(S_k))
        @inbounds for j in eachindex(obj.exp_sel_unfolded)
            obj.exp_sel_unfolded[j] += w_k * row[j]
        end
    end

    if obj.exterior
        w_neu = pre.s_ext_lo > 0.0 ? dfe_mass_below(obj.model, obj.param_spec, pre.s_ext_lo) : 0.0
        w_del = dfe_mass_above(obj.model, obj.param_spec, pre.s_ext_hi)
        @inbounds for j in eachindex(obj.exp_sel_unfolded)
            obj.exp_sel_unfolded[j] += w_neu * obj.neutral_row[j] + w_del * obj.del_row[j]
        end
    end

    @inbounds for j in eachindex(obj.exp_sel_unfolded)
        obj.exp_sel_unfolded[j] *= data.Ln_poly
    end

    if data.folded
        fold_sfs!(obj.exp_sel, obj.exp_sel_unfolded)
    else
        @inbounds for j in 1:data.n_bins
            obj.exp_sel[j] = obj.exp_sel_unfolded[j]
        end
    end

    lay = obj.layout
    if (!data.folded) && (!obj.no_syn_orient_error) && (lay.oer_idx != 0)
        oer = obj.full_param[lay.oer_idx]
        @inbounds for j in 1:data.n_bins
            comp = data.n - j
            obj.exp_neut[j] = data.neutral_base[j] * (1.0 - oer) + data.neutral_base[comp] * oer
        end
    else
        @inbounds for j in 1:data.n_bins
            obj.exp_neut[j] = data.neutral_base[j]
        end
    end

    denom0 = obj.exp_neut[1] + obj.exp_sel[1]
    theta = denom0 > 0 ? (data.specS[1] + data.specN[1]) / denom0 : 1.0
    obj.r[1] = 1.0
    if obj.free_rj
        @inbounds for j in 2:data.n_bins
            denom = obj.exp_neut[j] + obj.exp_sel[j]
            obj.r[j] = (denom > 0 && theta > 0) ? (data.specS[j] + data.specN[j]) / denom / theta : 1.0
        end
    else
        fill!(obj.r, 1.0)
    end

    @inbounds for j in 1:data.n_bins
        obj.exp_neut[j] *= theta * obj.r[j]
        obj.exp_sel[j] *= theta * obj.r[j]
    end
    return theta
end

@inline _kernel_disp(obj::LikelihoodObjective) =
    obj.layout.kernel_idx == 0 ? NaN : obj.full_param[obj.layout.kernel_idx]

# Log-likelihood at the parameters already set by `_set_params!`. Returns (ll, theta).
function _loglik!(obj::LikelihoodObjective)
    data, pre = obj.data, obj.pre
    theta = _expected_sfs!(obj)
    disp = _kernel_disp(obj)

    ll = 0.0
    @inbounds for j in 1:data.n_bins
        mu_s = obj.exp_neut[j]
        if mu_s > 0
            ll += log_kernel(obj.kernel, data.specS_int[j], mu_s, data.logfactS[j], disp)
        end
        mu_n = obj.exp_sel[j]
        if mu_n > 0
            ll += log_kernel(obj.kernel, data.specN_int[j], mu_n, data.logfactN[j], disp)
        end
    end

    if data.div_available
        if obj.constrained_div
            lambda = data.fixS > 0 && data.Ls_div > 0 ? (data.fixS / data.Ls_div) - (theta / data.n) : 1e-12
            lambda < 1e-12 && (lambda = 1e-12)
            exp_divS = lambda * data.Ls_div

            if obj.ne_ratio != 1.0
                # GRAPES `expected_divergence`: the DFE at the ancestral Ne, with the
                # theta * dN correction (the bare fixation integral would vanish for a
                # DFE pushed off the grid).
                @inbounds for i in 1:obj.layout.n_spec
                    obj.param_spec_div[i] = obj.param_spec[i] * obj.ne_mult[i]
                end
                fill_dfe_vals!(obj.dfe_vals_div, obj.quad, obj.model, obj.param_spec_div, pre, obj.scratch)
                exp_divN = expected_divergence_from_dfe(obj.dfe_vals_div, obj.model, obj.param_spec_div,
                                                        pre, data.Ln_div, lambda, theta * obj.ne_ratio, data.n)
            else
                omega_pred = integrate_fixation(obj.dfe_vals, pre, 1e12; adaptive = false)
                for (S_k, w_k) in discrete_components(obj.model, obj.param_spec)
                    omega_pred += w_k * fixation_probability(effective_s(pre.dominance, Float64(S_k)))
                end
                if obj.exterior
                    pre.s_ext_lo > 0.0 &&
                        (omega_pred += dfe_mass_below(obj.model, obj.param_spec, pre.s_ext_lo))
                    omega_pred += dfe_mass_above(obj.model, obj.param_spec, pre.s_ext_hi) * pre.fixprob[1]
                end
                exp_divN = lambda * data.Ln_div * omega_pred
            end

            tot_exp = exp_divS + exp_divN
            if tot_exp > 0
                coeff = (data.fixS + data.fixN) / tot_exp
                exp_divS *= coeff
                exp_divN *= coeff
            end
            exp_divS > 0 && (ll += log_kernel(obj.kernel, data.fixS, exp_divS, loggamma(data.fixS + 1), disp))
            exp_divN > 0 && (ll += log_kernel(obj.kernel, data.fixN, exp_divN, loggamma(data.fixN + 1), disp))
        else
            # Saturated: expectation = observation. Constant in the parameters; kept so
            # the reported log-likelihood is comparable across divergence models.
            data.fixS > 0 && (ll += log_kernel(obj.kernel, data.fixS, max(data.fixS, 1e-300), loggamma(data.fixS + 1), disp))
            data.fixN > 0 && (ll += log_kernel(obj.kernel, data.fixN, max(data.fixN, 1e-300), loggamma(data.fixN + 1), disp))
        end
    end
    return ll, theta
end

function (obj::LikelihoodObjective)(x_opt::AbstractVector{Float64})
    _set_params!(obj, x_opt)
    ll, _ = _loglik!(obj)
    return isfinite(ll) ? -ll : 1e300
end

# The log-likelihood with the expected SFS held fixed and only E[Dn] varied, for the
# alpha CI (a likelihood drop, so it must use the fit's kernel).
function loglik_from_expectations(data::DataInput, exp_neut::Vector{Float64}, exp_sel::Vector{Float64},
                                  exp_divS::Float64, exp_divN::Float64,
                                  kernel::AbstractLikelihoodKernel, disp::Float64)
    ll = 0.0
    @inbounds for j in 1:data.n_bins
        mu_s = exp_neut[j]
        mu_s > 0 && (ll += log_kernel(kernel, data.specS_int[j], mu_s, data.logfactS[j], disp))
        mu_n = exp_sel[j]
        mu_n > 0 && (ll += log_kernel(kernel, data.specN_int[j], mu_n, data.logfactN[j], disp))
    end
    if data.div_available
        exp_divS > 0 && (ll += log_kernel(kernel, data.fixS, exp_divS, loggamma(data.fixS + 1), disp))
        exp_divN > 0 && (ll += log_kernel(kernel, data.fixN, exp_divN, loggamma(data.fixN + 1), disp))
    end
    return ll
end

# GRAPES' alpha CI: move E[Dn] in steps of 0.1% until the log-likelihood drops by
# `threshold`. NaN for a bound not reached within (Dn/5, 3 Dn).
function alpha_ci_bounds(data::DataInput, exp_neut::Vector{Float64}, exp_sel::Vector{Float64},
                         exp_divS::Float64, opt_divN::Float64, exp_div_N3::Float64,
                         opt_lnL::Float64, threshold::Float64,
                         kernel::AbstractLikelihoodKernel, disp::Float64)
    increment = 0.001
    function bound(dir)
        prev_diff = 0.0
        cur_divN = opt_divN
        while true
            cur_divN *= (1.0 + dir * increment)
            cur_diff = opt_lnL - loglik_from_expectations(data, exp_neut, exp_sel, exp_divS, cur_divN, kernel, disp)
            (cur_diff >= threshold && prev_diff <= threshold) && break
            (dir < 0 ? cur_divN < opt_divN / 5.0 : cur_divN > opt_divN * 3.0) && return NaN
            prev_diff = cur_diff
        end
        return (cur_divN - exp_div_N3) / data.fixN
    end
    return bound(-1), bound(1)
end

function integrate_fixation(dfe_vals::Vector{Float64}, pre::PrecomputedData, threshold::Real; adaptive::Bool)
    n_s = length(pre.s_grid)
    idx = _last_index_leq(pre.s_grid, Float64(threshold))
    acc = 0.0
    if adaptive
        start = idx + 1
        start > n_s && return 0.0
        @inbounds for i in start:n_s
            acc += dfe_vals[i] * pre.fixprob_weighted[i]
        end
    else
        idx < 1 && return 0.0
        @inbounds for i in 1:idx
            acc += dfe_vals[i] * pre.fixprob_weighted[i]
        end
    end
    return acc
end

# -----------------------------------------------------------------------------
# Expected divergence given the DFE (GRAPES `expected_divergence`)
# -----------------------------------------------------------------------------

# Fixation rate at `Seff` plus the correction for fixed differences that are
# still segregating, per site.
function _divergence_weight(Seff::Float64, L::Float64, lambda::Float64, theta::Float64,
                            n::Int, pre::PrecomputedData)
    abs(Seff) <= QUASI_ZERO && return L * (lambda + theta / n)
    return L * (lambda * (Seff / (1.0 - exp(-Seff))) + theta * dn_correction(pre, Seff))
end

function _expected_divergence_grid(dfe_vals::Vector{Float64}, pre::PrecomputedData, L::Float64,
                                   lambda::Float64, theta::Float64, n::Int, minS::Float64)
    idx = _last_index_leq(pre.s_grid, minS)
    idx < 1 && return 0.0
    acc = 0.0
    @inbounds for i in 1:idx
        # The adaptive / non-adaptive split is on the real S; the fixation uses S_eff.
        acc += dfe_vals[i] * _divergence_weight(pre.s_eff[i], L, lambda, theta, n, pre) * pre.s_widths[i]
    end
    return acc
end

# The divergence counterpart of the exterior SFS terms: without it every mutation
# less deleterious than the innermost cell would be fitted to the SFS and then
# dropped from E[Dn]. Zero unless `uses_exterior(pre.expectation)`.
function _exterior_divergence(model::AbstractDFEModel, params::AbstractVector{Float64},
                              pre::PrecomputedData, L::Float64, lambda::Float64,
                              theta::Float64, n::Int, minS::Float64)
    uses_exterior(pre.expectation) || return 0.0
    acc = 0.0
    if pre.s_ext_lo > 0.0 && 0.0 <= minS
        acc += dfe_mass_below(model, params, pre.s_ext_lo) * _divergence_weight(0.0, L, lambda, theta, n, pre)
    end
    S1 = pre.s_grid[1]
    if S1 < 0.0 && S1 <= minS
        acc += dfe_mass_above(model, params, pre.s_ext_hi) * _divergence_weight(pre.s_eff[1], L, lambda, theta, n, pre)
    end
    return acc
end

"""
    expected_divergence_from_dfe(dfe_vals, model, params, pre, L, lambda, theta, n; minS = Inf)

Expected nonsynonymous divergence from the mutations with S <= `minS`: the grid, the
discrete components and (when the expectation uses them) the exterior masses.
"""
function expected_divergence_from_dfe(dfe_vals::Vector{Float64}, model::AbstractDFEModel,
                                      params::AbstractVector{Float64}, pre::PrecomputedData,
                                      L::Float64, lambda::Float64, theta::Float64, n::Int;
                                      minS::Float64 = Inf)
    acc = _expected_divergence_grid(dfe_vals, pre, L, lambda, theta, n, minS)
    for (S, w) in discrete_components(model, params)
        S <= minS || continue
        acc += w * _divergence_weight(effective_s(pre.dominance, Float64(S)), L, lambda, theta, n, pre)
    end
    return acc + _exterior_divergence(model, params, pre, L, lambda, theta, n, minS)
end

################################
######     MK methods     ######
################################

# Haploid sample size behind an SFS, the denominator of the frequency axis. It is only
# `size(sfs, 1) + 1` for a full spectrum: `parse_sfs` keeps the true DAC in column 1 when
# a cutoff removes rows, so a trimmed SFS needs `n`. A low cutoff (DACs not starting at 1)
# is detected and refused; a high-only cutoff cannot be detected, so pass `n` or use the
# alpha(x) table returned by `parse_sfs`.
function _alpha_x_n(sfs::AbstractMatrix, n::Union{Nothing,Int})
    if !isnothing(n)
        n > 0 || throw(ArgumentError("n must be positive, got $n"))
        sfs[end, 1] < n || throw(ArgumentError(
            "n = $n is not a valid haploid sample size: the SFS reaches DAC $(Int(sfs[end,1]))"))
        return n
    end
    m = size(sfs, 1)
    (sfs[1, 1] == 1 && sfs[end, 1] == m) || throw(ArgumentError(
        "the DAC column runs $(Int(sfs[1,1]))..$(Int(sfs[end,1])) over $m rows, so this SFS " *
        "was trimmed and the haploid sample size cannot be recovered from it. Pass " *
        "`n = <haploid sample size>`, or take alpha(x) from `parse_sfs`, which returns it " *
        "as [dac, frequency, alpha_x] with the frequency already correct."))
    return m + 1
end

"""
    α_x(sfs, divergence; cutoff, cumulative, n)

Estimate alpha(x) for each derived allele frequency.

# Arguments
 - `sfs::Matrix`: SFS parsed from `parse_sfs()`; column 1 is the derived allele count.
 - `divergence::Vector`: divergence data from `parse_sfs()`.
 - `cutoff::Vector{Float64}`: frequency window `[low, high]`. Default `[0.0, 1.0]`.
 - `cumulative::Bool`: estimate from the cumulative SFS. Default `true`.
 - `n::Union{Nothing,Int}`: haploid sample size. Default `nothing`, inferred as
   `size(sfs, 1) + 1`, with an error when the SFS was trimmed.

# Returns
 - `Matrix{Float64}`: `[dac, frequency, alpha_x]`, one row per retained frequency.
"""
function α_x(sfs::Matrix, divergence::Vector; cutoff::Vector{Float64}=[0.0, 1.0],
    cumulative::Bool=true, n::Union{Nothing,Int}=nothing)

    nn = _alpha_x_n(sfs, n)
    dac = Float64.(view(sfs, :, 1))
    f = dac ./ nn

    keep = (f .>= cutoff[1]) .& (f .<= cutoff[2])
    sfs_tmp = sfs[keep, :]
    cumulative && (sfs_tmp = cumulative_sfs(sfs_tmp))

    alpha_x = @. 1 - (sfs_tmp[:, 2] / sfs_tmp[:, 3] * divergence[2] / divergence[1])

    return hcat(dac[keep], f[keep], alpha_x)
end

"""
    α_x(sfs, divergence; kwargs...)

`α_x` for several SFS and divergence sets, in parallel.

# Returns
 - `Vector{Matrix{Float64}}`: one `[dac, frequency, alpha_x]` table per set.
"""
function α_x(sfs::Vector, divergence::Vector; kwargs...)
    return ThreadsX.mapi((x, y) -> α_x(x, y; kwargs...), sfs, divergence)
end

# p-value correction for the multi-gene methods. NaN p-values are left out of the
# correction (they would corrupt the ranking) and kept as NaN.
function _apply_correction(df::DataFrame, correction::Symbol)
    adjustment = if correction == :fdr
        BenjaminiHochberg()
    elseif correction == :bonferroni
        Bonferroni()
    else
        throw(ArgumentError("correction must be :fdr or :bonferroni, got :$correction"))
    end
    pvals = df.p_value
    n_total = length(pvals)
    valid = .!isnan.(pvals)
    adjusted = fill(NaN, n_total)
    if any(valid)
        adjusted[valid] = adjust(pvals[valid], n_total, adjustment)
    end
    pval_idx = columnindex(df, :p_value)
    insertcols!(df, pval_idx + 1, :pvalue_adjusted => adjusted)
    return df
end

################################
######   Asymptotic MK    ######
################################

"""
    AbstractAsymptoticModel

Curve fitted to alpha(x) by `asymptotic_mk`: [`ExponentialFit`](@ref) (default) or
[`LinearFit`](@ref). Only the chosen model is fitted; if it fails the result is NaN.
"""
abstract type AbstractAsymptoticModel end

"""
    ExponentialFit()

`alpha(x) = a + b·exp(-c·x)` (Messer & Petrov 2013). The asymptote is `a + b·exp(-c)`.
"""
struct ExponentialFit <: AbstractAsymptoticModel end

"""
    LinearFit()

`alpha(x) = a + b·x`. The asymptote is `a + b`.
"""
struct LinearFit <: AbstractAsymptoticModel end

_model_name(::ExponentialFit) = :exponential
_model_name(::LinearFit) = :linear

_amk_nan(model, f, α_i) = (asymptotic_mk=NaN, ci_low=NaN, ci_high=NaN,
    a=NaN, b=NaN, c=NaN, model=_model_name(model), f=f, alpha_x=α_i)

# Confidence interval of the asymptote. The asymptote is a function of the parameters,
# model(1, p), so its standard error comes from the delta method over the fit's
# covariance: se = sqrt(grad' * vcov * grad). `confint(fit)` would give one interval per
# parameter, not this one.
function _asymptote_ci(fit, grad::Vector{Float64}, asymp::Float64; level::Float64=0.95)
    return try
        se = sqrt(dot(grad, vcov(fit) * grad))
        isfinite(se) || return (NaN, NaN)
        t = quantile(TDist(dof(fit)), 1 - (1 - level) / 2)
        (asymp - t * se, asymp + t * se)
    catch
        (NaN, NaN)
    end
end

function _asymptotic_fit(model::ExponentialFit, f::Vector{Float64}, α_i::Vector{Float64})
    length(f) < 3 && return _amk_nan(model, f, α_i)

    exp_model(x, p) = @. p[1] + p[2] * exp(-x * p[3])
    fitted = try
        start = curve_fit(
            exp_model, f, α_i, [-1.0, -1.0, 1.0];
            lower=[-1.0, -1.0, 1.0], upper=[1.0, 1.0, 10.0],
        )
        curve_fit(exp_model, f, α_i, start.param)
    catch
        nothing
    end
    (isnothing(fitted) || !fitted.converged) && return _amk_nan(model, f, α_i)

    asymp = exp_model(1, fitted.param)
    a, b, c = fitted.param
    # gradient of a + b*exp(-c*x) at x = 1
    lo, hi = _asymptote_ci(fitted, [1.0, exp(-c), -b * exp(-c)], asymp)
    return (asymptotic_mk=asymp, ci_low=lo, ci_high=hi, a=a, b=b, c=c,
        model=_model_name(model), f=f, alpha_x=α_i)
end

function _asymptotic_fit(model::LinearFit, f::Vector{Float64}, α_i::Vector{Float64})
    length(f) < 3 && return _amk_nan(model, f, α_i)

    lin_model(x, p) = @. p[1] + p[2] * x
    fitted = try
        curve_fit(lin_model, f, α_i, [0.0, 0.0])
    catch
        nothing
    end
    (isnothing(fitted) || !fitted.converged) && return _amk_nan(model, f, α_i)

    a, b = fitted.param
    asymp = a + b
    # gradient of a + b*x at x = 1
    lo, hi = _asymptote_ci(fitted, [1.0, 1.0], asymp)
    return (asymptotic_mk=asymp, ci_low=lo, ci_high=hi, a=a, b=b, c=0.0,
        model=_model_name(model), f=f, alpha_x=α_i)
end

# Drop alpha(x) values no fit can use. alpha = 1 means Pn = 0, an absence of data.
function _amk_na_filter(f::Vector{Float64}, α_i::Vector{Float64})
    flt = .!isnan.(α_i) .&& .!isinf.(α_i) .&& α_i .!= 1
    return f[flt], α_i[flt]
end

"""
    asymptotic_mk(sfs, divergence; model, cutoff, cumulative, na_rm, n)
    asymptotic_mk(α_table; model, cutoff, na_rm)
    asymptotic_mk(α_x; model, f, cutoff, na_rm)
    asymptotic_mk(sfs::Vector, divergence::Vector; ...)

Asymptotic MK test (Messer & Petrov 2013): fit a curve to alpha(x) against the derived
allele frequency and take its value at x = 1.

alpha(x) is computed from `sfs`/`divergence`, taken from the `[dac, frequency, alpha_x]`
table returned by `parse_sfs` or `α_x`, or given as a vector. The `Vector` forms run per
gene or set in parallel and return a `DataFrame` with one row each.

# Arguments
 - `model::AbstractAsymptoticModel`: `ExponentialFit()` (default) or `LinearFit()`. Only
   this model is fitted; a fit that fails returns NaN.
 - `cutoff::Vector{Float64}`: frequency window `[low, high]`. Default `[0.0, 1.0]`. With
   `sfs`/`divergence` the SFS is cut before the cumulative sum; with a table or vector the
   alpha(x) values already computed are filtered. These differ under `cumulative = true`,
   so pass the cutoff where alpha(x) is computed.
 - `cumulative::Bool`: use the cumulative SFS (sfs/divergence only). Default `true`.
 - `n::Union{Nothing,Int}`: haploid sample size (sfs/divergence only). Needed when the SFS
   was trimmed by a cutoff; see `α_x`.
 - `f`: frequencies of the alpha(x) vector (vector form only). Default
   `(1:length(α_x)) / (length(α_x) + 1)`, i.e. a full spectrum.
 - `na_rm::Bool`: drop `NaN`, `Inf` and `1.0` alpha(x) values before fitting. Default `true`.

# Returns
 - `NamedTuple`: `(asymptotic_mk, ci_low, ci_high, a, b, c, model, f, alpha_x)`. The 95% CI
   is the interval of the asymptote. All NaN (except `model`, `f`, `alpha_x`) when fewer
   than three values remain or the fit fails. `c = 0` for `LinearFit`.
 - `DataFrame` with one such row per gene or set for the `Vector` forms.
"""
function asymptotic_mk(
    sfs::Matrix,
    divergence::Vector;
    model::AbstractAsymptoticModel=ExponentialFit(),
    cutoff::Vector{Float64}=[0.0, 1.0],
    cumulative::Bool=true,
    na_rm::Bool=true,
    n::Union{Nothing,Int}=nothing,
)
    tbl = α_x(sfs, divergence; cutoff=cutoff, cumulative=cumulative, n=n)
    f = Float64.(tbl[:, 2])
    α_i = Float64.(tbl[:, 3])

    na_rm && ((f, α_i) = _amk_na_filter(f, α_i))

    return _asymptotic_fit(model, f, α_i)
end

function asymptotic_mk(
    α_table::AbstractMatrix{<:Real};
    model::AbstractAsymptoticModel=ExponentialFit(),
    cutoff::Vector{Float64}=[0.0, 1.0],
    na_rm::Bool=true,
)
    size(α_table, 2) == 3 || throw(ArgumentError(
        "expected a 3-column [dac, frequency, alpha_x] table from `parse_sfs` or `α_x`, " *
        "got $(size(α_table, 2)) columns"))

    f = Float64.(α_table[:, 2])
    α_i = Float64.(α_table[:, 3])

    if cutoff != [0.0, 1.0]
        keep = (f .>= cutoff[1]) .& (f .<= cutoff[2])
        f = f[keep]
        α_i = α_i[keep]
    end

    na_rm && ((f, α_i) = _amk_na_filter(f, α_i))

    return _asymptotic_fit(model, f, α_i)
end

function asymptotic_mk(
    α_x_values::AbstractVector{<:Real};
    model::AbstractAsymptoticModel=ExponentialFit(),
    f::Union{Nothing,AbstractVector{<:Real}}=nothing,
    cutoff::Vector{Float64}=[0.0, 1.0],
    na_rm::Bool=true,
)
    α_i = collect(Float64.(α_x_values))
    _s = length(α_i)

    freq = if isnothing(f)
        collect((1:_s) ./ (_s + 1))
    else
        length(f) == _s || throw(ArgumentError(
            "f has $(length(f)) entries but alpha(x) has $_s; they must match"))
        collect(Float64.(f))
    end

    if cutoff != [0.0, 1.0]
        keep = (freq .>= cutoff[1]) .& (freq .<= cutoff[2])
        freq = freq[keep]
        α_i = α_i[keep]
    end

    na_rm && ((freq, α_i) = _amk_na_filter(freq, α_i))

    return _asymptotic_fit(model, freq, α_i)
end

# Several genes or sets: one row each, fitted in parallel.
function _asymptotic_mk_many(fit_one, n::Int)
    results = Vector{NamedTuple}(undef, n)
    Threads.@threads for i in 1:n
        results[i] = Logging.with_logger(fit_one(i), Logging.NullLogger())
    end
    df = DataFrame(results)
    n_nan = sum(isnan.(df.asymptotic_mk))
    n_nan > 0 && @warn "asymptotic_mk: fit failed for $n_nan gene(s). Returning NaN."
    return df
end

asymptotic_mk(sfs::Vector, divergence::Vector; kwargs...) =
    _asymptotic_mk_many(i -> () -> asymptotic_mk(sfs[i], divergence[i]; kwargs...), length(sfs))

asymptotic_mk(α_tables::AbstractVector{<:AbstractMatrix{<:Real}}; kwargs...) =
    _asymptotic_mk_many(i -> () -> asymptotic_mk(α_tables[i]; kwargs...), length(α_tables))

asymptotic_mk(α_x_values::AbstractVector{<:AbstractVector{<:Real}}; kwargs...) =
    _asymptotic_mk_many(i -> () -> asymptotic_mk(α_x_values[i]; kwargs...), length(α_x_values))

################################
######  Other MK tests    ######
################################

"""
    imputed_mk(sfs, divergence; cutoff)

Imputed MK test (Murga-Moreno et al. 2022, https://doi.org/10.1093/g3journal/jkac206):
slightly deleterious polymorphisms below `cutoff` are imputed and removed from Pn.

# Arguments
 - `sfs::Matrix`: SFS parsed from `parse_sfs()`.
 - `divergence::Vector`: divergence data from `parse_sfs()`.
 - `cutoff::Float64`: derived allele frequency separating low-frequency variants. Default `0.15`.

# Returns
 - `NamedTuple`: `(alpha, p_value, b, f, d, omega, omega_a, omega_d)`; `b`, `f`, `d` and
   the omegas need Ln and Ls in `divergence`.
"""
function imputed_mk(
    sfs::Matrix,
    divergence::Vector;
    cutoff::Float64=0.15,
)

    sfs = copy(sfs)
    if sum(@view sfs[:, 1]) > 1
        n = size(sfs, 1)
        sfs[:, 1] .= (1:n) ./ (n + 1)
    end

    pn = sum(sfs[:, 2])
    ps = sum(sfs[:, 3])
    dn = divergence[1]
    ds = divergence[2]

    deleterious = 0
    ### Estimating slightly deleterious with pn/ps ratio
    flt_low = (sfs[:, 1] .<= cutoff)
    pn_low = sum(sfs[flt_low, 2])
    ps_low = sum(sfs[flt_low, 3])

    flt_inter = (sfs[:, 1] .>= cutoff) .& (sfs[:, 1] .<= 1)
    pn_inter = sum(sfs[flt_inter, 2])
    ps_inter = sum(sfs[flt_inter, 3])

    ratio_ps = ps_inter > 0 ? ps_low / ps_inter : 0.0
    deleterious = ps_inter > 0 ? pn_low - (pn_inter * ratio_ps) : 0

    if (deleterious > pn) || (deleterious < 0)
        deleterious = 0
        pn_neutral = round(pn - deleterious, digits=3)
    else
        pn_neutral = round(pn - deleterious, digits=3)
    end

    alpha = round(1 - ((pn_neutral / ps) * (ds / dn)), digits=3)

    #  method = :minlike same results R, python two.sides
    p_value =
        pvalue(FisherExactTest(Int(ceil(ps)), Int(ceil(pn_neutral)), Int(ds), Int(dn)))

    b = f = d = omega = omega_a = omega_d = NaN
    if size(divergence, 1) > 2
        ln = divergence[3]
        ls = divergence[4]
        # ## Estimation of b: weakly deleterious
        b = (deleterious / ps) * (ls / ln)

        ## Estimation of f: neutral sites
        f = (ls * pn_neutral) / (ln * ps)

        ## Estimation of d, strongly deleterious sites
        d = 1 - (f + b)

        ka = dn / ln
        ks = ds / ls
        omega = ka / ks

        # omega_a and omega_d
        omega_a = omega * alpha
        omega_d = omega - omega_a
    end

    return (alpha=alpha, p_value=p_value, b=b, f=f, d=d, omega=omega, omega_a=omega_a, omega_d=omega_d)
end

"""
    imputed_mk(sfs::Vector, divergence::Vector; cutoff, correction)

`imputed_mk` for several genes or sets, in parallel.

# Arguments
 - `cutoff::Float64`: as for one set. Default `0.15`.
 - `correction::Union{Symbol,Nothing}`: p-value correction, `:fdr` or `:bonferroni`.
   Default `nothing`.

# Returns
 - `DataFrame` with one row per gene or set.
"""
function imputed_mk(
    sfs::Vector,
    divergence::Vector;
    cutoff::Float64=0.15,
    correction::Union{Symbol,Nothing}=nothing
)
    df = DataFrame(ThreadsX.mapi(sfs, divergence) do x, y
        Logging.with_logger(Logging.NullLogger()) do
            imputed_mk(x, y, cutoff=cutoff)
        end
    end)
    n_nan = sum(isnan.(df.alpha))
    n_nan > 0 && @warn "imputed_mk: insufficient data for $n_nan gene(s). Returning NaN."

    isnothing(correction) && return df

    return _apply_correction(df, correction)
end


"""
    fww_mk(sfs, divergence; cutoff)

Fay, Wyckoff and Wu MK test (https://doi.org/10.1038/4151024a): polymorphisms below
`cutoff` are removed.

# Arguments
 - `sfs::Matrix`: SFS parsed from `parse_sfs()`.
 - `divergence::Vector`: divergence data from `parse_sfs()`.
 - `cutoff::Float64`: minimum derived allele frequency. Default `0.15`.

# Returns
 - `NamedTuple`: `(alpha, p_value, omega, omega_a, omega_d)`.
"""
function fww_mk(
    sfs::Matrix,
    divergence::Vector;
    cutoff::Float64=0.15,
)

    # DAC to freqs
    sfs = copy(sfs)
    if sum(@view sfs[:, 1]) > 1
        n = size(sfs, 1)
        sfs[:, 1] .= (1:n) ./ (n + 1)
    end

    dn = divergence[1]
    ds = divergence[2]

    flt_inter = (sfs[:, 1] .>= cutoff) .& (sfs[:, 1] .<= 1)
    pn_inter = sum(sfs[flt_inter, 2])
    ps_inter = sum(sfs[flt_inter, 3])

    if ps_inter == 0 || dn == 0
        return (alpha=NaN, p_value=NaN, omega=NaN, omega_a=NaN, omega_d=NaN)
    end

    alpha = round(1 - ((pn_inter / ps_inter) * (ds / dn)), digits=3)
    #  method = :minlike same results R, python two.sides
    p_value =
        pvalue(FisherExactTest(Int(ps_inter), Int(ceil(pn_inter)), Int(ds), Int(dn)))

    omega = omega_a = omega_d = 0
    if size(divergence, 1) > 2
        ln = divergence[3]
        ls = divergence[4]
        ka = dn / ln
        ks = ds / ls

        omega = ka / ks

        # omega_a and omega_d
        omega_a = omega * alpha
        omega_d = omega - omega_a
    end

    return (alpha=alpha, p_value=p_value, omega=omega, omega_a=omega_a, omega_d=omega_d)
end

"""
    fww_mk(sfs::Vector, divergence::Vector; cutoff, correction)

`fww_mk` for several genes or sets, in parallel. `correction` is `:fdr`, `:bonferroni` or
`nothing` (default).

# Returns
 - `DataFrame` with one row per gene or set.
"""
function fww_mk(
    sfs::Vector,
    divergence::Vector;
    cutoff::Float64=0.15,
    correction::Union{Symbol,Nothing}=nothing
)
    df = DataFrame(ThreadsX.mapi(sfs, divergence) do x, y
        Logging.with_logger(Logging.NullLogger()) do
            fww_mk(x, y, cutoff=cutoff)
        end
    end)
    n_nan = sum(isnan.(df.alpha))
    n_nan > 0 && @warn "fww_mk: insufficient data for $n_nan gene(s). Returning NaN."

    isnothing(correction) && return df

    return _apply_correction(df, correction)
end

"""
    standard_mk(sfs, divergence)

Standard MK test: alpha = 1 - (Pn/Ps)(Ds/Dn), with a Fisher exact test.

# Arguments
 - `sfs::Matrix`: SFS parsed from `parse_sfs()`.
 - `divergence::Vector`: divergence data from `parse_sfs()`.

# Returns
 - `NamedTuple`: `(alpha, p_value, omega, omega_a, omega_d)`.
"""
function standard_mk(
    sfs::Matrix,
    divergence::Vector;
)

    # DAC to freqs
    sfs = copy(sfs)
    if (sfs[1, 1] >= 1)
        nn = size(view(sfs, :, 1), 1)
        sfs[:, 1] .= collect(1:nn) ./ (nn + 1)
    end

    pn = sum(sfs[:, 2])
    ps = sum(sfs[:, 3])
    dn = divergence[1]
    ds = divergence[2]

    if ps == 0 || dn == 0
        return (alpha=NaN, p_value=NaN, omega=NaN, omega_a=NaN, omega_d=NaN)
    end

    alpha = round(1 - ((pn / ps) * (ds / dn)), digits=3)
    #  method = :lnnlike same results R, python two.sides
    p_value = pvalue(FisherExactTest(Int(ps), Int(ceil(pn)), Int(ds), Int(dn)))

    omega = omega_a = omega_d = 0
    if size(divergence, 1) > 2
        ln = divergence[3]
        ls = divergence[4]

        ka = dn / ln
        ks = ds / ls
        omega = ka / ks

        # omega_a and omega_d
        omega_a = omega * alpha
        omega_d = omega - omega_a
    end
    return (alpha=alpha, p_value=p_value, omega=omega, omega_a=omega_a, omega_d=omega_d)
end

"""
    standard_mk(sfs::Vector, divergence::Vector; correction)

`standard_mk` for several genes or sets, in parallel. `correction` is `:fdr`,
`:bonferroni` or `nothing` (default).

# Returns
 - `DataFrame` with one row per gene or set.
"""
function standard_mk(
    sfs::Vector,
    divergence::Vector;
    correction::Union{Symbol,Nothing}=nothing
)
    df = DataFrame(ThreadsX.mapi(sfs, divergence) do x, y
        Logging.with_logger(Logging.NullLogger()) do
            standard_mk(x, y)
        end
    end)
    n_nan = sum(isnan.(df.alpha))
    n_nan > 0 && @warn "standard_mk: insufficient data for $n_nan gene(s). Returning NaN."

    isnothing(correction) && return df

    return _apply_correction(df, correction)
end



"""
    ni(sfs, divergence)
    ni(sfs::Vector, divergence::Vector; correction)

Neutrality index, NI = (Pn/Ps) / (Dn/Ds), with a Fisher exact test.

# Returns
 - `NamedTuple` `(ni, p_value)`, or a `DataFrame` with one row per gene or set.
"""
function ni(
    sfs::Matrix,
    divergence::Vector;
)

    pn = sum(@view sfs[:, 2])
    ps = sum(@view sfs[:, 3])
    dn = divergence[1]
    ds = divergence[2]

    if ds ≡ 0.0 || ps ≡ 0.0 || dn ≡ 0.0
        @warn "Cannot estimate neutrility index with ds, dn or ps = 0. Returning NaN values"
        _ni = NaN
        p_value = NaN
    else

        p_value = pvalue(FisherExactTest(Int(ps), Int(ceil(pn)), Int(ds), Int(dn)))

        _ni = (pn / ps) / (dn / ds)
    end
    return (ni=_ni, p_value=p_value)
end

function ni(
    sfs::Vector,
    divergence::Vector;
    correction::Union{Symbol,Nothing}=nothing
)
    df = DataFrame(ThreadsX.mapi(sfs, divergence) do x, y
        Logging.with_logger(Logging.NullLogger()) do
            ni(x, y)
        end
    end)
    n_nan = sum(isnan.(df.ni))
    n_nan > 0 && @warn "ni: insufficient data for $n_nan gene(s). Returning NaN."
    isnothing(correction) && return df
    return _apply_correction(df, correction)
end

"""
    ni_cmh(sfs, divergence)

Neutrality index over several genes with the Cochran-Mantel-Haenszel weighting of
Stoletzki & Eyre-Walker (2011).

"""
function ni_cmh(sfs::Vector, divergence::Vector)

    @assert length(sfs) ≡ length(divergence) "Must contains same number of genes"


    num = 0.0
    den = 0.0

    @inbounds for i in eachindex(sfs)
        @views begin
            pn = sum(sfs[i][:, 2])
            ps = sum(sfs[i][:, 3])
        end
        dn = divergence[i][1]
        ds = divergence[i][2]

        w = ps + ds
        if w > 0
            num += (ds * pn) / w
            den += (dn * ps) / w
        end
    end

    return (ni = num / den)

end

"""
    dos(sfs, divergence)
    dos(sfs::Vector, divergence::Vector; correction)

Direction of selection, DoS = Dn/(Dn+Ds) - Pn/(Pn+Ps), with a Fisher exact test.

# Returns
 - `NamedTuple` `(dos, p_value)`, or a `DataFrame` with one row per gene or set.
"""
function dos(
    sfs::Matrix,
    divergence::Vector;
)

    pn = sum(@view sfs[:, 2])
    ps = sum(@view sfs[:, 3])
    dn = divergence[1]
    ds = divergence[2]

    dn_ratio = dn / (dn + ds)
    pn_ratio = pn / (pn + ps)

    if dn_ratio == 0
        @warn "Cannot estimate DoS without divergence data. Returning NaN values"
        dos = NaN
        p_value = NaN
    else
        dos = dn_ratio - pn_ratio

        p_value = pvalue(FisherExactTest(Int(ps), Int(ceil(pn)), Int(ds), Int(dn)))
    end
    return (dos=dos, p_value=p_value)
end

function dos(
    sfs::Vector,
    divergence::Vector;
    correction::Union{Symbol,Nothing}=nothing
)
    df = DataFrame(ThreadsX.mapi(sfs, divergence) do x, y
        Logging.with_logger(Logging.NullLogger()) do
            dos(x, y)
        end
    end)
    n_nan = sum(isnan.(df.dos))
    n_nan > 0 && @warn "dos: missing divergence data for $n_nan gene(s). Returning NaN."
    isnothing(correction) && return df
    return _apply_correction(df, correction)
end

"""
    grapes(sfs,divergence,model)

Run grapes using SFS and divergence data. The function will install grapes using Conda from genomedk channel. 

# Arguments
 - `sfs::Vector{Matrix{Float64}}`: SFS data parsed from parse_sfs()
 - `divergence::Vector{Matrix{Int64}}`: divergence data parsed from parse_sfs()
 - `model::String`: grapes model (GammaZero, GammaExpo, DisplGamma, ScaledBeta, FGMBesselK and all)
 - `n::Int64`: smaller sample size to project the SFS
 - `nearly_neutral::Int64`: Sadp threshold of above which a mutation is considered adaptive
 - `FWW_threshold::Float64`: minimal allele frequency in FWW α
 - `nb_rand_start::Int64`: number of random starting values in model optimization (default=0); setting positive values will slow down the program but decrease the probability of being trapped in local optima.
 - `anc_to_rec_Ne_ratio::Float64`: divergence/polymorphism Ne ratio
 - `no_div_data::Bool`: only use divergence data to estimate DFE
 - `no_div_param::Bool`: implements so-called [-A] version in Galtier (2016); also called a DFE by Tataru et al. (2017), and indicated with * in Rousselle et al. (2018); irrelevant if model = GammaZero; (default=false)
 - `no_syn_orient_error::Bool`: force equal synonymous and non-synonymous mis-orientation rate(default=false)
 - `fold::Bool`: fold the SFS
 - `fixed_param::String`: this option should be used if one does not want to optimize every parameter, but rather have some parameters fixed to predefined values; parameter names and predefined values are passed via a control file; see example at the bottom of this file (default=none).
# Output
 - `DataFrame`: grapes model estimation.
"""
function grapes(
    sfs::Union{Vector,Matrix},
    divergence::Union{Vector,Matrix},
    model::String;
    n::Int64=0,
    nearly_neutral::Int64 = 5,
    FWW_threshold::Float64 = 0.15,
    nb_rand_start::Int64 = 0,
    anc_to_rec_Ne_ratio::Float64 = 1.0,
    no_div_data::Bool = false,
    no_div_param::Bool = false,
    no_syn_orient_error::Bool = false,
    fold::Bool = false,
    fixed_param::String = "",
)

    if any(0 .∈ divergence)
        throw(
            ArgumentError(
                "Your SFS contains 0 values at the selected DACs or the divergence is 0. Please consider to bin the SFS and re-estimate the rates using the selected bin as sample the new sample size.",
            ),
        )
    end
    grapes_bin = CondaPkg.which("grapes")

    if (isnothing(grapes_bin))
        CondaPkg.add("grapes-static", channel = "genomedk")
        grapes_bin = CondaPkg.which("grapes")
    end

    @assert model ∈ [
        "GammaZero",
        "GammaExpo",
        "GammaGamma",
        "DisplGamma",
        "ScaledBeta",
        "FGMBesselK",
        "all",
    ] "Please select a valid model: GammaZero GammaExpo GammaGamma DisplGamma ScaledBeta FGMBesselK all"

    # Change to Vectors of Matrix if needed
    if sfs isa Matrix
        sfs = [sfs]
        divergence = [divergence]
    end

    # Reduce to n-1 to input in Grapes as SFS = n - 1 entries
    if n != 0
        sfs_r = project.(sfs, n)
        n *= 2
    else
        sfs_r = sfs
        n     = size.(sfs,1) .+ 1
    end

    pn = map(x -> permutedims(x[:, 2]), sfs_r)
    ps = map(x -> permutedims(x[:, 3]), sfs_r)

    dn = map(x -> [x[1]], divergence)
    ds = map(x -> [x[2]], divergence)
    ln = map(x -> [x[3]], divergence)
    ls = map(x -> [x[4]], divergence)

    idx = string.(collect(1:length(sfs_r)))

    @info "Converting SFS to dofe file"
    dofe = @. sfs_to_dofe(pn, ps, dn, ds, ln, ls, idx, n)
    h = fill(DataFrame(["" ""; "#unfolded" ""], :auto), length(sfs_r))

    output_dofe   = @. tempname() * idx
    output_grapes = @. output_dofe * "." * model

    @. write_files(h, output_dofe, false)
    @. write_files(dofe, output_dofe, true)

    @info "Running Grapes"
    if model == "GammaZero" || no_div_param
        no_div_param = false
    end

    bool_options = ["-no_div_data", "-no_div_param", "-no_syn_orient_error", "-fold"]
    bool_options = join(
        bool_options[any(
            vcat(no_div_data, no_div_param, no_syn_orient_error, fold),
            dims = 2,
        )],
        " ",
    )

    if !isempty(fixed_param)
        fixed_param = "-fixed_param " * fixed_param
    end

    r(
        d,
        o,
        md = model,
        gr = grapes_bin,
        nt = nearly_neutral,
        fww = FWW_threshold,
        nb = nb_rand_start,
        ne = anc_to_rec_Ne_ratio,
        bopt = bool_options,
        fp = fixed_param,
    ) = run(
        `$gr -in $d -out $o -model $md -nearly_neutral $nt -FWW_threshold $fww -nb_rand_start $nb -anc_to_rec_Ne_ratio $ne $bopt $fp`,
    )

    @suppress_out begin
        ThreadsX.mapi((x,y) -> r(x,y), output_dofe, output_grapes; ntasks = Threads.nthreads())
    end

    df = map((x,y) -> read_clean_grapes(x,y,model),output_dofe,output_grapes)
    df = df[.!isnothing.(df)]

    return (vcat(df...))
end

function sfs_to_dofe(pn, ps, dn, ds, ln, ls, w, n)
    return(DataFrame(hcat("dofe_" * string(w), n, ln, pn, ls, ps, ln, dn, ls, ds...), :auto))
end

function read_clean_grapes(dofe::String,grapes_file::String,model::String)
    # Try to open file if grapes executed properly
    output = @suppress_err begin
        try
            if model == "all"
                df = CSV.read(grapes_file, DataFrame)
            else
                df = CSV.read(grapes_file, DataFrame, footerskip = 1, skipto = 3)

                idx_1 = findall(occursin.(model,names(df)))
                idx_2 = findfirst(occursin.("theta",names(df)))
                df = hcat(df[:,1:20],df[:,idx_1],df[:,idx_2:end])
            end

            rm(dofe)
            rm(grapes_file)
            rm(grapes_file * ".profile")
            rm(grapes_file * ".messages")

            df

        catch
            @warn "$dofe did not complete the Grapes estimation. Please check your data."
            nothing
        end
    end

    return(output)
end

"""
    abcmk(param::parameters,
          sfs::Vector,
          divergence::Vector;
          h5_file::String,
          summstat_size::Int64,
          output_folder::String,
          alpha::Union{Nothing,Vector{Float64}} = nothing,
          B_bins::Union{Nothing,Vector{Float64}} = nothing,
          tol::Int64=1000,
          abc_method::String="abcreg",
          rm_summaries::Bool=true,
          stat::String="median",
          sampling::AbstractPoissonSampling=RjSampling(),
          row_prior::AbstractRowPrior=FlatAlphaPrior()) -> DataFrame

Performs Approximate Bayesian Computation (ABC) using either regression adjustment (ABCreg) or a basic ABC method and returns a summary of the posterior distributions.

# Arguments
- `param::parameters`: The model parameters.
- `sfs::Vector`: Site frequency spectrum data.
- `divergence::Vector`: Divergence data.
- `h5_file::String`: Path to the HDF5 file containing the data.
- `summstat_size::Int64`: Size of the summary statistics.
- `output_folder::String`: Folder where outputs will be stored.
- `alpha::Union{Nothing,Vector{Float64}}`: Optional alpha values for the summary statistics function. Default is `nothing`.
- `B_bins::Union{Nothing,Vector{Float64}}`: Optional B bins for the summary statistics function. Default is `nothing`.
- `tol::Int64`: Number of reatained values in ABC. Defining the olerance level for the ABC method. Default is `1000`.
- `abc_method::String`: Method to use for ABC, either "abcreg" or another method. Default is `"abcreg"`.
- `rm_summaries::Bool`: Whether to remove summaries after computation. Default is `true`.
- `stat::String`: The statistic to summarize the posteriors, such as "mode". Default is `"median"`.
- `sampling::AbstractPoissonSampling`: how the summary statistics are sampled from the rates
  table. `RjSampling()` (default) or `IndependentSampling()`.
- `row_prior::AbstractRowPrior`: which table rows are drawn. `FlatAlphaPrior()` (default) or
  `TableRowPrior()`.

The defaults (`RjSampling()`, `FlatAlphaPrior()`, `"median"`) are the configuration measured
against SLiM in the two-epoch and Tennessen benchmarks (see [`RjSampling`](@ref)). ABC-MK
v1.0.0 is `sampling = IndependentSampling(), row_prior = TableRowPrior(), stat = "mode"`.

# Returns
- `DataFrame`: A DataFrame summarizing the posterior distributions.
"""
function abcmk(
    param::parameters,
    sfs::Vector,
    divergence::Vector;
    h5_file::String,
    summstat_size::Int64,
    output_folder::String,
    alpha::Union{Nothing,Vector{Float64}} = nothing,
    B_bins::Union{Nothing,Vector{Float64}} = nothing,
    tol::Int64=1000,
    abc_method::String="abcreg",
    rm_summaries::Bool=true,
    stat::String="median",
    sampling::AbstractPoissonSampling=RjSampling(),
    row_prior::AbstractRowPrior=FlatAlphaPrior(),
)
    @assert abc_method ∈ ["abcreg","abc"] "Please select ABCreg or abc from R package"

    summ_stats = summary_statistics(param,sfs,divergence,h5_file=h5_file,output_folder=output_folder,summstat_size=summstat_size,alpha=alpha,B_bins=B_bins,sampling=sampling,row_prior=row_prior);

    if lowercase(abc_method) == "abcreg"
        posteriors = ABCreg(output_folder=output_folder,S=length(param.dac),tol=tol,rm_summaries=rm_summaries);
    else
        posteriors_adj_unadj = abc(output_folder=output_folder,S=length(param.dac),tol=tol,rm_summaries=rm_summaries);
        posteriors = posteriors_adj_unadj["loclinear"]
        if all(isnan.(posteriors[1][:,1]))
            @. replace!(posteriors, NaN=>0)
        end
    end

    df_abcmk = summary_abc(posteriors,stat=stat)

    return df_abcmk
end

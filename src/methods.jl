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
    abcmk(param::parameters,
          sfs::Vector,
          divergence::Vector;
          h5_file::String,
          summstat_size::Real,
          output_folder::String,
          alpha::Union{Nothing,Vector{Float64}} = nothing,
          B_bins::Union{Nothing,Vector{Float64}} = nothing,
          tol::Real=1000,
          abc_method::String="abcreg",
          rm_summaries::Bool=true,
          stat::String="median",
          sampling::AbstractPoissonSampling=RjSampling(),
          row_prior::AbstractRowPrior=FlatAlphaPrior()) -> DataFrame

Performs Approximate Bayesian Computation (ABC) using either regression adjustment (ABCreg) or a basic ABC method and returns a summary of the posterior distributions.

# Arguments
- `param::parameters`: The model parameters.
- `sfs::Vector`: Site frequency spectrum data, one matrix per set. A single `Matrix` (one set,
  as `parse_sfs` returns it) is also accepted.
- `divergence::Vector`: Divergence data, one vector per set, or a single `Vector{<:Real}` with a
  single `Matrix` SFS.
- `h5_file::String`: Path to the HDF5 file containing the data.
- `summstat_size::Real`: Number of summary statistics rows to sample, a positive whole number (`1e5` works).
- `output_folder::String`: Folder where outputs will be stored.
- `alpha::Union{Nothing,Vector{Float64}}`: Optional alpha values for the summary statistics function. Default is `nothing`.
- `B_bins::Union{Nothing,Vector{Float64}}`: Optional B bins for the summary statistics function. Default is `nothing`.
- `tol::Real`: Number of retained values in ABC, a positive whole number (`1e3` works). Default is `1000`.
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
    summstat_size::Real,
    output_folder::String,
    alpha::Union{Nothing,Vector{Float64}} = nothing,
    B_bins::Union{Nothing,Vector{Float64}} = nothing,
    tol::Real=1000,
    abc_method::String="abcreg",
    rm_summaries::Bool=true,
    stat::String="median",
    sampling::AbstractPoissonSampling=RjSampling(),
    row_prior::AbstractRowPrior=FlatAlphaPrior(),
)
    @assert abc_method ∈ ["abcreg","abc"] "Please select ABCreg or abc from R package"
    summstat_size = _as_count(summstat_size, "summstat_size")
    tol = _as_count(tol, "tol")

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

# One set, as `parse_sfs` returns it for a single gene list: `sfs::Matrix`, `divergence::Vector`.
abcmk(param::parameters, sfs::AbstractMatrix{<:Real}, divergence::AbstractVector{<:Real}; kwargs...) =
    abcmk(param, [Matrix{Float64}(sfs)], [Vector{Float64}(divergence)]; kwargs...)

# >>> dfe_ml entry points (added with the engine; test/verbatim_gate.sh skips this block)
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
    dfe_ml(sfs, divergence, model::AbstractDFEModel; kwargs...) -> DFEResult

Maximum-likelihood fit of the DFE `model` and alpha to the SFS and divergence of one
or several datasets, with GRAPES' model (Galtier 2016). One fit per dataset.

# Arguments
 - `sfs`: a matrix `[DAC Pn Ps]` with every bin `1..n-1` (unfolded; `n` chromosomes),
   as `parse_sfs` returns it, or a Vector of such matrices.
 - `divergence`: `[Dn, Ds, Ln, Ls]`, or a Vector of those (one per SFS).
 - `model`: an `AbstractDFEModel` instance (`GammaZero()`, `GammaExpo()`,
   `WithLethal(GammaZero())`, a user model, ...).

# Keywords (data)
 - `project_sfs::Union{Nothing,Int} = nothing`: project the SFS to this many diploid
   individuals (`2 project_sfs` chromosomes) with `project` first.
 - `fold::Bool = false`: fold the SFS (`misorient` is then not fitted).
 - `labels = nothing`: names of the datasets (the `dataset` column); default `1:N`.
 - `Ln_poly`, `Ls_poly`: site counts of the polymorphism data, when they differ from
   the divergence ones (default: `Ln`, `Ls`).

# Keywords (model axes; each is independent of the others)
 - `ne_model::AbstractNeModel = ConstantNe()`: `AncestralNeRatio(r)` or
   `AncestralNeRatio(:free)` (adds the parameter `ne_ratio`).
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
 - `summaries = (DFEMAlpha(),)`: the `AbstractAlphaSummary` column sets of `inference`.

# Keywords (fit)
 - `fixed::NamedTuple = (;)`: pin parameters by name, e.g. `(ben_prop = 0.0,)`.
 - `bounds::NamedTuple = (;)`: replace bounds by name, e.g. `(ben_prop = (1e-6, 0.99),)`.
 - `no_syn_orient_error::Bool = false`: do not fit `misorient`.
 - `nb_rand_start::Int = 0`: extra random starting points (seeded, reproducible).
 - `nearly_neutral::Real = 5`: S threshold of adaptive substitutions in alpha.
 - `alpha_ci_threshold::Real = 1.92`: log-likelihood drop of the alpha CI (half the
   chi-square(1) 95% quantile; GRAPES' interval). 0 disables it and the CI columns are NaN.
   The CI exists only under `SaturatedDivergence`, as in GRAPES.
 - `cache::DFECache = DFECache()`: precomputed tables; pass one to reuse them across
   calls (safe across threads).

# Returns
A [`DFEResult`](@ref): `inference` (alpha, omegas, lnL, convergence: one row per fit),
`params` (every fitted parameter, long), `sfs` (observed and fitted SFS, `r_j`),
`divergence` (expected Dn by S class), `config` (the settings) and `extra`. Parameter names
follow `<component>_<quantity>` (`del_shape`, `del_mean`, `ben_prop`, `ben_mean`, ...,
`misorient`); `from_legacy_name` maps the GRAPES / pre-2026-10-15 names.
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
                summaries = (DFEMAlpha(),),
                fixed::NamedTuple = NamedTuple(),
                bounds::NamedTuple = NamedTuple(),
                no_syn_orient_error::Bool = false,
                nb_rand_start::Int = 0,
                nearly_neutral::Real = 5,
                alpha_ci_threshold::Real = 1.92,
                labels = nothing,
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

    datasets = labels === nothing ? collect(1:length(sfs_l)) : collect(labels)
    length(datasets) == length(sfs_l) ||
        throw(ArgumentError("got $(length(datasets)) labels for $(length(sfs_l)) datasets"))

    pieces = Vector{Any}(undef, length(sfs_l))
    for i in eachindex(sfs_l)
        n_i = Int(round(maximum(sfs_l[i][:, 1]))) + 1
        pre = get_precomputed(cache, n_i; integrator = integrator, grid = grid,
                              dominance = dominance, expectation = expectation)
        data = build_data_input(sfs_l[i], div_l[i], fold, uses_divergence(div_model);
                                Ln_poly = Float64(Ln_poly), Ls_poly = Float64(Ls_poly),
                                neutral_shape = pre.neutral_shape)
        pieces[i] = fit_model(data, pre, model; kernel = kernel, quad = quadrature, optimizer = optimizer,
                              ne_model = ne_model, div_model = div_model, demography = dem,
                              no_syn_orient_error = no_syn_orient_error, fixed = fixed, bounds = bounds,
                              nb_rand_start = nb_rand_start, nearly_neutral = Float64(nearly_neutral),
                              alpha_ci_threshold = Float64(alpha_ci_threshold), summaries = summaries)
    end
    config = OrderedDict{Symbol,Any}(
        :expectation => _axis_string(expectation), :integrator => _axis_string(integrator),
        :quadrature => _axis_string(quadrature), :s_grid => _axis_string(grid),
        :dominance => _axis_string(dominance), :demography => _axis_string(dem),
        :ne_model => _axis_string(ne_model), :div_model => _axis_string(div_model),
        :kernel => _axis_string(kernel), :optimizer => _axis_string(optimizer),
        :nearly_neutral => Float64(nearly_neutral), :alpha_ci_threshold => Float64(alpha_ci_threshold),
        :ci_method => "likelihood_ratio", :nb_rand_start => nb_rand_start,
        :fixed => _axis_string(fixed), :bounds => _axis_string(bounds),
        :no_syn_orient_error => no_syn_orient_error, :fold => fold,
        :project_sfs => project_sfs === nothing ? 0 : project_sfs,
        :mktest_version => string(pkgversion(@__MODULE__)))
    return _assemble_result(pieces, datasets, config)
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


"""
    dfe_model_comparison(r::DFEResult) -> (aic = DataFrame, lrt = DataFrame)
    dfe_model_comparison(sfs, divergence; models, cache = DFECache(), kwargs...)
        -> (aic = DataFrame, lrt = DataFrame, result = DFEResult)

AIC (with Akaike weights) and likelihood-ratio tests between DFE models fitted to the same
data under the same likelihood function: fits are grouped by `dataset` and by the settings
that define the likelihood (expectation, integrator, quadrature, grid, dominance,
demography, Ne model, divergence model, kernel, fold, projection), so models fitted with
different optimizers, starts, `fixed` values or `bounds` are compared, and models fitted
under different integrators are not. Model `i` is nested in `j` when its free parameters
are a strict subset of `j`'s (GammaZero in GammaExpo, ExpoZero in GammaZero, GammaZero in
`WithLethal(GammaZero())`); non-converged fits stay in the AIC table but not in the LRTs.
Both tables carry `dataset` and `group` (the likelihood setting) and are sorted by
log-likelihood within a group: pick the AIC-best model with `argmin(aic)`, not row 1.

The second form fits `models` (default: the eleven continuous shapes) with `dfe_ml`,
sharing one `DFECache`, and forwards `kwargs` to it. `fixed` and `bounds` are given to
each model restricted to the parameters it has, so `bounds = (ben_prop = (1e-6, 0.99),)`
widens the beneficial class of the models that have one.
"""
function dfe_model_comparison(r::DFEResult)
    inf = r.inference
    keycfg = select(r.config, :config, intersect(collect(_LIKELIHOOD_SETTINGS), propertynames(r.config)))
    lik = Dict(row.config => Tuple(row[Not(:config)]) for row in eachrow(keycfg))
    groups = unique([lik[c] for c in inf.config])
    free = Dict(f => Set(r.params.parameter[(r.params.fit .== f) .& r.params.free]) for f in inf.fit)
    aics, lrts = DataFrame[], DataFrame[]
    for ds in unique(inf.dataset), (gi, g) in enumerate(groups)
        sub = inf[(inf.dataset .== ds) .& ([lik[c] == g for c in inf.config]), :]
        nrow(sub) == 0 && continue
        names = String.(sub.model)
        allunique(names) || throw(ArgumentError("dataset $ds: model names repeat within one likelihood setting: $names"))
        fs = [free[f] for f in sub.fit]
        nested(i, j) = fs[i] ⊊ fs[j]
        extra = Pair{Symbol,Vector}[:fit => collect(sub.fit), :converged => collect(Bool, sub.converged)]
        for c in (:alpha, :omega_a)
            hasproperty(sub, c) && push!(extra, c => collect(sub[!, c]))
        end
        t = _aic_lrt_tables(names, collect(Int, sub.n_params), collect(Float64, sub.loglik), nested;
                            extra_cols = extra, skip = .!collect(Bool, sub.converged))
        insertcols!(t.aic, 1, :dataset => fill(ds, nrow(t.aic)), :group => fill(gi, nrow(t.aic)))
        insertcols!(t.lrt, 1, :dataset => fill(ds, nrow(t.lrt)), :group => fill(gi, nrow(t.lrt)))
        push!(aics, t.aic); push!(lrts, t.lrt)
    end
    return (aic = vcat(aics...; cols = :union), lrt = vcat(lrts...; cols = :union))
end

const _DEFAULT_COMPARISON_MODELS = (ExpoZero(), GammaZero(), ExpoExpo(), GammaExpo(), ExpoGamma(),
                                    GammaGamma(), DisplGamma(), ScaledBeta(), LogNormalZero(),
                                    NeuGamma(), FGMBesselK())

function dfe_model_comparison(sfs, divergence;
                              models = collect(AbstractDFEModel, _DEFAULT_COMPARISON_MODELS),
                              fixed::NamedTuple = NamedTuple(), bounds::NamedTuple = NamedTuple(),
                              cache::DFECache = DFECache(), kwargs...)
    # `fixed` and `bounds` reach each model restricted to the names it has (a bound on
    # `ben_prop` applies to GammaExpo, not to GammaZero); a name no model has is an error.
    names_of(m) = (n = MKtest.dfe_param_names(m); push!(n, :misorient);
                   haskey(kwargs, :kernel) && append!(n, collect(keys(kernel_parameters(kwargs[:kernel])))); n)
    allnames = union((names_of(m) for m in models)...)
    for k in (keys(fixed)..., keys(bounds)...)
        k in allnames || throw(ArgumentError("$k is not a parameter of any model in the comparison"))
    end
    pick(nt, m) = (n = names_of(m); NamedTuple(k => v for (k, v) in pairs(nt) if k in n))
    fits = DFEResult[]
    for m in models
        @info "Fitting $(dfe_name(m))"
        push!(fits, dfe_ml(sfs, divergence, m; cache = cache, fixed = pick(fixed, m),
                           bounds = pick(bounds, m), kwargs...))
    end
    result = vcat(fits...)
    t = dfe_model_comparison(result)
    return (aic = t.aic, lrt = t.lrt, result = result)
end

"""
    write_dfe_result(prefix, r::DFEResult) -> Vector{String}

Write every table of `r` as CSV: `<prefix>_inference.csv`, `_params.csv`, `_sfs.csv`,
`_divergence.csv`, `_config.csv`, and `_extra_<name>.csv` per extra table. Interval
columns (`(lo, hi)` tuples) are written as `<name>_low` / `<name>_high`. Returns the paths.
For one file with the tuples intact, `JLD2.save("r.jld2", "result", r)`.
"""
function write_dfe_result(prefix::AbstractString, r::DFEResult)
    paths = String[]
    for k in (:inference, :params, :sfs, :divergence, :config)
        f = "$(prefix)_$(k).csv"
        CSV.write(f, _split_intervals(getfield(r, k)))
        push!(paths, f)
    end
    for (name, t) in r.extra
        f = "$(prefix)_extra_$(name).csv"
        CSV.write(f, _split_intervals(t))
        push!(paths, f)
    end
    return paths
end

"""
    read_dfe_result(prefix) -> DFEResult

Read the tables written by `write_dfe_result`, joining `<name>_low` / `<name>_high` back
into `(lo, hi)` tuples.
"""
function read_dfe_result(prefix::AbstractString)
    rd(k) = _join_intervals(CSV.read("$(prefix)_$(k).csv", DataFrame))
    dir, base = splitdir(prefix)
    dir = isempty(dir) ? "." : dir
    extra = OrderedDict{Symbol,DataFrame}()
    for f in sort(readdir(dir))
        pre, suf = base * "_extra_", ".csv"
        (startswith(f, pre) && endswith(f, suf)) || continue
        extra[Symbol(f[length(pre)+1:end-length(suf)])] = _join_intervals(CSV.read(joinpath(dir, f), DataFrame))
    end
    return DFEResult(rd(:inference), rd(:params), rd(:sfs), rd(:divergence), rd(:config), extra)
end

# <<< dfe_ml entry points

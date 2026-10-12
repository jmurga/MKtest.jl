module MKtestPlotsExt

# plot_α_x and plot_asymptotic_mk, loaded when Plots is.

using MKtest, Plots
import MKtest: plot_α_x, plot_asymptotic_mk, α_x, asymptotic_mk, project, quantile, DataFrame

# Formats savefig writes with the GR backend.
const _PLOT_EXTS = (".png", ".svg", ".pdf", ".ps", ".html")

# Anything alpha(x)-shaped as a Vector{Vector{Float64}}.
function _ax_curves(x::AbstractMatrix{<:Real})
    size(x, 2) == 3 || throw(ArgumentError(
        "expected a 3-column [dac, frequency, alpha_x] table from `parse_sfs` or `α_x`, " *
        "got $(size(x, 2)) columns"))
    return Vector{Float64}[Float64.(x[:, 3])]
end

_ax_curves(x::AbstractVector{<:AbstractMatrix{<:Real}}) =
    reduce(vcat, (_ax_curves(m) for m in x); init=Vector{Float64}[])

_ax_curves(x::AbstractVector{<:Real}) = Vector{Float64}[Float64.(x)]
_ax_curves(x::AbstractVector{<:AbstractVector{<:Real}}) = Vector{Float64}[Float64.(c) for c in x]

function _ax_curves(x::NamedTuple)
    haskey(x, :alpha_x) || throw(ArgumentError(
        "NamedTuple has no `alpha_x` field (got $(keys(x))); pass an `asymptotic_mk` result"))
    return _ax_curves(x.alpha_x)
end

function _ax_curves(x::DataFrame)
    hasproperty(x, :alpha_x) || throw(ArgumentError(
        "DataFrame has no `alpha_x` column (got $(names(x))); pass a multi-gene `asymptotic_mk` result"))
    out = Vector{Float64}[]
    for v in x.alpha_x
        append!(out, _ax_curves(v))
    end
    return out
end

function _ax_curves(x)
    applicable(iterate, x) || throw(ArgumentError(
        "cannot interpret $(typeof(x)) as alpha(x) data; expected a vector of numbers, a " *
        "vector of such vectors, an `asymptotic_mk` result, or an iterable of either"))
    els = collect(x)
    isempty(els) && throw(ArgumentError("alpha(x) input is empty"))
    if all(e -> e isa Real, els)
        return Vector{Float64}[Float64.(els)]
    elseif all(e -> e isa AbstractVector, els)
        return Vector{Float64}[Float64.(collect(e)) for e in els]
    end
    throw(ArgumentError(
        "iterable of $(eltype(els)) is neither a curve nor a collection of curves"))
end

# Frequencies of each curve of `_ax_curves(x)`, or `nothing` where the input does not carry
# them (the position i is then DAC i, frequency i / n).
_ax_freqs(x::AbstractMatrix{<:Real}) = Union{Nothing,Vector{Float64}}[Float64.(x[:, 2])]
_ax_freqs(x::AbstractVector{<:AbstractMatrix{<:Real}}) =
    Union{Nothing,Vector{Float64}}[Float64.(m[:, 2]) for m in x]
_ax_freqs(x::NamedTuple) =
    Union{Nothing,Vector{Float64}}[haskey(x, :f) ? Float64.(x.f) : nothing]
_ax_freqs(x::DataFrame) = hasproperty(x, :f) ?
    Union{Nothing,Vector{Float64}}[Float64.(v) for v in x.f] :
    Union{Nothing,Vector{Float64}}[nothing for _ in 1:size(x, 1)]
_ax_freqs(x) = Union{Nothing,Vector{Float64}}[nothing for _ in _ax_curves(x)]

# x axis: the DAC of a [dac, frequency, alpha_x] table, otherwise the row index.
_ax_coords(x::AbstractMatrix{<:Real}) = size(x, 2) == 3 ? Float64.(x[:, 1]) : nothing
_ax_coords(x) = nothing

# Percentile band of the control curves at each DAC. Non-finite values (Ps = 0 in a small
# set) are skipped; a DAC with no finite value is left as a gap.
function _ax_band(controls::Vector{Vector{Float64}}, n::Int, percentil)
    length(percentil) == 2 || throw(ArgumentError(
        "percentil must hold exactly two values (low, high); got $(length(percentil))"))
    p_lo, p_hi = minimum(percentil), maximum(percentil)
    (0.0 <= p_lo && p_hi <= 1.0) ||
        throw(ArgumentError("percentil values must lie in [0, 1]; got $percentil"))
    isempty(controls) && throw(ArgumentError("control set is empty"))

    short = findfirst(c -> length(c) < n, controls)
    isnothing(short) || throw(ArgumentError(
        "control curve $short has length $(length(controls[short])) but the focal curve has " *
        "length $n; every control must be at least as long as the focal curve"))

    m = length(controls)
    low = Vector{Float64}(undef, n)
    high = Vector{Float64}(undef, n)
    buf = Vector{Float64}(undef, m)

    @inbounds for i in 1:n
        k = 0
        for j in 1:m
            v = controls[j][i]
            isfinite(v) || continue
            k += 1
            buf[k] = v
        end
        if k == 0
            low[i] = NaN
            high[i] = NaN
            continue
        end
        v = view(buf, 1:k)
        sort!(v)
        low[i] = quantile(v, p_lo, sorted=true)
        high[i] = quantile(v, p_hi, sorted=true)
    end
    return low, high
end

# Legend labels by position, padded with α1, α2, ... when there are fewer labels than curves.
function _ax_labels(labels, k::Int)
    isnothing(labels) && return ["α$i" for i in 1:k]
    labs = String[String(l) for l in labels]
    length(labs) >= k && return labs[1:k]
    return vcat(labs, ["α$i" for i in (length(labs)+1):k])
end

_ax_compat_labels(labels, l1, l2) =
    !isnothing(labels) ? labels :
    (isnothing(l1) && isnothing(l2)) ? nothing :
    [something(l1, "α1"), something(l2, "α2")]

function _check_output_path(output_path)
    ext = lowercase(splitext(output_path)[2])
    ext in _PLOT_EXTS || throw(ArgumentError(
        "cannot infer a figure format from \"$output_path\": extension " *
        (isempty(ext) ? "is missing" : "`$ext` is not supported") *
        ". Use one of " * join(_PLOT_EXTS, ", ")))
end

function _plot_α_x(
    curves::Vector{Vector{Float64}};
    control=nothing,
    percentil::AbstractVector{<:Real}=[0.1, 0.9],
    labels=nothing,
    x=nothing,
    freqs=nothing,
    cutoff::AbstractVector{<:Real}=[0.0, 1.0],
    n::Union{Nothing,Integer}=nothing,
    output_path::Union{Nothing,AbstractString}=nothing,
    kwargs...,
)
    isempty(curves) && throw(ArgumentError("no alpha(x) curve to plot"))
    length(cutoff) == 2 && 0 <= cutoff[1] <= cutoff[2] <= 1 || throw(ArgumentError(
        "cutoff must be [low, high] with 0 <= low <= high <= 1; got $cutoff"))
    n_max = maximum(length, curves)
    xc = isnothing(x) ? nothing : Float64.(x)
    labs = _ax_labels(labels, length(curves))
    has_band = !isnothing(control)

    # Positions of each curve inside the frequency cutoff. Without known frequencies,
    # position i is DAC i of a haploid sample of `n` (default: longest curve + 1).
    nn = isnothing(n) ? n_max + 1 : Int(n)
    nn > n_max || throw(ArgumentError(
        "n = $nn is not a valid haploid sample size for an alpha(x) curve of length $n_max"))
    fs = isnothing(freqs) ? fill(nothing, length(curves)) : freqs
    _freq(len, f) = isnothing(f) ? (1:len) ./ nn : view(f, 1:len)
    _keep(len, f) = findall(v -> cutoff[1] <= v <= cutoff[2], _freq(len, f))
    keeps = [_keep(length(c), fs[i]) for (i, c) in enumerate(curves)]
    all(isempty, keeps) && throw(ArgumentError(
        "no alpha(x) value lies inside the frequency cutoff $cutoff"))
    _xs(k) = isnothing(xc) ? Float64.(k) : xc[k]
    n_last = maximum(k -> isempty(k) ? 0 : last(k), keeps)

    # Keywords from the caller override these defaults.
    kw = values(kwargs)
    defaults = (;
        title="αₓ values",
        lw=2,
        ylabel="α",
        xlabel="DAC",
        legend=(length(curves) == 1 && !has_band) ? false : :bottomright,
    )
    if isnothing(xc) && !any(k -> haskey(kw, k), (:xlim, :xlims, :xlimit, :xlimits))
        defaults = merge(defaults, (; xticks=0:max(1, round(Int, n_last / 8)):n_last))
    end
    opts = merge(defaults, kw)

    p = nothing

    # The control band is drawn first, behind the curves.
    if has_band
        low, high = _ax_band(_ax_curves(control), n_max, percentil)
        # The band follows the frequencies of the longest curve.
        kb = _keep(n_max, fs[argmax(length.(curves))])
        lbl = "control " * string(round(100 * minimum(percentil), digits=1)) * "-" *
              string(round(100 * maximum(percentil), digits=1)) * "%"
        band = (; fillrange=high[kb], fillalpha=0.25, linealpha=0, color=:grey, label=lbl)
        p = plot(_xs(kb), low[kb]; merge(opts, band)...)
    end

    colors = palette(:tab10)
    for (i, c) in enumerate(curves)
        k = keeps[i]
        series = (; label=labs[i], color=colors[mod1(i, 10)])
        if isnothing(p)
            p = plot(_xs(k), c[k]; merge(opts, series)...)
        else
            plot!(p, _xs(k), c[k]; lw=opts.lw, series...)
        end
    end

    if !isnothing(output_path)
        _check_output_path(output_path)
        savefig(p, output_path)
    end

    return p
end

function plot_α_x(
    a_x::AbstractVector{<:Real}...;
    labels=nothing,
    label_1::Union{Nothing,AbstractString}=nothing,
    label_2::Union{Nothing,AbstractString}=nothing,
    kwargs...,
)
    curves = Vector{Float64}[]
    for v in a_x
        append!(curves, _ax_curves(v))
    end
    return _plot_α_x(curves; labels=_ax_compat_labels(labels, label_1, label_2), kwargs...)
end

plot_α_x(curves::AbstractVector{<:AbstractVector{<:Real}}; kwargs...) =
    _plot_α_x(_ax_curves(curves); kwargs...)

plot_α_x(α_table::AbstractMatrix{<:Real}; kwargs...) =
    _plot_α_x(_ax_curves(α_table); x=_ax_coords(α_table), freqs=_ax_freqs(α_table), kwargs...)

plot_α_x(α_tables::AbstractVector{<:AbstractMatrix{<:Real}}; kwargs...) =
    _plot_α_x(_ax_curves(α_tables); x=_ax_coords(first(α_tables)), freqs=_ax_freqs(α_tables),
              kwargs...)

plot_α_x(result::NamedTuple; kwargs...) =
    _plot_α_x(_ax_curves(result); freqs=_ax_freqs(result), kwargs...)
plot_α_x(result::DataFrame; kwargs...) =
    _plot_α_x(_ax_curves(result); freqs=_ax_freqs(result), kwargs...)

plot_α_x(
    a_x::AbstractVector{<:Real},
    control::AbstractVector{<:AbstractVector{<:Real}};
    kwargs...,
) = _plot_α_x(_ax_curves(a_x); control=control, kwargs...)

function plot_α_x(
    sfs::Matrix,
    divergence::Vector;
    sample_size::Union{Nothing,Int}=nothing,
    cumulative::Bool=true,
    cutoff::AbstractVector{<:Real}=[0.0, 1.0],
    n::Union{Nothing,Int}=nothing,
    kwargs...,
)
    if !isnothing(sample_size)
        sfs = project(sfs, sample_size)
    end
    # Cut the SFS before the cumulative sum, as `asymptotic_mk` does.
    tbl = α_x(sfs, divergence; cutoff=Float64.(cutoff), cumulative=cumulative, n=n)
    return _plot_α_x(_ax_curves(tbl); x=_ax_coords(tbl), freqs=_ax_freqs(tbl), kwargs...)
end

function plot_asymptotic_mk(
    result::NamedTuple;
    output_path::Union{Nothing,AbstractString}=nothing,
)
    p = plot(result.f, result.alpha_x,
        label="Observed α(x)",
        color=:navy,
        linewidth=2,
        xlabel="Derived allele frequency x",
        ylabel="MK α(x)",
        title="Asymptotic MK Test: α(x) vs Frequency",
        legend=:bottomright,
        xlims=(0, 1),
    )

    if !isnan(result.asymptotic_mk)
        if !isnan(result.ci_low) && !isnan(result.ci_high)
            plot!(p, [0.0, 1.0], [result.ci_low, result.ci_low],
                fillrange=result.ci_high,
                fillalpha=0.3,
                alpha=0,
                color=:gray,
                label="95% CI ($(round(result.ci_low, digits=2)) - $(round(result.ci_high, digits=2)))",
            )
        end

        hline!(p, [result.asymptotic_mk],
            linestyle=:dash,
            color=:red,
            linewidth=2,
            label="α_asym = $(round(result.asymptotic_mk, digits=2))",
        )

        x_curve = range(0.05, 1.0, length=100)
        if result.model == :exponential
            y_curve = result.a .+ result.b .* exp.(-result.c .* x_curve)
            fit_label = "Fit: a + b·exp(−cx)"
        else
            y_curve = result.a .+ result.b .* x_curve
            fit_label = "Fit: a + b·x"
        end
        plot!(p, x_curve, y_curve, color=:red, linewidth=2, label=fit_label)
    end

    if !isnothing(output_path)
        _check_output_path(output_path)
        savefig(p, output_path)
    end
    return p
end

function plot_asymptotic_mk(
    sfs::Matrix,
    divergence::Vector;
    model::MKtest.AbstractAsymptoticModel=MKtest.ExponentialFit(),
    cutoff::Vector{Float64}=[0.0, 1.0],
    cumulative::Bool=true,
    n::Union{Nothing,Int}=nothing,
    output_path::Union{Nothing,AbstractString}=nothing,
)
    result = asymptotic_mk(sfs, divergence; model=model, cutoff=cutoff, cumulative=cumulative, n=n)
    return plot_asymptotic_mk(result; output_path=output_path)
end

end

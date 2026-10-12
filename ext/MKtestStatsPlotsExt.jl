module MKtestStatsPlotsExt

# plot_violins_empirical, loaded when StatsPlots is.

using MKtest, StatsPlots
import MKtest: plot_violins_empirical, empirical_test, _empirical_params, DataFrame, median

const _PLOT_EXTS = (".png", ".svg", ".pdf", ".ps", ".html")

# Row titles of Figure 4: Total/Strong/Weak.
function _row_label(c::Symbol)
    s = string(c)
    base = startswith(s, "ωₐ") ? "ωₐ" : startswith(s, "α") ? "α" : return s
    endswith(s, "_strong") && return "Strong " * base
    endswith(s, "_weak") && return "Weak " * base
    s == base && return "Total " * base
    return s
end

function plot_violins_empirical(
    tbl::DataFrame;
    omega::Bool = false,
    params::AbstractVector{Symbol} = _empirical_params(omega),
    alternative::Symbol = :greater,
    pvalues = nothing,
    output_path::Union{Nothing,AbstractString} = nothing,
    kwargs...,
)
    pv = pvalues === false ? nothing :
         isnothing(pvalues) ? empirical_test(tbl; params = params, alternative = alternative) :
         pvalues
    analyses = unique(tbl.analysis)
    k = length(analyses)
    is_case = tbl.type .== "case"

    rows = map(enumerate(params)) do (r, c)
        ctrl = [filter(isfinite, Float64.(tbl[(tbl.analysis .== a) .& .!is_case, c]))
                for a in analyses]
        case = [Float64(tbl[(tbl.analysis .== a) .& is_case, c][1]) for a in analyses]

        lo = min(minimum(minimum, ctrl), minimum(case))
        hi = max(maximum(maximum, ctrl), maximum(case))
        pad = 0.08 * max(hi - lo, eps())
        bottom = r == length(params)

        p = plot(;
            legend = false,
            ylabel = _row_label(c),
            xlims = k == 1 ? (0.0, 2.0) : (0.4, k + 0.6),
            ylims = (lo - 2.5pad, hi + pad),
            xticks = (1:k, bottom ? analyses : fill("", k)),
            xrotation = k > 1 ? 30 : 0,
            grid = :y,
            framestyle = :box,
        )
        for i in 1:k
            violin!(p, fill(i, length(ctrl[i])), ctrl[i];
                    color = :grey70, linecolor = :black, linewidth = 1, label = "")
            m = median(ctrl[i])
            plot!(p, [i - 0.35, i + 0.35], [m, m]; color = :black, lw = 3, label = "")
            scatter!(p, [i], [case[i]]; color = :black, markersize = 6,
                     markerstrokewidth = 0, label = "")
            if !isnothing(pv)
                s = pv[(pv.analysis .== analyses[i]) .& (pv.param .== string(c)), :stars]
                isempty(s) || isempty(s[1]) ||
                    annotate!(p, i, lo - 1.5pad, text(s[1], 14, :black))
            end
        end
        p
    end

    n = length(rows)
    defaults = (; layout = (n, 1), size = (max(360, 90k + 200), 250n),
                link = :x, left_margin = 5StatsPlots.mm, bottom_margin = 5StatsPlots.mm)
    fig = plot(rows...; merge(defaults, values(kwargs))...)

    if !isnothing(output_path)
        ext = lowercase(splitext(output_path)[2])
        ext in _PLOT_EXTS || throw(ArgumentError(
            "cannot infer a figure format from \"$output_path\". Use one of " *
            join(_PLOT_EXTS, ", ")))
        savefig(fig, output_path)
    end
    return fig
end

end

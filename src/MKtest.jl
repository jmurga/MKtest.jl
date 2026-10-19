module MKtest

using Parameters,
    SparseArrays,
    CSV,
    JLD2,
    DataFrames,
    StatsBase,
    Unzip,
    ThreadsX,
    Suppressor,
    CondaPkg,
    Random,
    Distributions,
    KernelDensity,
    GLM

# Analytical solutions
import NLsolve: nlsolve
import SpecialFunctions: polygamma, zeta, loggamma
import PoissonRandom: pois_rand
import LinearAlgebra.BLAS: get_num_threads, set_num_threads
import LinearAlgebra: mul!, svd

# Parse data
import GZip: open
import Parsers: parse
import OrderedCollections: OrderedDict
import Random: randstring
import Tables: table,matrix
import DelimitedFiles: readdlm

# MK-approaches
import LsqFit: curve_fit, vcov, dof
import MultipleTesting: BenjaminiHochberg, Bonferroni, adjust
import LinearAlgebra: dot
import Logging
import HypothesisTests: pvalue, FisherExactTest
import StatsBase: ordinalrank

# dfe_ml engine (src/ml/)
import Optim
import LineSearches
import ADTypes
import SpecialFunctions: gamma_inc, beta_inc, erfc, besselk

include("abcmk/parameters.jl")
include("abcmk/fixations.jl")
include("abcmk/sfs.jl")
include("abcmk/summary_statistics.jl")
include("abcmk/rates.jl")
include("data/parse.jl")
include("abcmk/abc.jl")
include("utils.jl")
include("ml/interface.jl")
include("ml/models.jl")
include("ml/integrators.jl")
include("ml/grids.jl")
include("ml/cache.jl")
include("ml/likelihood.jl")
include("ml/fit.jl")
include("ml/result.jl")
include("methods.jl")
include("bootstrap.jl")
include("alpha_b.jl")

# Plotting: the methods are in ext/MKtestPlotsExt.jl and load with `using Plots`.
"""
    plot_α_x(a_x::AbstractVector{<:Real}...; control, percentil, cutoff, n, labels, output_path, kwargs...)
    plot_α_x(curves::AbstractVector{<:AbstractVector{<:Real}}; ...)
    plot_α_x(α_table::AbstractMatrix; ...)
    plot_α_x(α_tables::AbstractVector{<:AbstractMatrix}; ...)
    plot_α_x(result::Union{NamedTuple,DataFrame}; ...)
    plot_α_x(a_x::AbstractVector{<:Real}, control::AbstractVector{<:AbstractVector{<:Real}}; ...)
    plot_α_x(sfs::Matrix, divergence::Vector; sample_size, cumulative, cutoff, n, ...)

Plot one or more alpha(x) curves against the derived allele count (DAC), optionally over a
percentile band of control curves. Requires `using Plots`.

# Arguments
 - `a_x`, `curves`: alpha(x) vectors; position i is DAC i (the x axis is the index).
 - `α_table`, `α_tables`: `[dac, frequency, alpha_x]` tables from [`α_x`](@ref); the DAC
   column is the x axis (of the first table for a vector of tables).
 - `result`: an [`asymptotic_mk`](@ref) result, a NamedTuple (one curve) or a DataFrame
   (one curve per row); its `f` gives the frequencies.
 - `sfs::Matrix`, `divergence::Vector`: compute alpha(x) with [`α_x`](@ref) first.
 - `control`: control curves (any of the curve inputs above), drawn as a grey band behind
   the curves; each must be at least as long as the longest curve. Non-finite values are
   skipped. Default `nothing`.
 - `percentil::AbstractVector{<:Real}`: low and high quantiles of the band. Default `[0.1, 0.9]`.
 - `cutoff::AbstractVector{<:Real}`: derived allele frequency window `[low, high]`; curves
   and band are drawn only inside it. Default `[0.0, 1.0]`. Frequencies come from the
   table's frequency column or the `asymptotic_mk` result's `f`; for plain vectors position
   i is DAC i of `n` haploid samples. With `sfs`/`divergence` the SFS is cut before the
   cumulative sum, as in `asymptotic_mk`.
 - `n::Union{Nothing,Int}`: haploid sample size behind plain vectors (default: longest
   curve + 1), or of a trimmed SFS (sfs/divergence only).
 - `labels`: legend labels, in curve order; missing ones are `α1`, `α2`, .... Default
   `nothing`. The vararg method also takes the older `label_1`, `label_2`.
 - `output_path`: save the figure; the format comes from the extension (`.png`, `.svg`,
   `.pdf`, `.ps`, `.html`; any other throws `ArgumentError`). Default `nothing`.
 - `sample_size::Union{Nothing,Int}`: project the SFS first (sfs/divergence only).
 - `cumulative::Bool`: use the cumulative SFS (sfs/divergence only). Default `true`.
 - Other keywords (`xlim`, `ylim`, `size`, `title`, ...) are passed to `Plots.plot` and
   override the defaults (`title = "αₓ values"`, `lw = 2`, `xlabel = "DAC"`, `ylabel = "α"`).

# Returns
 - `Plots.Plot`.
"""
function plot_α_x end

"""
    plot_asymptotic_mk(result::NamedTuple; output_path = nothing)
    plot_asymptotic_mk(sfs::Matrix, divergence::Vector; model = ExponentialFit(),
                       cutoff = [0.0, 1.0], cumulative = true, n = nothing, output_path = nothing)

Plot an asymptotic MK fit: the observed alpha(x) against the derived allele frequency, and,
when the asymptote is finite, the fitted curve of the model used ([`ExponentialFit`](@ref)
or [`LinearFit`](@ref), drawn over x in [0.05, 1]), the asymptote as a dashed line and its
95% CI as a band (if finite). Requires `using Plots`.

# Arguments
 - `result`: a single-dataset [`asymptotic_mk`](@ref) result (fields `f`, `alpha_x`,
   `asymptotic_mk`, `ci_low`, `ci_high`, `model`, `a`, `b`, `c`).
 - `sfs`, `divergence`, `model`, `cutoff`, `cumulative`, `n`: passed to
   [`asymptotic_mk`](@ref), whose result is then plotted.
 - `output_path`: save the figure; format from the extension (`.png`, `.svg`, `.pdf`,
   `.ps`, `.html`; any other throws `ArgumentError`). Default `nothing`.

# Returns
 - `Plots.Plot`.
"""
function plot_asymptotic_mk end

"""
    plot_violins_empirical(tbl; omega, params, alternative, pvalues, output_path, kwargs...)

Case against control sets per analysis, as Figure 4 of bioRxiv
doi:10.64898/2025.12.12.693836: one row per parameter (α, α_strong, α_weak by default), one column per
analysis. A grey violin shows the control sets, a black bar their median, a black dot the
case and stars its empirical p-value. Requires `using StatsPlots`.

# Arguments
 - `tbl::DataFrame`: output of [`case_control_table`](@ref); analyses keep their order.
 - `omega::Bool`: plot ωₐ, ωₐ_strong, ωₐ_weak instead. Default `false`.
 - `params`: parameter columns, one row each; overrides `omega`.
 - `alternative::Symbol`: passed to [`empirical_test`](@ref). Default `:greater`.
 - `pvalues`: an `empirical_test` result to annotate; `nothing` (default) computes it from
   `tbl` with `params` and `alternative`, `false` draws no stars.
 - `output_path`: save the figure; format from the extension (`.png`, `.svg`, `.pdf`,
   `.ps`, `.html`; any other throws `ArgumentError`). Default `nothing`.
 - Other keywords (`size`, `dpi`, ...) are passed to the final `plot` and override the
   defaults (layout one row per parameter, size scaled with the number of analyses).

# Returns
 - `Plots.Plot`.
"""
function plot_violins_empirical end

end

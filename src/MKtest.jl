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
import LinearAlgebra: mul!

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
include("methods.jl")
include("bootstrap.jl")
include("alpha_b.jl")

# Plotting: the methods are in ext/MKtestPlotsExt.jl and load with `using Plots`.
"""
    plot_α_x(a_x...; control, percentil, labels, output_path, kwargs...)
    plot_α_x(α_table; ...)
    plot_α_x(a_x, control; percentil, ...)
    plot_α_x(sfs, divergence; sample_size, cumulative, ...)
    plot_α_x(result; ...)

Plot one or more alpha(x) curves against the derived allele count, optionally over a
percentile band of control curves. Requires `using Plots`.

# Arguments
 - `a_x`: alpha(x) vectors, `[dac, frequency, alpha_x]` tables from `parse_sfs`/`α_x`
   (the DAC is then the x axis), or `asymptotic_mk` results (NamedTuple or DataFrame).
 - `sfs::Matrix`, `divergence::Vector`: compute alpha(x) with `α_x`.
 - `control`: control curves, drawn as a band. Default `nothing`.
 - `percentil::AbstractVector{<:Real}`: low and high quantiles of the band. Default `[0.1, 0.9]`.
 - `cutoff::AbstractVector{<:Real}`: derived allele frequency window `[low, high]`; curves
   and band are drawn only inside it. Default `[0.0, 1.0]`. Frequencies come from the
   table's frequency column or the `asymptotic_mk` result's `f`; for plain vectors position
   i is DAC i of `n` haploid samples. With `sfs`/`divergence` the SFS is cut before the
   cumulative sum, as in `asymptotic_mk`.
 - `n::Union{Nothing,Int}`: haploid sample size behind plain vectors (default: longest
   curve + 1), or of a trimmed SFS (sfs/divergence only).
 - `labels`: legend labels, in curve order. Default `nothing` (`α1`, `α2`, ...).
 - `output_path`: save the figure; the format comes from the extension (`.png`, `.svg`,
   `.pdf`, `.ps`, `.html`). Default `nothing`.
 - `sample_size::Union{Nothing,Int}`: project the SFS first (sfs/divergence only).
 - `cumulative::Bool`: use the cumulative SFS (sfs/divergence only). Default `true`.
 - Other keywords (`xlim`, `ylim`, `size`, `title`, ...) are passed to `Plots.plot`.

# Returns
 - `Plots.Plot`.
"""
function plot_α_x end

"""
    plot_asymptotic_mk(result; output_path)
    plot_asymptotic_mk(sfs, divergence; model, cutoff, cumulative, output_path)

Plot an asymptotic MK fit: the observed alpha(x), the fitted curve of the model used
(`ExponentialFit` or `LinearFit`), the asymptote and its 95% CI. Requires `using Plots`.

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
 - `pvalues`: an `empirical_test` result to annotate; computed from `tbl` by default, or
   `false` for no stars.
 - `output_path`: save the figure; format from the extension (`.png`, `.svg`, `.pdf`, ...).
 - Other keywords (`size`, `dpi`, ...) are passed to the final `plot`.

# Returns
 - `Plots.Plot`.
"""
function plot_violins_empirical end

end

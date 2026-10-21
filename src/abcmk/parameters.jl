################################
####### Define parameters ######
################################

"""
    parameters(; N = 1000, n = 500, dac = [2, 4, 5, 10, 20, 50, 200, 500, 700], ...)

Mutable keyword struct (`Parameters.@with_kw`) holding the ABC-MK analytical model. Every
model function ([`rates`](@ref), [`analytical_alpha`](@ref), [`sfs_neg`](@ref), ...)
reads these fields, and some set them. Set `N`, `n`, `dac` and `cutoff` to match the data; most
of the others are drawn from the priors by `rates`.

All selection coefficients (`γ_weak`, `γ_strong`, `γ_deleterious`, `γ_flanking`) are `2N` times the
heterozygote selection coefficient. `θ_*` and `ρ` are population-scaled per site
(`4Nμ`, `4Nr`).

# Fields
 - `N::Int64 = 1000`: population size of the model.
 - `n::Int64 = 500`: sample size in diploid individuals.
 - `dac::Vector{Int64} = [2, 4, 5, 10, 20, 50, 200, 500, 700]`: derived allele counts at
   which the α(x) tails are summarised. Ascending, `< nn`, within `cutoff`, and without
   singletons (`summary_statistics` rejects `1`).
 - `cutoff::Vector{Float64} = [0.0, 1.0]`: derived-allele frequency range of the data.
   Also part of the reference-table key `"N/n/cutoff"`.
 - `γ_weak::Int64 = 10`, `γ_strong::Int64 = 500`: weakly and strongly beneficial coefficients.
 - `α::Float64 = 0.4`: α, the fraction of nonsynonymous fixations that are
   beneficial (at `B = 0.999`).
 - `α_weak::Float64 = 0.2`: α_W, the part of `α` due to weakly beneficial fixations
   (`0 <= α_weak <= α`).
 - `shape::Float64 = 0.184`: shape of the gamma DFE of deleterious mutations.
 - `γ_deleterious::Int64 = -457`: mean of that DFE.
 - `scale::Float64 = shape / abs(γ_deleterious)`: rate of the gamma DFE. Computed at
   construction only: update it if you change `shape` or `γ_deleterious` afterwards.
 - `B::Float64 = 0.999`: background selection value (`Ne = B N`).
 - `γ_flanking::Int64 = -500`: selection on the linked deleterious flanking sites.
 - `L_flanking::Int64 = 200000`: length of the flanking region.
 - `ρ::Float64 = 0.001`: recombination rate.
 - `θ_flanking::Float64 = 1e-3`: mutation rate of the flanking sites; [`set_θ!`](@ref)
   sets it from `B`.
 - `θ_coding::Float64 = 1e-3`: mutation rate of the coding sites.
 - `p_weak::Float64 = 0`, `p_strong::Float64 = 0`: probabilities of weakly and strongly
   beneficial mutations; [`set_ppos!`](@ref) solves them from `α` and `α_weak`.
 - `f_syn::Float64 = 0.25`: fraction of new coding mutations that are
   synonymous (`1 - f_syn` are nonsynonymous).
 - `T_div::Float64 = 5.0`: divergence time used only by `analytical_alpha`, where it cancels
   out of α.
 - `isolines::Bool = false`: inbred lines; `parse_private_share` then takes `n`
   chromosomes instead of `2n`.
 - `NN::Int64 = 2N`, `nn::Int64 = 2n`: numbers of chromosomes in the population and in the
   sample.
"""
@with_kw mutable struct parameters
    γ_weak::Int64 = 10
    γ_strong::Int64 = 500
    α_weak::Float64 = 0.2
    α::Float64 = 0.4
    θ_flanking::Float64 = 1e-3
    γ_flanking::Int64 = -500
    θ_coding::Float64 = 1e-3
    shape::Float64 = 0.184
    γ_deleterious::Int64 = -457
    scale::Float64 = shape / abs(γ_deleterious)
    B::Float64 = 0.999
    p_weak::Float64 = 0
    p_strong::Float64 = 0
    N::Int64 = 1000
    n::Int64 = 500
    L_flanking::Int64 = 2 * 10^5
    ρ::Float64 = 0.001
    T_div::Float64 = 5.0
    cutoff::Vector{Float64} = [0.0, 1.0]
    isolines::Bool = false
    f_syn::Float64 = 0.25

    NN::Int64 = 2 * N
    nn::Int64 = 2 * n
    dac::Array{Int64,1} = [2, 4, 5, 10, 20, 50, 200, 500, 700]
end

"""
    copy_param(p::parameters)

Copy of the parameters struct, with its own `cutoff` and `dac` vectors. Faster than
`deepcopy`; `rates()` gives each task its own copy.
"""
function copy_param(p::parameters)
    parameters(
        γ_weak=p.γ_weak,
        γ_strong=p.γ_strong,
        α_weak=p.α_weak,
        α=p.α,
        θ_flanking=p.θ_flanking,
        γ_flanking=p.γ_flanking,
        θ_coding=p.θ_coding,
        shape=p.shape,
        γ_deleterious=p.γ_deleterious,
        scale=p.scale,
        B=p.B,
        p_weak=p.p_weak,
        p_strong=p.p_strong,
        N=p.N,
        n=p.n,
        L_flanking=p.L_flanking,
        ρ=p.ρ,
        T_div=p.T_div,
        cutoff=copy(p.cutoff),
        isolines=p.isolines,
        f_syn=p.f_syn,
        NN=p.NN,
        nn=p.nn,
        dac=copy(p.dac),
    )
end

function assertion_params(param::parameters)
    @assert (param.cutoff[1] >= 0.0) & (param.cutoff[end] <= 1.0) "Frequency cutoff must be at the interval [0,1]"
    @assert all(param.dac .< param.nn) "Selected DAC is larger than sample size"

    freq_dac = hcat(round.(collect(1:(param.nn-1)) ./ param.nn, digits=4), 1:(param.nn-1))
    freq_dac = freq_dac[
        (freq_dac[:, 1].>=param.cutoff[1]).&(freq_dac[:, 1].<=param.cutoff[end]),
        :,
    ]

    @assert all(in(freq_dac[:, 2]).(param.dac)) "Please select a DAC according to your frequency cutoff"

    # `summary_statistics` fills rows in `dac` order and `RjSampling` splits consecutive
    # tail sums into DAC intervals, so `dac` must be ascending.
    @assert issorted(param.dac) "param.dac must be in ascending order, got $(param.dac)"

end

# Residuals of the 2x2 system at a candidate root. Used by set_ppos! to decide whether a
# zero component is the answer or a failure.
function _ppos_residual(param::parameters, p_weak::Float64, p_strong::Float64)
    rt = alpha_exp_sim_tot(param, p_weak, p_strong) - param.α
    rl = alpha_exp_sim_low(param, p_weak, p_strong) - param.α_weak
    return (isfinite(rt) ? abs(rt) : Inf), (isfinite(rl) ? abs(rl) : Inf)
end

################################
###### Solving parameters ######
################################
"""
    Br(param, theta) -> Float64

Expected reduction in neutral diversity `π/π₀` caused by `param.L_flanking` linked deleterious
sites with mutation rate `theta` ([Charlesworth 1994](https://doi.org/10.1017/S0016672300032365)):

```math
\\frac{\\pi}{\\pi_{0}} = e^{-\\frac{4\\mu L_f}{2 r L_f + t}}
```

with `μ = theta/(2NN)`, `r = ρ/(2NN)` and `t = -γ_flanking/NN`.
"""
function Br(param::parameters, theta::Float64)
    ρ::Float64 = param.ρ
    t::Float64 = -1.0 * param.γ_flanking / (param.NN + 0.0)
    μ::Float64 = theta / (2.0 * param.NN)
    r::Float64 = ρ / (2.0 * param.NN)

    out::Float64 = ℯ^(-4 * μ * param.L_flanking / (2 * param.L_flanking * r + t))
    return out
end

# Set mutation rate given the expected reduction in nucleotide diversity (B value) in a locus.
"""
    set_θ!(param) -> Float64

Set `param.θ_flanking` to the flanking mutation rate that gives the background selection
value `param.B`, and return it.

`Br(θ) = exp(-2θ L_flanking / (L_flanking ρ - γ_flanking))` inverts in closed form:

```
θ_flanking = -log(B) (L_flanking ρ - γ_flanking) / (2 L_flanking)
```

Uses `B`, `L_flanking`, `ρ` and `γ_flanking`: set `ρ` before calling it.
"""
function set_θ!(param::parameters)
    # Br(θ) = exp(-4μLf / (2Lf·r + t)) with μ=θ/(2NN), r=ρ/(2NN), t=-γ_flanking/NN.
    # Numerator and denominator both carry a 1/NN factor that cancels, giving:
    #   Br(θ) = exp(-2θ·L_flanking / (L_flanking·ρ - γ_flanking))
    # Solving Br(θ) = B yields the formula below.
    param.θ_flanking = -log(param.B) * (param.L_flanking * param.ρ - param.γ_flanking) / (2.0 * param.L_flanking)
    return param.θ_flanking
end

function alpha_exp_sim_low(param::parameters, p_weak::Float64, p_strong::Float64)
    return fix_pos_sim(param, param.γ_weak, 0.5 * p_weak) / (
        fix_pos_sim(param, param.γ_weak, 0.5 * p_weak) +
        fix_pos_sim(param, param.γ_strong, 0.5 * p_strong) +
        fix_neg(param, 0.5 * p_weak + 0.5 * p_strong)
    )
end

function alpha_exp_sim_tot(param::parameters, p_weak::Float64, p_strong::Float64)
    return (
        fix_pos_sim(param, param.γ_weak, 0.5 * p_weak) +
        fix_pos_sim(param, param.γ_strong, 0.5 * p_strong)
    ) / (
        fix_pos_sim(param, param.γ_weak, 0.5 * p_weak) +
        fix_pos_sim(param, param.γ_strong, 0.5 * p_strong) +
        fix_neg(param, 0.5 * p_weak + 0.5 * p_strong)
    )
end

function solv_eqns(param, config)
    p_weak, p_strong = config
    return (
        alpha_exp_sim_tot(param, p_weak, p_strong) - param.α,
        alpha_exp_sim_low(param, p_weak, p_strong) - param.α_weak,
    )
end

"""
    set_ppos!(param; warm_start = [0.0, 0.0], atol = 1e-6, iterations = 50)

Solve `param.p_weak` and `param.p_strong`, the probabilities of weakly and strongly
beneficial mutations, so that the model's fixation α equals `param.α` and its weak
part equals `param.α_weak` at the current `B` and `θ_flanking`.

`NLsolve.nlsolve` runs from `warm_start`, then `[0.001, 0.001]`, `[0, 0.001]` and
`[0.001, 0]`, at most `iterations` steps each; negative roots are clamped to 0 and the
first root whose residuals are `<= atol` is kept. A zero component is accepted when it
satisfies the system (the right root when `α_weak == 0` or `α_weak == α`).

If no start reaches `atol`, the exact root of the system is tried: every fixation rate is
linear in `p_weak` and `p_strong`, so it has a closed form (`_ppos_linear`). nlsolve
fails where the deleterious DFE fixes almost nothing (shape above ~0.8 at the default
coefficients): the root is then `p ~ 1e-7` or smaller, far below every start. The closed
form runs only after nlsolve has failed, so a model nlsolve solves keeps its bits.
Converging starts take <= 12 steps at the default prior; the cap stops failing ones,
which would otherwise run nlsolve's default 1000 steps.

Returns the tuple `(p_weak, p_strong)`. Throws an `ErrorException` if neither route
reaches `atol` (a non-finite fixation rate, e.g. `shape == 1`, the pole of the Hurwitz
zeta); `rates` then drops the model.
"""
function set_ppos!(param::parameters; kw...)
    _set_ppos!(param; kw...)
    return param.p_weak, param.p_strong
end

# set_ppos!, returning the route that solved it: `:nlsolve` or `:linear`. `rates` reuses
# a root as the next draw's starting point only if nlsolve found it.
function _set_ppos!(param::parameters; warm_start::Vector{Float64} = [0.0, 0.0],
                    atol::Float64 = 1e-6, iterations::Int = 50)
    function f!(F, x, param=param)
        F[1] = alpha_exp_sim_tot(param, x[1], x[2]) - param.α
        F[2] = alpha_exp_sim_low(param, x[1], x[2]) - param.α_weak
    end

    best_l, best_h, best_err = 0.0, 0.0, Inf
    @inbounds for start in (warm_start, [0.001, 0.001], [0.0, 0.001], [0.001, 0.0])
        l, h = try
            nlsolve(f!, collect(float.(start)); iterations = iterations).zero
        catch
            continue
        end
        l = l < 0.0 ? 0.0 : l          # the system has no meaning at negative ppos
        h = h < 0.0 ? 0.0 : h
        rt, rl = _ppos_residual(param, l, h)
        err = max(rt, rl)
        if err < best_err
            best_l, best_h, best_err = l, h, err
        end
        err <= atol && break           # a zero component is fine if it satisfies the system
    end

    route = :nlsolve
    if !(best_err <= atol)
        l, h = _ppos_linear(param)
        rt, rl = _ppos_residual(param, l, h)
        err = max(rt, rl)
        if err < best_err
            best_l, best_h, best_err = l, h, err
            route = :linear
        end
    end

    best_err <= atol || throw(ErrorException(
        "set_ppos! did not reach α=$(param.α), α_weak=$(param.α_weak) " *
        "(best residual $(best_err)); model rejected"))

    param.p_weak, param.p_strong = best_l, best_h
    return route
end

"""
	_ppos_linear(param) -> (p_weak, p_strong)

Exact root of `set_ppos!`'s system. Every fixation rate is linear in ppos:
`F_weak = a·p_weak/2`, `F_strong = b·p_strong/2`, `F_neg = c·(1 - (p_weak + p_strong)/2)`
with `a = fix_pos_sim(γ_weak, 1)`, `b = fix_pos_sim(γ_strong, 1)`, `c = fix_neg(0)`.
Setting `F_weak = α_weak·D`, `F_strong = (α - α_weak)·D`, `F_neg = (1 - α)·D` gives

```
D = c / ((1 - α) + c·(α_weak/a + (α - α_weak)/b)),
p_weak = 2·α_weak·D/a,   p_strong = 2·(α - α_weak)·D/b.
```

Negative components (only if `α_weak > α`) are clipped to 0, as in `set_ppos!`; the
caller checks the residual.
"""
function _ppos_linear(param::parameters)
    a = fix_pos_sim(param, param.γ_weak, 1.0)
    b = fix_pos_sim(param, param.γ_strong, 1.0)
    c = fix_neg(param, 0.0)
    αt, αw = param.α, param.α_weak
    D = c / ((1.0 - αt) + c * (αw / a + (αt - αw) / b))
    return max(2.0 * αw * D / a, 0.0), max(2.0 * (αt - αw) * D / b, 0.0)
end

# Binomial matrices from binom_op, so they are built once per (N, n, B).
const _BINOM_CACHE = Dict{Tuple{Int,Int,UInt64},Dict{Float64,SparseMatrixCSC{Float64,Int64}}}()
const _BINOM_CACHE_SINGLE = Dict{Tuple{Int,Int,Float64},SparseMatrixCSC{Float64,Int64}}()

# binom * (1/j) for each (N, n, B). It does not depend on the model, so sfs_neut only
# scales it by B * θ_coding * f_syn.
const _NEUT_CACHE = Dict{Tuple{Int,Int,Float64},Vector{Float64}}()

"""
    TailWeights

Projection weights used by `rates()` (see `_sfs_dac!` in `src/abcmk/rates.jl`).

`rates()` only needs the cumulative SFS at the `dac` values, so the binomial projection
followed by `cumulative_vector!` reduces to one dot product per DAC against

```
W[k, j] = sum_{i >= dac[k]} binom[i, j]
```

`neut_tail[k]` is the same sum for the neutral SFS, which does not depend on the model.
This is faster than projecting the full spectrum and gives the same values.
"""
struct TailWeights
    W::Matrix{Float64}
    neut_tail::Vector{Float64}
end

# Keyed by (NN, nn, B, hash(dac)): W depends on the DACs.
const _TAIL_CACHE = Dict{Tuple{Int,Int,Float64,UInt64},TailWeights}()

"""
	clear_binom_cache!()

Clear cached binomial convolution matrices. Call this if you change N or n between runs
in the same Julia session.
"""
function clear_binom_cache!()
    empty!(_BINOM_CACHE)
    empty!(_BINOM_CACHE_SINGLE)
    empty!(_NEUT_CACHE)
    empty!(_TAIL_CACHE)
end

"""
    _tail_weights(binom, neut_base, dac) -> TailWeights

Build `TailWeights` from a binomial matrix: sum the rows of `binom` from the bottom up
and keep the running sum at each `dac` row.
"""
function _tail_weights(
    binom::SparseMatrixCSC{Float64,Int64},
    neut_base::Vector{Float64},
    dac::Vector{Int64},
)
    nrows, ncols = size(binom)
    @assert all(d -> 1 <= d <= nrows, dac) "DAC out of range for the binom matrix"

    W = zeros(Float64, length(dac), ncols)
    tail = zeros(Float64, ncols)
    neut_tail = zeros(Float64, length(dac))
    neut_run = 0.0

    # columns of `brow` are the rows of `binom`
    brow = sparse(transpose(binom))
    pos = Dict{Int,Int}(d => k for (k, d) in enumerate(dac))

    @inbounds for r in nrows:-1:1
        for p in brow.colptr[r]:(brow.colptr[r+1]-1)
            tail[brow.rowval[p]] += brow.nzval[p]
        end
        neut_run += neut_base[r]
        k = get(pos, r, 0)
        if k != 0
            copyto!(view(W, k, :), tail)
            neut_tail[k] = neut_run
        end
    end
    return TailWeights(W, neut_tail)
end

"""
    tail_weights_op(param, binom, B_bins) -> Dict{Float64,TailWeights}

One cached `TailWeights` per B value, from the matrices `binom_op` returned.
"""
function tail_weights_op(
    param::parameters,
    binom::Dict{Float64,SparseMatrixCSC{Float64,Int64}},
    B_bins::Vector{Float64},
)
    dac = param.dac
    dac_hash = hash(dac)
    out = Dict{Float64,TailWeights}()
    for b in B_bins
        A = binom[b]
        key = (param.NN, size(A, 1), b, dac_hash)
        out[b] = get!(_TAIL_CACHE, key) do
            neut_key = (param.NN, size(A, 1), b)
            neut_base = get(_NEUT_CACHE, neut_key, nothing)
            if neut_base === nothing
                # binom_op normally fills this; recompute if a caller skipped it.
                NN2 = size(A, 2) - 1
                h = Vector{Float64}(undef, NN2 + 1)
                @inbounds for j in 0:NN2
                    h[j+1] = j == 0 ? 0.0 : 1.0 / j
                end
                neut_base = A * h
                _NEUT_CACHE[neut_key] = neut_base
            end
            _tail_weights(A, neut_base, dac)
        end
    end
    return out
end

"""
    binom_op(param) -> SparseMatrixCSC{Float64,Int64}
    binom_op(param, B_bins) -> Dict{Float64,SparseMatrixCSC{Float64,Int64}}

Binomial sampling matrix that projects population frequencies onto the sample.

With `NN2 = ceil(NN B)` chromosomes (the population size reduced by background
selection), entry `[i, j+1]` is the probability of `i` derived alleles in `nn` sampled
chromosomes when the population frequency is `j/NN2`, for `i = 1:nn-1` and `j = 0:NN2`
(rounded to 10 digits, zeros dropped). The one-argument method uses `param.B`; the second
returns one matrix per value of `B_bins`. Results are cached per `(NN, nn, B)`; call
`clear_binom_cache!()` to free them.
"""
function binom_op(param::parameters, B_bins::Vector{Float64})
    key = (param.NN, param.nn, hash(B_bins))
    haskey(_BINOM_CACHE, key) && return _BINOM_CACHE[key]

    bn = Dict(
        B_bins[i] => spzeros(param.nn + 1, param.NN) for i = 1:length(B_bins)
    )

    for b in B_bins
        NN2 = convert(Int64, ceil(param.NN * b))
        samples = collect(1:(param.nn-1))
        p_size = collect(0:NN2)
        samples_freqs = permutedims(p_size / NN2)
        neutral_sfs = @. 1 / p_size
        replace!(neutral_sfs, Inf => 0.0)

        f(x, y=param.nn) = Binomial(y, x)
        z = f.(samples_freqs)

        out = pdf.(z, samples)
        out = round.(out, digits=10)
        out_sparse = SparseArrays.dropzeros(SparseArrays.sparse(out))
        bn[b] = out_sparse

        # binom * (1/j), reused by sfs_neut for every model.
        neut_key = (param.NN, size(out_sparse, 1), b)
        if !haskey(_NEUT_CACHE, neut_key)
            _NEUT_CACHE[neut_key] = out_sparse * neutral_sfs
        end
    end

    _BINOM_CACHE[key] = bn
    return (bn)
end

function binom_op(param)
    @unpack NN, nn, B = param
    key = (NN, nn, B)
    haskey(_BINOM_CACHE_SINGLE, key) && return _BINOM_CACHE_SINGLE[key]

    NN2 = convert(Int64, ceil(NN * B))
    samples = collect(1:(nn-1))
    p_size = collect(0:NN2)
    samples_freqs = permutedims(p_size / NN2)
    neutral_sfs = @. 1 / p_size
    replace!(neutral_sfs, Inf => 0.0)

    f(x, y=nn) = Binomial(y, x)
    z = f.(samples_freqs)

    out = pdf.(z, samples)
    out = round.(out, digits=10)
    out_sparse = SparseArrays.dropzeros(SparseArrays.sparse(out))

    _BINOM_CACHE_SINGLE[key] = out_sparse

    # binom * (1/j), reused by sfs_neut.
    neut_key = (NN, size(out_sparse, 1), B)
    if !haskey(_NEUT_CACHE, neut_key)
        _NEUT_CACHE[neut_key] = out_sparse * neutral_sfs
    end

    return (out_sparse)
end

"""
    Φ(param, s::Int64) -> Float64

Factor (in `(0, 1]`) by which background selection and linkage reduce the fixation
probability of a beneficial mutation with coefficient `s` (`2N` times the heterozygote
coefficient), from the `param.L_flanking` linked deleterious sites
([Charlesworth 1994](https://doi.org/10.1017/S0016672300032365),
[Barton 1995](https://www.genetics.org/content/140/2/821)).

The product over linked sites of the per-site reduction gives

```math
\\Phi = e^{-\\frac{2 t \\mu \\left(\\psi_1\\left(\\frac{r + 2s + t}{r}\\right) - \\psi_1\\left(\\frac{r (L_f + 1) + 2s + t}{r}\\right)\\right)}{r^2}}
```

with `t = |γ_flanking|/NN`, `μ = θ_flanking/(2NN)`, `r = ρ/(2NN)`, `s = s/NN` and
`ψ₁` the trigamma function: Uricchio et al. 2019, Supplementary eq. 14, and
Murga-Moreno et al. 2024, eq. 9.
"""
function Φ(param::parameters, s::Int64)
    S::Float64 = abs(param.γ_flanking / (1.0 * param.NN))
    r::Float64 = param.ρ / (2.0 * param.NN)
    μ::Float64 = param.θ_flanking / (2.0 * param.NN)
    s::Float64 = s / (param.NN * 1.0)

    Ψ0::Float64 = polygamma(1, (r + 2.0 * s + S) / r)
    Ψ1::Float64 = polygamma(1, (r * (param.L_flanking + 1) + 2.0 * s + S) / r)

    out::Float64 = (ℯ^(-2.0 * S * μ * (Ψ0 - Ψ1) / (r^2)))
    return out
end

"""
    analytical_alpha(param; cumulative = false) -> (α, α_nopos)

Expected α(x) of the model in `param`, at every derived allele count `1:nn-1`:

```math
\\mathbb{E}[\\alpha_{x}] = 1 - \\frac{\\mathbb{E}[D_{s}]}{\\mathbb{E}[D_{n}]}\\frac{\\mathbb{E}[P_{n}(x)]}{\\mathbb{E}[P_{s}(x)]}
```

`p_weak` and `p_strong` are solved at `B = 0.999` ([`set_ppos!`](@ref)), then the fixation
rates and spectra are computed at `param.B`. `param` is modified: `p_weak`, `p_strong` and
`θ_flanking` (the value for `param.B`) are left set.

# Keywords
 - `cumulative::Bool = false`: use the tails (DAC >= x) of the spectra instead of the
   per-frequency counts.

# Returns
 - `α::Vector{Float64}`: α(x) with all nonsynonymous polymorphism (deleterious and
   beneficial).
 - `α_nopos::Vector{Float64}`: α(x) with deleterious polymorphism only, whose asymptote
   is the true α.
"""
function analytical_alpha(param::parameters; cumulative::Bool=false)
    binom = binom_op(param)
    ################################################################
    # Solve the model similarly to original python mktest  #
    ################################################################
    B = param.B

    set_θ!(param)
    θ_flanking = param.θ_flanking
    # Solve the probabilities of fixations without BGS
    ## Set non-bgs
    param.B = 0.999
    ## Solve the mutation rate
    set_θ!(param)
    ## Solve the probabilities
    set_ppos!(param)
    # Return to the original values
    param.θ_flanking = θ_flanking
    param.B = B

    ################################################################
    # Accounting for positive alleles segregating due to linkage #
    ################################################################

    # Fixation
    fN = param.B * fix_neut(param)
    fNeg = param.B * fix_neg(param, 0.5 * param.p_strong + 0.5 * param.p_weak)
    fPosL = fix_pos_sim(param, param.γ_weak, 0.5 * param.p_weak)
    fPosH = fix_pos_sim(param, param.γ_strong, 0.5 * param.p_strong)

    ds = fN
    dn = fNeg + fPosL + fPosH

    ## Polymorphism
    neut = sfs_neut(param, binom)

    selH = sfs_pos(param, param.γ_strong, param.p_strong, binom)

    selL = sfs_pos(param, param.γ_weak, param.p_weak, binom)
    selN = sfs_neg(param, param.p_strong + param.p_weak, binom)

    function split_columns(matrix::Array{Float64,2})
        (view(matrix, :, i) for i = 1:size(matrix, 2))
    end

    tmp = hcat(neut, selH, selL, selN)
    if cumulative
        tmp = cumulative_sfs(tmp, false)
    end

    neut, selH, selL, selN = split_columns(tmp)
    sel = (selH + selL) + selN

    # ps = @. neut / (sel+neut)
    # pn = @. sel / (sel+neut)

    ## Outputs
    α = @. 1 - ((ds / dn) * (sel / neut))

    #############################################################
    #		Accounting for neutral and deleterious alleles		#	#############################################################
    ## Fixation
    fN_nopos = fN * (param.θ_coding / 2.0) * param.T_div * param.NN
    fNeg_nopos = fNeg * (param.θ_coding / 2.0) * param.T_div * param.NN
    fPosL_nopos = fPosL * (param.θ_coding / 2.0) * param.T_div * param.NN
    fPosH_nopos = fPosH * (param.θ_coding / 2.0) * param.T_div * param.NN

    ds_nopos = fN_nopos
    dn_nopos = fNeg_nopos + fPosL_nopos + fPosH_nopos

    ## Polymorphism
    sel_nopos = selN
    ps_nopos = @. neut / (sel_nopos + neut)
    pn_nopos = @. sel_nopos / (sel_nopos + neut)

    α_nopos = 1 .- (ds_nopos / dn_nopos) .* (sel_nopos ./ neut)

    ##########
    # Output #
    ##########
    return (α, α_nopos)
end

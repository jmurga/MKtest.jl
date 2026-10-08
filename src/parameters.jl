################################
####### Define parameters ######
################################

"""
Mutable structure containing the variables required to solve the analytical approach. All the functions are solve using the internal values of the structure. You should declare a mutable structure to the perform the analytical estimations.

# Parameters
 - `gam_flanking::Int64`: Selection coefficient for deleterious alleles
 - `gL::Int64`: Selection coefficient for weakly benefical alleles
 - `gH::Int64`: Selection coefficient for strongly benefical alleles
 - `al_low::Float64`: Proportion of α due to weak selection
 - `al_tot::Float64`: α
 - `θ_flanking::Float64`: Mutation rate defining BGS strength
 - `θ_coding::Float64`: Mutation rate on coding region
 - `al::Float64`: DFE shape parameter
 - `be::Float64`: DFE scale parameter
 - `B::Float64`: BGS strength
 - `B_bins::Array{Float64,1}`: BGS values to simulate
 - `ppos_l::Float64`: Fixation probabily of weakly beneficial alleles
 - `ppos_h::Float64`: Fixation probabily of strongly beneficial alleles
 - `N::Int64`: Population size
 - `n::Int64`: Sample size
 - `Lf::Int64`: Flanking region length
 - `ρ::Float64`: Recombination rate
 - `TE::Float64`
"""
@with_kw mutable struct parameters
    gL::Int64 = 10
    gH::Int64 = 500
    al_low::Float64 = 0.2
    al_tot::Float64 = 0.4
    θ_flanking::Float64 = 1e-3
    gam_flanking::Int64 = -500
    θ_coding::Float64 = 1e-3
    shape::Float64 = 0.184
    gam_dfe::Int64 = -457
    scale::Float64 = shape / abs(gam_dfe)
    B::Float64 = 0.999
    ppos_l::Float64 = 0
    ppos_h::Float64 = 0
    N::Int64 = 1000
    n::Int64 = 500
    Lf::Int64 = 2 * 10^5
    ρ::Float64 = 0.001
    TE::Float64 = 5.0
    diploid::Bool = false
    cutoff::Vector{Float64} = [0.0, 1.0]
    isolines::Bool = false
    syn_nonsyn_ratio::Float64 = 0.25

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
        gL=p.gL,
        gH=p.gH,
        al_low=p.al_low,
        al_tot=p.al_tot,
        θ_flanking=p.θ_flanking,
        gam_flanking=p.gam_flanking,
        θ_coding=p.θ_coding,
        shape=p.shape,
        gam_dfe=p.gam_dfe,
        scale=p.scale,
        B=p.B,
        ppos_l=p.ppos_l,
        ppos_h=p.ppos_h,
        N=p.N,
        n=p.n,
        Lf=p.Lf,
        ρ=p.ρ,
        TE=p.TE,
        diploid=p.diploid,
        cutoff=copy(p.cutoff),
        isolines=p.isolines,
        syn_nonsyn_ratio=p.syn_nonsyn_ratio,
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
function _ppos_residual(param::parameters, ppos_l::Float64, ppos_h::Float64)
    rt = alpha_exp_sim_tot(param, ppos_l, ppos_h) - param.al_tot
    rl = alpha_exp_sim_low(param, ppos_l, ppos_h) - param.al_low
    return (isfinite(rt) ? abs(rt) : Inf), (isfinite(rl) ? abs(rl) : Inf)
end

################################
###### Solving parameters ######
################################
"""
	Br(Lmax,theta)

Expected reduction in nucleotide diversity. Explored at [Charlesworth B., 1994](https://doi.org/10.1017/S0016672300032365):

```math
\\frac{\\pi}{\\pi_{0}} = e^{\\frac{4\\mu L}{2rL+t}}
```

# Arguments
 - `param::parameters`
 - `theta::Float64`
# Returns
 - `Float64`: expected reduction in diversity given a non-coding length, mutation rate and defined recombination.
"""
function Br(param::parameters, theta::Float64)
    ρ::Float64 = param.ρ
    t::Float64 = -1.0 * param.gam_flanking / (param.NN + 0.0)
    μ::Float64 = theta / (2.0 * param.NN)
    r::Float64 = ρ / (2.0 * param.NN)

    out::Float64 = ℯ^(-4 * μ * param.Lf / (2 * param.Lf * r + t))
    return out
end

# Set mutation rate given the expected reduction in nucleotide diversity (B value) in a locus.
"""
	set_θ!(param)

Find the mutation rate `θ_flanking` such that `Br(param, θ_flanking) == param.B`.

`Br(θ) = exp(-2·θ·Lf / (Lf·ρ - gam_flanking))` inverts analytically:

```
θ_flanking = -log(B) · (Lf·ρ - gam_flanking) / (2·Lf)
```

The closed form gives the same value as the previous root search (`Roots.find_zero`)
to float precision, without allocations, which makes `rates()` faster.

# Returns
 - mutates `param.θ_flanking`.
"""
function set_θ!(param::parameters)
    # Br(θ) = exp(-4μLf / (2Lf·r + t)) with μ=θ/(2NN), r=ρ/(2NN), t=-gam_flanking/NN.
    # Numerator and denominator both carry a 1/NN factor that cancels, giving:
    #   Br(θ) = exp(-2θ·Lf / (Lf·ρ - gam_flanking))
    # Solving Br(θ) = B yields the formula below.
    param.θ_flanking = -log(param.B) * (param.Lf * param.ρ - param.gam_flanking) / (2.0 * param.Lf)
    return param.θ_flanking
end

function alpha_exp_sim_low(param::parameters, ppos_l::Float64, ppos_h::Float64)
    return fix_pos_sim(param, param.gL, 0.5 * ppos_l) / (
        fix_pos_sim(param, param.gL, 0.5 * ppos_l) +
        fix_pos_sim(param, param.gH, 0.5 * ppos_h) +
        fix_neg(param, 0.5 * ppos_l + 0.5 * ppos_h)
    )
end

function alpha_exp_sim_tot(param::parameters, ppos_l::Float64, ppos_h::Float64)
    return (
        fix_pos_sim(param, param.gL, 0.5 * ppos_l) +
        fix_pos_sim(param, param.gH, 0.5 * ppos_h)
    ) / (
        fix_pos_sim(param, param.gL, 0.5 * ppos_l) +
        fix_pos_sim(param, param.gH, 0.5 * ppos_h) +
        fix_neg(param, 0.5 * ppos_l + 0.5 * ppos_h)
    )
end

function solv_eqns(param, config)
    ppos_l, ppos_h = config
    return (
        alpha_exp_sim_tot(param, ppos_l, ppos_h) - param.al_tot,
        alpha_exp_sim_low(param, ppos_l, ppos_h) - param.al_low,
    )
end

"""
	set_ppos!(param; warm_start = [0.0, 0.0], atol = 1e-6)

Solve `(ppos_l, ppos_h)`, the probabilities of weakly and strongly beneficial mutations,
so the model's FIXATION alpha equals `param.al_tot` and its weak share `param.al_low`.

Tries `warm_start`, then `[0.001, 0.001]`, `[0, 0.001]`, `[0.001, 0]`, and accepts the
first root whose residuals are `<= atol`. A zero component is accepted when it satisfies
the system (it is the correct root when `al_low == 0` or `al_low == al_tot`). If no start
reaches the tolerance it throws, and `simulate_models!` drops the model instead of
writing a row whose alpha label does not match its fixations.
"""
function set_ppos!(param::parameters; warm_start::Vector{Float64} = [0.0, 0.0],
                   atol::Float64 = 1e-6)
    function f!(F, x, param=param)
        F[1] = alpha_exp_sim_tot(param, x[1], x[2]) - param.al_tot
        F[2] = alpha_exp_sim_low(param, x[1], x[2]) - param.al_low
    end

    best_l, best_h, best_err = 0.0, 0.0, Inf
    @inbounds for start in (warm_start, [0.001, 0.001], [0.0, 0.001], [0.001, 0.0])
        l, h = try
            nlsolve(f!, collect(float.(start))).zero
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

    best_err <= atol || throw(ErrorException(
        "set_ppos! did not reach al_tot=$(param.al_tot), al_low=$(param.al_low) " *
        "(best residual $(best_err)); model rejected"))

    param.ppos_l, param.ppos_h = best_l, best_h
end

# Binomial matrices from binom_op, so they are built once per (N, n, B).
const _BINOM_CACHE = Dict{Tuple{Int,Int,UInt64},Dict{Float64,SparseMatrixCSC{Float64,Int64}}}()
const _BINOM_CACHE_SINGLE = Dict{Tuple{Int,Int,Float64},SparseMatrixCSC{Float64,Int64}}()

# binom * (1/j) for each (N, n, B). It does not depend on the model, so sfs_neut only
# scales it by B * θ_coding * syn_nonsyn_ratio.
const _NEUT_CACHE = Dict{Tuple{Int,Int,Float64},Vector{Float64}}()

"""
    TailWeights

Projection weights used by `rates()` (see `_sfs_dac!` in `src/rates.jl`).

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
	binom_op(param)

Binomial convolution to sample the allele frequencies probabilites depending on background selection values and sample size.

# Arguments
 - `param::parameters`
 - `B_bins::Vector{Float64}`

# Returns
 - `Array{Float64,2}`: convoluted SFS for each B value defined in the model (B_bins). The estimations are saved at *convolutedBn.bn*.
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

	Φ(param,s)


Reduction in fixation probabilty due to background selection and linkage. The formulas used have been subjected to several theoretical works ([Charlesworth B., 1994](https://doi.org/10.1017/S0016672300032365), [Hudson et al., 1995](https://www.genetics.org/content/141/4/1605), [Nordborg et al. 1995](https://doi.org/10.1017/S0016672300033619), [Barton NH., 1995](https://www.genetics.org/content/140/2/821)).

The fixation probabilty of selected alleles are reduce by a factor ``\\phi``:

```math
\\phi(t,s) = e^{[\\frac{-2\\mu}{t(1+\\frac{rL}{t}+\\frac{2s}{t})}]}
```

Multiplying across all deleterious linkes sites, we find:

```math
\\Phi = \\prod_{1}^{L} = \\phi(t,s) = e^{\\frac{-2t\\mu(\\psi[1,\\frac{r+2s+L}{r}] - \\psi[1,\\frac{r(L+1)+2s+t}{r}])}{r^2}}
```

```math
\\phi(t,s) = e^{\\frac{-2t\\mu(\\psi[1,\\frac{r+2s+L}{r}] - \\psi[1,\\frac{r(L+1)+2s+t}{r}])}{r^2}}
```

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model
 - `gamma::Int64`: selection coefficient.

The expression is Uricchio et al. 2019 Supplementary eqn. 14 / Murga-Moreno et al. 2024
(G3) eqn. 9. Previous versions omitted the factor of two on the beneficial coefficient
(and the `r +` in the first polygamma), which made the reduction of strongly beneficial
fixations about twice too large.

# Returns
 - `Float64`: expected rate of positive fixations under background selection.

"""
function Φ(param::parameters, s::Int64)
    S::Float64 = abs(param.gam_flanking / (1.0 * param.NN))
    r::Float64 = param.ρ / (2.0 * param.NN)
    μ::Float64 = param.θ_flanking / (2.0 * param.NN)
    s::Float64 = s / (param.NN * 1.0)

    Ψ0::Float64 = polygamma(1, (r + 2.0 * s + S) / r)
    Ψ1::Float64 = polygamma(1, (r * (param.Lf + 1) + 2.0 * s + S) / r)

    out::Float64 = (ℯ^(-2.0 * S * μ * (Ψ0 - Ψ1) / (r^2)))
    return out
end

"""
	analytical_alpha(param)

Analytical α(x) estimation. Solve α(x) generally. We used the expected fixation rates and SFSto approach the asympotic value accouting for background selection, weakly and strong positive selection. α(x) can be estimated taking into account the role of positive selected alleles or not. In this way we explore the role of linkage to deleterious alleles in the coding region.

```math
\\mathbb{E}[\\alpha_{x}] =  1 - \\left(\\frac{\\mathbb{E}[D_{s}]}{\\mathbb{E}[D_{N}]}\\frac{\\mathbb{E}[P_{N}]}{\\mathbb{E}[P_{S}]}\\right)
```

# Arguments
 - `param::parameters`
# Returns
 - `Tuple{Vector{Float64},Vector{Float64}}`: α(x) accounting for weak adaptation, α(x) non-accouting for weak adaptation.
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
    fNeg = param.B * fix_neg(param, 0.5 * param.ppos_h + 0.5 * param.ppos_l)
    fPosL = fix_pos_sim(param, param.gL, 0.5 * param.ppos_l)
    fPosH = fix_pos_sim(param, param.gH, 0.5 * param.ppos_h)

    ds = fN
    dn = fNeg + fPosL + fPosH

    ## Polymorphism
    neut = sfs_neut(param, binom)

    selH = sfs_pos(param, param.gH, param.ppos_h, binom)

    selL = sfs_pos(param, param.gL, param.ppos_l, binom)
    selN = sfs_neg(param, param.ppos_h + param.ppos_l, binom)

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
    fN_nopos = fN * (param.θ_coding / 2.0) * param.TE * param.NN
    fNeg_nopos = fNeg * (param.θ_coding / 2.0) * param.TE * param.NN
    fPosL_nopos = fPosL * (param.θ_coding / 2.0) * param.TE * param.NN
    fPosH_nopos = fPosH * (param.θ_coding / 2.0) * param.TE * param.NN

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

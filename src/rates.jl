################################################################################
# SFS projection: tail weights                                                 #
#                                                                              #
# The per-DAC cumulative sample SFS is projected straight from the population- #
# frequency source vector with precomputed tail-summed binomial weights        #
# `W[k,j] = sum_{i >= dac[k]} binom[i,j]` (`TailWeights`, parameters.jl). It is #
# the full binomial projection + `cumulative_vector!` + indexing at `dac`,     #
# reassociated: agreement 3-4e-15 (dev docs/abcmk_performance.md), ~3x faster.#
################################################################################

"""
    _tail_ops(param, binom, B_bins) -> Vector{TailWeights}

One `TailWeights` per B bin, positionally aligned to `B_bins` (no Float64 hash in the
inner loop).
"""
function _tail_ops(param::parameters,
                   binom::Dict{Float64,SparseMatrixCSC{Float64,Int64}},
                   B_bins::Vector{Float64})
    tw = tail_weights_op(param, binom, B_bins)
    return TailWeights[tw[b] for b in B_bins]
end

################################################################################
# Recombination-rate prior                                                     #
################################################################################

"""
    _rho_prior(rho, param) -> (u_rho, rho_fixed)

Resolve the `rho` keyword of `rates()` into either a `Uniform`
distribution to draw from (one draw per model) or a fixed scalar.

| `rho`                  | result                                        |
|------------------------|-----------------------------------------------|
| `nothing` (default)    | `(nothing, param.ρ)` -- today's behaviour     |
| `0.002::Real`          | `(nothing, 0.002)`                            |
| `[1e-4, 5e-3]::Vector` | `(Uniform(1e-4, 5e-3), NaN)`                  |
| `[0.001]` (length 1)   | `(nothing, 0.001)`                            |

`rho[1]` / `rho[end]` are used for the bounds, mirroring how the `gH`, `gL`,
`gam_flanking` and `gam_dfe` ranges are read.

When the first element is `nothing` the caller must NOT consume any RNG, so that
a seeded `rates(...; rho = nothing)` reproduces the pre-`rho` random stream
bit-for-bit.
"""
_rho_prior(::Nothing, param::parameters) = (nothing, param.ρ)

function _rho_prior(rho::Real, param::parameters)
    rho > 0 || throw(ArgumentError("rho must be positive, got $rho"))
    return (nothing, Float64(rho))
end

function _rho_prior(rho::AbstractVector{<:Real}, param::parameters)
    isempty(rho) && throw(ArgumentError("rho range is empty"))
    all(>(0), rho) ||
        throw(ArgumentError("rho must be positive, got $rho"))
    lo, hi = Float64(rho[1]), Float64(rho[end])
    lo <= hi ||
        throw(ArgumentError("rho range must be ordered [lo, hi], got [$lo, $hi]"))
    lo == hi && return (nothing, lo)
    return (Uniform(lo, hi), NaN)
end

################################################################################
# Provenance and storage                                                       #
################################################################################

# What `_read_rates` reports when a file predates option provenance.
const _NO_PROVENANCE = "unrecorded (file predates option provenance)"

"""
    _rates_options() -> String

One-line record of the forward model a rates table was built under, stored in the file
next to the table. A table without it was built by an earlier MKtest (a different model).
"""
_rates_options() = "model=best-2026-09-28 (Phi published (eq 9); beneficial SFS at B level with sojourn at gamma*B; Phi flat on the SFS; guarded ppos solve); projection=TailWeight"

_jld2_group_key(N::Int, n::Int, cutoff_lo, cutoff_hi) =
    string(N) * "/" * string(n) * "/cutoff=[" * string(cutoff_lo) * "," * string(cutoff_hi) * "]"

"""
    _read_rates(input_file, param) -> NamedTuple

`(; models, neut, sel, dsdn, dac, options)` of the `(N, n, cutoff)` block of a rates file.
Reads tables written by this and earlier versions; `options` says which model built them.
"""
function _read_rates(input_file::String, param::parameters)
    grp = _jld2_group_key(param.N, param.n, param.cutoff[1], param.cutoff[end])
    return jldopen(input_file, "r") do h
        tmp = h[grp]
        (;
            models = tmp["models"],
            neut   = tmp["neut"],
            sel    = tmp["sel"],
            dsdn   = tmp["dsdn"],
            dac    = tmp["dac"],
            options = ("options" in keys(tmp)) ? tmp["options"] : _NO_PROVENANCE,
        )
    end
end

################################################################################
# Buffered hot path (S-tier #1 + #2)                                           #
#                                                                              #
# `SolveRatesBuffers` holds the four SFS output vectors + a scratch buffer     #
# (sized for the largest NN2 across B_bins). `solve_rates!` reuses them across #
# all 40 B-bins per model and writes results directly into preallocated row    #
# views of the global result matrices (`m_view`, `ps_view`, `pn_view`,         #
# `f_view`). This eliminates per-iteration SFS vector allocations, the         #
# `cumulative_vector!` intermediates, and the `vcat(m...)` reduction at the    #
# end of `rates()`.                                                            #
#                                                                              #
################################################################################

"""
    SolveRatesBuffers(binom, n_dac)

Preallocated work buffers for the ABC forward model. Sized once per task; reused across
the whole B-bin loop and every iteration assigned to that task. `scratch` holds the
population-frequency source vector (sized for the largest `NN2` in `binom`); `acc_tmp`
receives each `mul!(., W, src)` and `acc_sel` sums the selected components.
"""
struct SolveRatesBuffers
    scratch::Vector{Float64}
    acc_sel::Vector{Float64}
    acc_tmp::Vector{Float64}
    # Warm-start cache for `set_ppos!` -- `prev_ppos = [ppos_l, ppos_h]` from the
    # previous successful iteration in this task's chunk. nlsolve converges
    # noticeably faster (~2-3 Newton steps vs ~10 from cold start) when the
    # previous prior draw is similar to the next.
    prev_ppos::Vector{Float64}
end

function SolveRatesBuffers(binom::Dict{Float64,SparseMatrixCSC{Float64,Int64}},
                           n_dac::Int)
    scratch_len = maximum(size(b, 2) for b in values(binom))
    SolveRatesBuffers(
        Vector{Float64}(undef, scratch_len),
        Vector{Float64}(undef, n_dac),
        Vector{Float64}(undef, n_dac),
        zeros(Float64, 2),
    )
end

"""
    _sfs_dac!(ps_row, pn_row, buf, param, op::TailWeights)

Write the `n_dac` neutral and selected cumulative SFS values for one B bin.
"""
function _sfs_dac!(ps_row::AbstractVector, pn_row::AbstractVector,
                   buf::SolveRatesBuffers, param::parameters, op::TailWeights)
    @unpack ppos_h, ppos_l, gL, gH, dac = param
    has_h = ppos_h > 0.0
    has_l = ppos_l > 0.0
    W = op.W
    ncols = size(W, 2)
    nd = length(dac)

    # Neutral: no mat-vec, the tail sum is model-independent and precomputed.
    nsc = param.B * param.θ_coding * param.syn_nonsyn_ratio
    @inbounds @simd for k in 1:nd
        ps_row[k] = nsc * op.neut_tail[k]
    end

    # Selected, summed as (sel_h + sel_l) + sel_neg. Beneficial classes: B level, sojourn
    # at γ·B, Φ as one scalar (see `sfs_pos`).
    acc = buf.acc_sel
    fill!(acc, 0.0)
    if has_h
        src = _src_pos!(buf.scratch, param, gH * param.B, ppos_h, ncols)
        mul!(buf.acc_tmp, W, src)
        fh = param.B * param.θ_coding * Φ(param, gH) * (1 - param.syn_nonsyn_ratio)
        @inbounds @simd for k in 1:nd
            acc[k] += fh * buf.acc_tmp[k]
        end
    end
    if has_l
        src = _src_pos!(buf.scratch, param, gL * param.B, ppos_l, ncols)
        mul!(buf.acc_tmp, W, src)
        fl = param.B * param.θ_coding * Φ(param, gL) * (1 - param.syn_nonsyn_ratio)
        @inbounds @simd for k in 1:nd
            acc[k] += fl * buf.acc_tmp[k]
        end
    end
    let src = _src_neg!(buf.scratch, param, ppos_l + ppos_h, ncols)
        mul!(buf.acc_tmp, W, src)
        fneg = param.B * param.θ_coding * (1 - param.syn_nonsyn_ratio)
        @inbounds @simd for k in 1:nd
            acc[k] += fneg * buf.acc_tmp[k]
        end
    end

    @inbounds @simd for k in 1:nd
        pn_row[k] = acc[k]
    end
    return nothing
end

"""
    solve_rates!(m_row, ps_row, pn_row, f_row, buf, param, op)

Writes the 10-element model row, the per-DAC neutral & selected rows, and the 4-element
fixation row for one B bin directly into the caller-provided views. No allocations on
the hot path (buffers are reused).
"""
function solve_rates!(
    m_row::AbstractVector,
    ps_row::AbstractVector,
    pn_row::AbstractVector,
    f_row::AbstractVector,
    buf::SolveRatesBuffers,
    param::parameters,
    op::TailWeights,
)
    @unpack B, ppos_h, ppos_l, gL, gH, al_low, al_tot,
            gam_flanking, shape, scale, dac = param

    has_h = ppos_h > 0.0
    has_l = ppos_l > 0.0

    # Fixation
    f_n = B * fix_neut(param)
    f_neg = B * fix_neg(param, 0.5 * ppos_h + 0.5 * ppos_l)
    f_pos_l = has_l ? fix_pos_sim(param, gL, 0.5 * ppos_l) : 0.0
    f_pos_h = has_h ? fix_pos_sim(param, gH, 0.5 * ppos_h) : 0.0

    _sfs_dac!(ps_row, pn_row, buf, param, op)

    # Write directly into preallocated output row views
    @inbounds m_row[1] = B * 1000
    @inbounds m_row[2] = al_low
    @inbounds m_row[3] = al_tot
    @inbounds m_row[4] = gam_flanking
    @inbounds m_row[5] = gL
    @inbounds m_row[6] = gH
    @inbounds m_row[7] = shape
    @inbounds m_row[8] = -(shape / scale)
    @inbounds m_row[9] = param.ρ
    # Column 10 is the divergence/polymorphism Ne ratio. `rates()` has none, so it is 1;
    # kept so tables stay interchangeable with the development branch's.
    @inbounds m_row[10] = 1.0

    @inbounds f_row[1] = f_n
    @inbounds f_row[2] = f_neg + f_pos_l + f_pos_h
    @inbounds f_row[3] = f_pos_l
    @inbounds f_row[4] = f_pos_h

    return nothing
end

"""
    simulate_models!(param, buf, ops, B_bins, var_params,
                     m_view, ps_view, pn_view, f_view)

One prior draw: solve `ppos` at B = 0.999, then walk the B-bin loop writing rows
directly into the provided sub-matrices. `ops` is `_tail_ops(...)`, aligned to `B_bins`.
On `set_ppos!` failure (the model is refused) all rows remain zero; the caller's matrix is
zero-initialised and `rates()` drops them.
"""
function simulate_models!(
    param::parameters,
    buf::SolveRatesBuffers,
    ops::Vector,
    B_bins::Vector{Float64},
    var_params::Vector{Float64},
    m_view::AbstractMatrix,
    ps_view::AbstractMatrix,
    pn_view::AbstractMatrix,
    f_view::AbstractMatrix,
)
    (gH, gL, gam_flanking, gam_dfe, shape, alpha, alpha_low, rho) = var_params

    param.shape = shape
    param.gam_dfe = gam_dfe
    param.scale = abs(shape / gam_dfe)
    param.gH = Int(gH)
    param.gL = Int(gL)
    param.gam_flanking = Int(gam_flanking)
    param.al_tot = alpha
    param.al_low = alpha * alpha_low
    # Must precede set_θ!, which consumes ρ.
    param.ρ = rho

    param.B = 0.999
    set_θ!(param)
    try
        set_ppos!(param; warm_start=buf.prev_ppos)
    catch
        return nothing
    end
    # Cache the successful (ppos_l, ppos_h) for the next model in this task chunk.
    # Only warm-start when we got a non-degenerate solution; otherwise the next
    # call falls back to the [0, 0] / [0.001, 0.001] fallback chain inside set_ppos!.
    if param.ppos_l > 0.0 && param.ppos_h > 0.0
        buf.prev_ppos[1] = param.ppos_l
        buf.prev_ppos[2] = param.ppos_h
    end

    @inbounds for j in 1:length(B_bins)
        param.B = B_bins[j]
        set_θ!(param)
        solve_rates!(
            view(m_view, j, :),
            view(ps_view, j, :),
            view(pn_view, j, :),
            view(f_view, j, :),
            buf, param, ops[j],
        )
    end
    return nothing
end

"""
	rates(param; iterations, output_file, ...)

Solve `iterations` random models across `B_bins` and write the ABC reference table:
models, per-DAC neutral and selected cumulative SFS, and fixation rates. One forward
model (dev docs/abcmk_validation.md 4v, the measured best): Phi published (eq 9); beneficial SFS at B level with sojourn at gamma*B; Phi flat on the SFS; guarded ppos solve; TailWeight projection.

# Arguments
 - `param::parameters`: mutable structure containing the model (`N`, `n`, `dac`, `shape`).
 - `iterations::Int64`: number of prior draws; each is solved in every B bin.
 - `output_file::String`: rates file to write (JLD2; appends a block under `N/n/cutoff`).
 - `gH`, `gL`: ranges of strong / weak beneficial `2N·s` (defaults `[100, 1000]`,
   `[1, 10]`; drawn as `floor(Uniform)`, so `gL` reaches 9, not 10).
 - `gam_flanking`, `gam_dfe`: ranges of flanking / coding deleterious `2N·s`
   (paper Table 2 defaults `[-1000, -500]`, `[-1000, -200]`).
 - `alpha::Vector{Float64}`: range of `al_tot` (paper default `[0.1, 0.9]`).
 - `rho`: `nothing` (default) holds `param.ρ`; a scalar holds that value; `[lo, hi]` draws
   `ρ ~ Uniform(lo, hi)` per model. Recorded in the `rho` column.
 - `B_bins`: background-selection bins (default `0.1:0.025:0.975` plus `0.999`).

The DFE shape is drawn `Uniform(shape·2^-2, shape·2^2)` (uniform on the published range,
kept by decision -- not the paper's log-uniform). A model whose `ppos` cannot be solved is
refused and dropped; the count is logged.

# Returns
 - `Dict` with `:models`, `:neut`, `:sel`, `:dsdn`, `:dac`; the same is written to
   `output_file`.
"""
function rates(
    param::parameters;
    iterations::Int64,
    output_file::String,
    gH::Vector{Int64}=[100,1000],
    gL::Vector{Int64}=[1,10],
    gam_flanking::Vector{Int64}=[-1000,-500],
    gam_dfe::Vector{Int64}=[-1000,-200],
    alpha::Vector{Float64}=[0.1, 0.9],
    rho::Union{Nothing,Real,AbstractVector{<:Real}}=nothing,
    B_bins::Vector{Float64}=vcat(collect(0.1:0.025:0.975), 0.999),
    weak::Bool=true,
    basesize::Int64=10,
    distributed::Bool=false,
)

    @info "Solving binomial convolution to downsample the SFS"
    binom = binom_op(param, B_bins)

    assertion_params(param)

    @unpack shape = param

    @info "Creating priors distributions"

    # Parameters ranges
    u_gh = Uniform(gH[1], gH[end])
    u_gl = Uniform(gL[1], gL[end])
    u_gam_flanking = Uniform(gam_flanking[1], gam_flanking[end])
    u_gam_dfe = Uniform(gam_dfe[1], gam_dfe[end])

    # afac = Uniform(0.5 * 2^-4, 0.5 * 2^4)
    afac = Uniform(shape * 2^-2, shape * 2^2)

    # afac  = Uniform(0.5,4)

    u_tot = Uniform(alpha[1], alpha[end])
    u_low = ifelse(weak, 0.0:0.05:1, [0])

    # Recombination rate: either drawn per model or held at a fixed scalar.
    # The `isnothing(u_rho)` branch below consumes NO RNG, so the default
    # (`rho = nothing`) leaves the random stream identical to the pre-`rho`
    # implementation. Unlike the other priors, the draw is not rounded: ρ is
    # recorded rather than inferred, so the column holds the exact value used.
    u_rho, rho_fixed = _rho_prior(rho, param)

    priors = Vector{Float64}[]
    for i::Int64 = 1:iterations
        push!(
            priors,
            [
                floor(rand(u_gh)),
                floor(rand(u_gl)),
                floor(rand(u_gam_flanking)),
                floor(rand(u_gam_dfe)),
                rand(afac),
                round(rand(u_tot), digits=2),
                round(rand(u_low), digits=2),
                isnothing(u_rho) ? rho_fixed : rand(u_rho),
            ],
        )
    end

    # Preallocate result matrices: one row per (iteration × B_bin)
    @info "Solving models in multiple threads"
    n_b = length(B_bins)
    total_rows = iterations * n_b
    n_dac = length(param.dac)
    m  = zeros(Float64, total_rows, 10)
    r_ps = zeros(Float64, total_rows, n_dac)
    r_pn = zeros(Float64, total_rows, n_dac)
    r_f  = zeros(Float64, total_rows, 4)

    # Position-indexed operator vector aligned to B_bins -- avoids per-call Float64
    # hash on `binom[param.B]` inside the inner loop.
    ops = _tail_ops(param, binom, B_bins)

    # Chunked parallelism: one param struct + one SolveRatesBuffers per task,
    # reused across all iterations assigned to that task.
    nt = max(1, Threads.nthreads())
    chunks = collect(Iterators.partition(1:iterations, max(1, cld(iterations, nt))))
    @sync for chunk in chunks
        Threads.@spawn begin
            p_local = copy_param(param)
            buf = SolveRatesBuffers(binom, n_dac)
            for i in chunk
                offset = (i - 1) * n_b
                rows = offset+1:offset+n_b
                simulate_models!(
                    p_local, buf, ops, B_bins, priors[i],
                    view(m,    rows, :),
                    view(r_ps, rows, :),
                    view(r_pn, rows, :),
                    view(r_f,  rows, :),
                )
            end
        end
    end

    @info "Saving solved models in $output_file"
    # Filter rows where simulate_models! returned all zeros (set_ppos! failed
    # or solve_rates! threw for every B_bin).
    idx = vec(sum(view(m, :, 1:8), dims=2) .!= 0)
    # F6/F9: a rejected row is a model `set_ppos!` could not solve -- refused rather than
    # written mis-specified.
    let n_drop = count(!, idx)
        n_drop > 0 && @info "$(n_drop) of $(length(idx)) model x B rows rejected " *
            "($(round(100n_drop/length(idx), digits=2))%)"
    end

    models = @view DataFrame(
        m,
        [:B, :al_low, :al_tot, :gam_flanking, :gL, :gH, :shape, :gam_dfe, :rho, :ne_ratio],
    )[
        idx,
        :,
    ]
    neut = @view r_ps[idx, :]
    sel  = @view r_pn[idx, :]
    dsdn = @view r_f[idx,  :]

    # Saving multiple summary statistics
    n = OrderedDict{Int,Array}()
    s = OrderedDict{Int,Array}()
    for i in eachindex(param.dac)
        n[param.dac[i]] = neut[:, i]
        s[param.dac[i]] = sel[:, i]
    end

    grp = _jld2_group_key(param.N, param.n, param.cutoff[1], param.cutoff[end])
    JLD2.jldopen(output_file, "a+") do file
        file[grp * "/models"]  = DataFrame(models; copycols=true)
        file[grp * "/neut"]    = n
        file[grp * "/sel"]     = s
        file[grp * "/dsdn"]    = Matrix{Float64}(dsdn)
        file[grp * "/dac"]     = Vector{Int}(param.dac)
        file[grp * "/options"] = _rates_options()
    end

    return (Dict(
        :models => models,
        :neut => n,
        :sel => s,
        :dsdn => dsdn,
        :dac => param.dac,
    ))
end

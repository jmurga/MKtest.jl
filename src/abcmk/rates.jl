################################
######   SFS projection   ######
################################
# The cumulative SFS at each DAC is projected with `TailWeights` (parameters.jl) instead
# of the full binomial projection: the same values, about 3x faster.

"""
    _tail_ops(param, binom, B_bins) -> Vector{TailWeights}

One `TailWeights` per B bin, in `B_bins` order.
"""
function _tail_ops(param::parameters,
                   binom::Dict{Float64,SparseMatrixCSC{Float64,Int64}},
                   B_bins::Vector{Float64})
    tw = tail_weights_op(param, binom, B_bins)
    return TailWeights[tw[b] for b in B_bins]
end

################################
###### Recombination prior ######
################################

"""
    _rho_prior(rho, param) -> (u_rho, rho_fixed)

Turn the `ρ` keyword of `rates()` into a `Uniform` prior (one draw per model) or a
fixed value.

| `rho`                  | result                                        |
|------------------------|-----------------------------------------------|
| `nothing` (default)    | `(nothing, param.ρ)`                          |
| `0.002::Real`          | `(nothing, 0.002)`                            |
| `[1e-4, 5e-3]::Vector` | `(Uniform(1e-4, 5e-3), NaN)`                  |
| `[0.001]` (length 1)   | `(nothing, 0.001)`                            |

The bounds are `rho[1]` and `rho[end]`, as for the other ranges. A fixed value draws no
random numbers, so a seeded run with a fixed `ρ` gives the same models as before.
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

################################
######     Rates file     ######
################################

# Reported by `_read_rates` for files that do not record their model.
const _NO_PROVENANCE = "unrecorded (file predates option provenance)"

# DFE shape prior of `rates`: `Uniform(shape * 2^-f, shape * 2^f)`, f = `shape_factor`,
# the half-width of the paper's exponent afac (a = a0 * 2^afac, afac on [-2, 2]). f = 2
# gives exactly `Uniform(shape/4, 4 shape)` (powers of two are exact), so tables built
# before the keyword keep their bits.
function _shape_prior(shape::Real, shape_factor::Real)
    (isfinite(shape_factor) && shape_factor > 0) ||
        throw(ArgumentError("shape_factor must be finite and > 0, got $shape_factor"))
    f = Float64(shape_factor)
    return Uniform(shape * 2.0^-f, shape * 2.0^f)
end

"""
    _rates_options() -> String

Description of the forward model, stored in the rates file next to the table. Tables
without it were built by earlier versions of MKtest, with a different model.
"""
_rates_options() = "model=best-2026-09-28 (Phi published (eq 9); beneficial SFS at B level with sojourn at gamma*B; Phi flat on the SFS; guarded ppos solve); projection=TailWeight"

_jld2_group_key(N::Int, n::Int, cutoff_lo, cutoff_hi) =
    string(N) * "/" * string(n) * "/cutoff=[" * string(cutoff_lo) * "," * string(cutoff_hi) * "]"

"""
    _read_rates(input_file, param) -> NamedTuple

`(; models, neut, sel, dsdn, dac, options)` of the `(N, n, cutoff)` block of a rates file.
Reads files written by this and earlier versions.
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

"""
    RatesTable

One `(N, n, cutoff)` block of a rates file, loaded by [`load_rates`](@ref) and passed to
[`abcmk`](@ref) or [`summary_statistics`](@ref) as `rates = t` instead of a file path.
Fields: `file`, `N`, `n`, `cutoff`, `models`, `neut`, `sel`, `dsdn`, `dac`, `options` (as
[`rates`](@ref) wrote them), plus a cache of row-prior weights, computed by the first call
that needs them and reused by every later call with the same `row_prior`, `alpha` and
`B_bins`.
"""
struct RatesTable{TM,TN,TS,TD}
    file::String
    N::Int
    n::Int
    cutoff::Tuple{Float64,Float64}
    models::TM
    neut::TN
    sel::TS
    dsdn::TD
    dac::Vector{Int}
    options::String
    weights::Dict{Any,Union{Nothing,Vector{Float64}}}
    lock::ReentrantLock
end

function RatesTable(file::String, param::parameters, tmp)
    return RatesTable(file, param.N, param.n, (Float64(param.cutoff[1]), Float64(param.cutoff[end])),
                      tmp.models, tmp.neut, tmp.sel, tmp.dsdn, Vector{Int}(tmp.dac), String(tmp.options),
                      Dict{Any,Union{Nothing,Vector{Float64}}}(), ReentrantLock())
end

"""
    load_rates(h5_file, param) -> RatesTable

Read the `(param.N, param.n, param.cutoff)` block of a rates file once, for several
[`abcmk`](@ref) or [`summary_statistics`](@ref) calls in the same Julia session:

```julia
t = MKtest.load_rates("rates.jld2", param)
cases    = MKtest.abcmk(param, sfs_cases, div_cases; rates = t, summstat_size = 1e5, in_memory = true)
controls = MKtest.abcmk(param, sfs_ctl, div_ctl; rates = t, summstat_size = 1e5, in_memory = true)
```

Each call then skips the read and, from the second call with the same `row_prior`,
`alpha` and `B_bins`, the row-prior weights. The results are the same as with
`h5_file = "rates.jld2"`: the table is not changed and the random draws are the same.
The table stays in memory as long as `t` does.
"""
load_rates(h5_file::String, param::parameters) = RatesTable(h5_file, param, _read_rates(h5_file, param))

Base.show(io::IO, t::RatesTable) =
    print(io, "RatesTable(\"", t.file, "\", N = ", t.N, ", n = ", t.n, ", cutoff = ",
          collect(t.cutoff), ", ", size(t.dsdn, 1), " rows, dac = ", t.dac, ")")

_check_rates_source(h5_file, rates) =
    isnothing(h5_file) == isnothing(rates) && throw(ArgumentError(
        "give the rates table as h5_file = \"file.jld2\" or rates = (a path or a " *
        "load_rates table), exactly one of them"))

# The table to use: a loaded `rates` table (checked against `param`) or a file read now,
# from `rates` or `h5_file` (exactly one of the two).
function _rates_table(h5_file, rates, param::parameters)
    _check_rates_source(h5_file, rates)
    src = isnothing(rates) ? h5_file : rates
    src isa Union{AbstractString,RatesTable} || throw(ArgumentError(
        "rates must be a path or a load_rates table, got $(typeof(src))"))
    if src isa RatesTable
        (src.N, src.n, src.cutoff) == (param.N, param.n, (Float64(param.cutoff[1]), Float64(param.cutoff[end]))) ||
            throw(ArgumentError("the loaded rates table is N = $(src.N), n = $(src.n), cutoff = " *
                                "$(collect(src.cutoff)); param is N = $(param.N), n = $(param.n), " *
                                "cutoff = $(param.cutoff)"))
        missing_dac = setdiff(param.dac, src.dac)
        isempty(missing_dac) || throw(ArgumentError(
            "the loaded rates table has no DAC $(missing_dac) (it has $(src.dac))"))
        @info "Random sampling expected SFS and fixation rates from the loaded table $(src.file)"
        return src
    end
    @info "Opening and random sampling expected SFS and fixation rates"
    return load_rates(String(src), param)
end

# Row-prior weights over `idx`, computed once per (row_prior, alpha, B_bins) and table.
function _table_weights(t::RatesTable, row_prior, alpha, B_bins, idx)
    key = (row_prior, isnothing(alpha) ? nothing : copy(alpha), isnothing(B_bins) ? nothing : copy(B_bins))
    lock(t.lock) do
        get!(() -> _row_weights(row_prior, t.dsdn, idx; B = t.models.B), t.weights, key)
    end
end

################################
######   Solving models   ######
################################
# Each task reuses one `SolveRatesBuffers` for all its models and B bins, and writes the
# results straight into the output matrices, so `rates()` barely allocates and runs faster.

"""
    SolveRatesBuffers(binom, n_dac)

Work buffers for one task of `rates()`. `scratch` holds the source vectors of
`_src_pos!`/`_src_neg!`, `acc_tmp` each projection and `acc_sel` their sum.
"""
struct SolveRatesBuffers
    scratch::Vector{Float64}
    acc_sel::Vector{Float64}
    acc_tmp::Vector{Float64}
    # ppos of the previous model of this task: the first starting point of set_ppos!,
    # which then needs fewer iterations.
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
    @unpack p_strong, p_weak, γ_weak, γ_strong, dac = param
    has_h = p_strong > 0.0
    has_l = p_weak > 0.0
    W = op.W
    ncols = size(W, 2)
    nd = length(dac)

    # Neutral: precomputed, it does not depend on the model.
    nsc = param.B * param.θ_coding * param.f_syn
    @inbounds @simd for k in 1:nd
        ps_row[k] = nsc * op.neut_tail[k]
    end

    # Selected: (sel_h + sel_l) + sel_neg, the beneficial classes as in `sfs_pos`.
    acc = buf.acc_sel
    fill!(acc, 0.0)
    if has_h
        src = _src_pos!(buf.scratch, param, γ_strong * param.B, p_strong, ncols)
        mul!(buf.acc_tmp, W, src)
        fh = param.B * param.θ_coding * Φ(param, γ_strong) * (1 - param.f_syn)
        @inbounds @simd for k in 1:nd
            acc[k] += fh * buf.acc_tmp[k]
        end
    end
    if has_l
        src = _src_pos!(buf.scratch, param, γ_weak * param.B, p_weak, ncols)
        mul!(buf.acc_tmp, W, src)
        fl = param.B * param.θ_coding * Φ(param, γ_weak) * (1 - param.f_syn)
        @inbounds @simd for k in 1:nd
            acc[k] += fl * buf.acc_tmp[k]
        end
    end
    let src = _src_neg!(buf.scratch, param, p_weak + p_strong, ncols)
        mul!(buf.acc_tmp, W, src)
        fneg = param.B * param.θ_coding * (1 - param.f_syn)
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

Write the model row, the neutral and selected SFS at each DAC, and the fixation rates of
one B bin into the given rows, without allocating.
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
    @unpack B, p_strong, p_weak, γ_weak, γ_strong, α_weak, α,
            γ_flanking, shape, scale, dac = param

    has_h = p_strong > 0.0
    has_l = p_weak > 0.0

    # Fixation
    f_n = B * fix_neut(param)
    f_neg = B * fix_neg(param, 0.5 * p_strong + 0.5 * p_weak)
    f_pos_l = has_l ? fix_pos_sim(param, γ_weak, 0.5 * p_weak) : 0.0
    f_pos_h = has_h ? fix_pos_sim(param, γ_strong, 0.5 * p_strong) : 0.0

    _sfs_dac!(ps_row, pn_row, buf, param, op)

    @inbounds m_row[1] = B * 1000
    @inbounds m_row[2] = α_weak
    @inbounds m_row[3] = α
    @inbounds m_row[4] = γ_flanking
    @inbounds m_row[5] = γ_weak
    @inbounds m_row[6] = γ_strong
    @inbounds m_row[7] = shape
    @inbounds m_row[8] = -(shape / scale)
    @inbounds m_row[9] = param.ρ
    # Column 10 (ne_ratio) is always 1 in `rates()`.
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

Solve one prior draw: `ppos` at B = 0.999, then every B bin, written into the given
rows. If `set_ppos!` fails the rows stay zero and `rates()` drops the model.
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
    (γ_strong, γ_weak, γ_flanking, γ_deleterious, shape, alpha, alpha_low, rho) = var_params

    param.shape = shape
    param.γ_deleterious = γ_deleterious
    param.scale = abs(shape / γ_deleterious)
    param.γ_strong = Int(γ_strong)
    param.γ_weak = Int(γ_weak)
    param.γ_flanking = Int(γ_flanking)
    param.α = alpha
    param.α_weak = alpha * alpha_low
    # Before set_θ!, which uses ρ.
    param.ρ = rho

    param.B = 0.999
    set_θ!(param)
    route = try
        _set_ppos!(param; warm_start=buf.prev_ppos)
    catch
        return nothing
    end
    # Starting point for the next model of this task (only if both are non-zero). A
    # `:linear` root (nlsolve failed, p ~ 1e-7) is not reused: the next draw then starts
    # where it did before the closed form existed, so its bits do not change.
    if route === :nlsolve && param.p_weak > 0.0 && param.p_strong > 0.0
        buf.prev_ppos[1] = param.p_weak
        buf.prev_ppos[2] = param.p_strong
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
    rates(param; iterations, output_file, γ_strong = [100, 1000], γ_weak = [1, 10],
          γ_flanking = [-1000, -500], γ_deleterious = [-1000, -200], α = [0.1, 0.9],
          ρ = nothing, shape_factor = 2, B_bins = vcat(0.1:0.025:0.975, 0.999),
          weak = true)

Build the ABC-MK reference table: draw `iterations` models from the priors, solve each in
every `B` bin, and append the table to the JLD2 file `output_file`.

For each draw, [`set_ppos!`](@ref) solves the beneficial probabilities at `B = 0.999`;
draws it cannot solve are dropped (the count is logged). Then, for each `B`, the model's
expected synonymous and nonsynonymous SFS tails (DAC >= k) at every `param.dac` and its
fixation rates are computed. The draws run on all Julia threads (`julia -t auto`); with a
seed, the table reproduces only at the same thread count.

# Arguments
 - `param::parameters`: the model; `N`, `n`, `dac`, `cutoff` and `shape` are used as
   set, the other drawn fields are overwritten in per-task copies.

# Keywords
 - `iterations::Int64`: number of prior draws; the table has up to
   `iterations * length(B_bins)` rows.
 - `output_file::String`: JLD2 file, opened in append mode; the block is written under
   the group `"N/n/cutoff=[lo,hi]"`, which must not exist yet.
 - `γ_strong`, `γ_weak::Vector{Int64}`: ranges of the strongly and weakly beneficial coefficients,
   drawn as `floor(Uniform(lo, hi))` (so `γ_weak` reaches 9, not 10).
 - `γ_flanking`, `γ_deleterious::Vector{Int64}`: ranges of the flanking and coding
   deleterious coefficients (mean of the DFE), drawn the same way.
 - `α::Vector{Float64}`: range of `param.α`, drawn uniformly and rounded to 2 digits.
 - `ρ`: `nothing` keeps `param.ρ`; a positive number fixes ρ; `[lo, hi]` draws
   `ρ ~ Uniform(lo, hi)` per model. A fixed ρ draws no random numbers.
 - `shape_factor::Real`: width of the DFE shape prior, as the exponent half-width
   (default 2, must be > 0); see below.
 - `B_bins::Vector{Float64}`: background selection values (37 by default).
 - `weak::Bool`: draw the weak share of α from `0:0.05:1` (`α_weak = α * share`);
   `false` sets `α_weak = 0`.
 - `basesize`, `distributed`: not used.

The DFE shape is drawn from `Uniform(shape 2^-f, shape 2^f)` with `f = shape_factor`
(default 2: `[shape/4, 4 shape]`, the paper's range for `afac` in `a = a0 2^afac`; the
paper draws `afac` uniformly, which makes its shape prior log-uniform). Above shape ~1 at the default coefficients almost no
deleterious mutation fixes, so such a draw keeps a realistic Dn/Ds only at an α close to 1.
All coefficients are `2N` times the heterozygote selection coefficient.

# Returns
A `Dict` with the same content as the file block:
 - `:models`: `DataFrame`, one row per (draw, B). Its columns keep their stored names, so
   earlier tables stay readable: `B` (B x 1000), `al_low` (`α_weak`), `al_tot` (`α`),
   `gam_flanking` (`γ_flanking`), `gL` (`γ_weak`), `gH` (`γ_strong`), `shape`, `gam_dfe`
   (`γ_deleterious`), `rho` (`ρ`), `ne_ratio` (always 1).
 - `:neut`, `:sel`: `OrderedDict{Int,Array}` from each DAC to the expected synonymous and
   nonsynonymous tails of every row.
 - `:dsdn`: matrix of fixation rates, columns `Ds`, `Dn`, and the parts of `Dn` from
   weakly and from strongly beneficial mutations.
 - `:dac`: `param.dac`.

The file block also stores `options`, a description of the forward model.
"""
function rates(
    param::parameters;
    iterations::Int64,
    output_file::String,
    γ_strong::Vector{Int64}=[100,1000],
    γ_weak::Vector{Int64}=[1,10],
    γ_flanking::Vector{Int64}=[-1000,-500],
    γ_deleterious::Vector{Int64}=[-1000,-200],
    α::Vector{Float64}=[0.1, 0.9],
    ρ::Union{Nothing,Real,AbstractVector{<:Real}}=nothing,
    shape_factor::Real=2,
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
    u_gh = Uniform(γ_strong[1], γ_strong[end])
    u_gl = Uniform(γ_weak[1], γ_weak[end])
    u_gam_flanking = Uniform(γ_flanking[1], γ_flanking[end])
    u_gam_dfe = Uniform(γ_deleterious[1], γ_deleterious[end])

    # DFE shape prior
    afac = _shape_prior(shape, shape_factor)

    u_tot = Uniform(α[1], α[end])
    u_low = ifelse(weak, 0.0:0.05:1, [0])

    # Recombination rate: drawn per model or fixed. A fixed value draws no random
    # numbers. The draw is not rounded: the `rho` column stores the value used.
    u_rho, rho_fixed = _rho_prior(ρ, param)

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

    # One TailWeights per B bin, in B_bins order.
    ops = _tail_ops(param, binom, B_bins)

    # One parameter copy and one SolveRatesBuffers per task.
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
    # Drop the models set_ppos! could not solve (their rows are all zero).
    idx = vec(sum(view(m, :, 1:8), dims=2) .!= 0)
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

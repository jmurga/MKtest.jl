################################
###    Summary statistics    ###
################################

# The Poisson samplers write into preallocated `SamplingBuffers` in single loops instead
# of allocating intermediate arrays, which makes `summary_statistics` faster.
# Layout: neut, sel, sampled_pn, sampled_ps are (n_dac × N); α_mat, ω_mat are (N × 3) and
# (N × 4); α_summaries is (N × n_dac).

"""
    SamplingBuffers(N, n_dac)

Preallocated work buffers for `sampling_summaries!`, one per observed dataset.
"""
struct SamplingBuffers
    N::Int
    n_dac::Int
    sampled_dn::Vector{Float64}
    sampled_ds::Vector{Float64}
    α_mat::Matrix{Float64}        # N × 3 -- α_weak, α_strong, α_total
    ω_mat::Matrix{Float64}        # N × 4 -- ωₐ_weak, ωₐ_strong, ωₐ, ωₙₐ
    sampled_pn::Matrix{Float64}   # n_dac × N
    sampled_ps::Matrix{Float64}   # n_dac × N
    α_summaries::Matrix{Float64}  # N × n_dac
    mask::Vector{Bool}            # N -- single-pass filter mask
    expected_values::Matrix{Float64}  # N × (12 + n_dac) -- column-major assembly
end

function SamplingBuffers(N::Int, n_dac::Int)
    SamplingBuffers(
        N, n_dac,
        Vector{Float64}(undef, N),       # sampled_dn
        Vector{Float64}(undef, N),       # sampled_ds
        Matrix{Float64}(undef, N, 3),    # α_mat
        Matrix{Float64}(undef, N, 4),    # ω_mat
        Matrix{Float64}(undef, n_dac, N),# sampled_pn
        Matrix{Float64}(undef, n_dac, N),# sampled_ps
        Matrix{Float64}(undef, N, n_dac),# α_summaries
        Vector{Bool}(undef, N),          # mask
        Matrix{Float64}(undef, N, 12 + n_dac),  # expected_values (max width)
    )
end

################################
######  Poisson sampling  ######
################################

"""
    AbstractPoissonSampling

How `sampling_summaries!` draws the Poisson realisations that become the ABC summary
statistics: [`RjSampling`](@ref) (the default), which accounts for demography in the
polymorphism, or [`IndependentSampling`](@ref) (ABC-MK v1.0.0).
"""
abstract type AbstractPoissonSampling end

"""
    IndependentSampling()

Poisson sampler of ABC-MK v1.0.0: every class is drawn as an independent Poisson whose
mean is the model's share of the observed total.

Let `O_k` be the observed polymorphism at DAC >= k (Pn + Ps) and `D = Dn + Ds` the
observed divergence. For a table row with tail rates `λps_k`, `λpn_k` and fixation rates
`λds`, `λdn`, `λW`, `λS`:

    ps_k ~ Poisson(λps_k / (λps_k + λpn_k) O_k)    pn_k ~ Poisson(λpn_k / (λps_k + λpn_k) O_k)
    ds   ~ Poisson(λds / (λds + λdn) D)            dn   ~ Poisson(λdn / (λds + λdn) D)
    sw   ~ Poisson(λW / (λds + λdn) D)             ss   ~ Poisson(λS / (λds + λdn) D)

`sw + ss > dn` can happen, so a sampled α above 1 is kept as a sampling artifact; only
non-finite rows are dropped. Pass it as `sampling` to [`summary_statistics`](@ref) or
[`abcmk`](@ref).
"""
struct IndependentSampling <: AbstractPoissonSampling end

"""
    poisson_fixation!(α_mat, ω_mat, sampled_dn, sampled_ds, observed_values, λds, λdn, λweak, λstrong)

In-place Poisson sampling for fixation rates. Writes results directly into the
caller-provided buffers. When `λds + λdn == 0` the row is marked NaN so the downstream
filter drops it.
"""
function poisson_fixation!(
    α_mat::AbstractMatrix,
    ω_mat::AbstractMatrix,
    sampled_dn::AbstractVector,
    sampled_ds::AbstractVector,
    observed_values::AbstractVector,
    λds::AbstractVector,
    λdn::AbstractVector,
    λweak::AbstractVector,
    λstrong::AbstractVector,
)
    return poisson_fixation!(α_mat, ω_mat, sampled_dn, sampled_ds, observed_values,
                             λds, λdn, λweak, λstrong, IndependentSampling())
end

function poisson_fixation!(
    α_mat::AbstractMatrix,
    ω_mat::AbstractMatrix,
    sampled_dn::AbstractVector,
    sampled_ds::AbstractVector,
    observed_values::AbstractVector,
    λds::AbstractVector,
    λdn::AbstractVector,
    λweak::AbstractVector,
    λstrong::AbstractVector,
    ::IndependentSampling,
)
    observed_d, ln, ls = observed_values
    @inbounds for i in eachindex(λds)
        denom = λds[i] + λdn[i]
        if denom == 0.0
            sampled_ds[i] = NaN
            sampled_dn[i] = NaN
            α_mat[i, 1] = NaN; α_mat[i, 2] = NaN; α_mat[i, 3] = NaN
            ω_mat[i, 1] = NaN; ω_mat[i, 2] = NaN; ω_mat[i, 3] = NaN; ω_mat[i, 4] = NaN
            continue
        end
        inv_denom = 1.0 / denom
        sampled_ds[i] = pois_rand(λds[i] * inv_denom * observed_d)
        sampled_dn[i] = pois_rand(λdn[i] * inv_denom * observed_d)
        sw = pois_rand(λweak[i] * inv_denom * observed_d)
        ss = pois_rand(λstrong[i] * inv_denom * observed_d)

        dn_i = sampled_dn[i]
        ds_i = sampled_ds[i]

        # α columns
        α_mat[i, 1] = sw / dn_i              # α_weak
        α_mat[i, 3] = (sw + ss) / dn_i       # α_total
        α_mat[i, 2] = α_mat[i, 3] - α_mat[i, 1]   # α_strong

        # ω columns -- guard against ds_i == 0 => NaN
        dₛ_val = ds_i / ls
        ω_mat[i, 1] = (sw / ln) / dₛ_val
        ω_mat[i, 2] = (ss / ln) / dₛ_val
        ω_mat[i, 3] = ((sw + ss) / ln) / dₛ_val
        ω_mat[i, 4] = ((dn_i - ss - sw) / ln) / dₛ_val
    end
    return nothing
end

"""
    poisson_polymorphism!(sampled_pn, sampled_ps, observed_values, λps, λpn)

In-place Poisson sampling for polymorphism. `λps, λpn` are `(n_dac × N)`; outputs
match. A degenerate row (`λps + λpn == 0`) gets 0, so the downstream ratio becomes NaN
and the row is filtered.
"""
function poisson_polymorphism!(
    sampled_pn::AbstractMatrix,
    sampled_ps::AbstractMatrix,
    observed_values::AbstractVector,
    λps::AbstractMatrix,
    λpn::AbstractMatrix,
)
    return poisson_polymorphism!(sampled_pn, sampled_ps, observed_values, λps, λpn,
                                 IndependentSampling())
end

function poisson_polymorphism!(
    sampled_pn::AbstractMatrix,
    sampled_ps::AbstractMatrix,
    observed_values::AbstractVector,
    λps::AbstractMatrix,
    λpn::AbstractMatrix,
    ::IndependentSampling,
)
    n_dac, N = size(λps)
    @inbounds for j in 1:N
        for i in 1:n_dac
            denom = λps[i, j] + λpn[i, j]
            if denom == 0.0
                sampled_pn[i, j] = 0.0
                sampled_ps[i, j] = 0.0
            else
                ov = observed_values[i]
                inv_denom = 1.0 / denom
                sampled_ps[i, j] = pois_rand(λps[i, j] * inv_denom * ov)
                sampled_pn[i, j] = pois_rand(λpn[i, j] * inv_denom * ov)
            end
        end
    end
    return nothing
end

"""
    _prepare_sampling(sampling, sfs, dac, nn) -> AbstractPoissonSampling

Per-dataset form of the sampler, built once per observed dataset in `summary_statistics`
from the RAW per-frequency SFS `[DAC pn ps]`; `nn` is the number of sampled chromosomes
(`param.nn = 2n`). The default returns `sampling` itself. A sampler that needs the
observed spectrum ([`RjSampling`](@ref)) returns a filled copy.
"""
_prepare_sampling(sampling::AbstractPoissonSampling, sfs::AbstractMatrix, dac::Vector{Int64}, nn::Integer) = sampling

"""
    RjSampling(; w = (0.0, 2.5), anchor = (0.15, 0.70), tilt = 0.05)

[`IndependentSampling`](@ref) with a demographic nuisance on the polymorphism: the
data's synonymous SFS distortion `r_i` (GRAPES' r_j), normalised on an anchor band, is
applied to the model's synonymous class and, raised to a random power `w`, to its
nonsynonymous class. With `tilt > 0` the nonsynonymous class also gets a random frequency
tilt. Fixations are drawn as in `IndependentSampling`.

The default sampler of [`summary_statistics`](@ref) and [`abcmk`](@ref), used there with
[`FlatAlphaPrior`](@ref) and `stat = "median"`.

# Keywords
 - `w = (0.0, 2.5)`: range of `w ~ Uniform(w...)`, drawn per table row. `w = 0`:
   nonsynonymous sites ignore the distortion; `w = 1`: the same distortion on both
   classes (GRAPES); `w > 1`: stronger on nonsynonymous sites. Needs `0 <= w[1] <= w[2]`.
 - `anchor = (0.15, 0.70)`: band of derived-allele frequency on which `r` is normalised,
   so that `r` changes the shape of the spectra but never the Pn/Ps level, which only α
   sets. Needs `0 < anchor[1] < anchor[2] <= 1`.
 - `tilt = 0.05`: range of `τ ~ Uniform(-tilt, tilt)`; `0` draws no `τ`.

# Steps, per dataset and table row
1. The DACs of `param.dac` split the spectrum into intervals `[dac_i, dac_(i+1) - 1]`, the
   last one up to `nn - 1`. The row's interval masses `e_s,i`, `e_n,i` are differences of
   its tails.
2. `r_i = (Ps_obs,i / sum Ps_obs) / (e_s,i / sum e_s)`, divided by its `e_s`-weighted mean
   over the intervals lying wholly inside the anchor band (`anchor[1] <= dac_i/nn` and
   `stop_i/nn <= anchor[2]`).
3. The synonymous masses become `e_s,i r_i` and the nonsynonymous masses `e_n,i r_i^w`.
4. If `tilt > 0`, the nonsynonymous masses are also multiplied by `f_i^τ / c_τ`, with
   `f_i = sqrt(dac_i stop_i) / nn` and `c_τ` chosen so the Pn mass on the anchor band is
   unchanged: a zero-mean change of shape aimed at expansions.
5. `ps_k` and `pn_k` are drawn from the distorted tails as in `IndependentSampling`.

`w` and `τ` are marginalised: they are not written to the summaries.

# Requirements
`param.dac` ascending, Ps > 0 in every DAC interval of each dataset, and at least one
whole interval inside the anchor band (at n = 50, DACs 20 and 50 give the interval
20-49). Each violation throws an `ArgumentError`; [`valid_data`](@ref) finds the datasets
that would fail. Pass the RAW per-frequency SFS: `summary_statistics` builds the tails.

The defaults were selected on SLiM simulations (two-epoch histories at n = 50 and the
Tennessen history at n = 661), where they gave an α RMSE of 0.029 and 0.030 against
0.055 and 0.039 for the v1.0.0 sampler and prior.
"""
struct RjSampling <: AbstractPoissonSampling
    w_lo::Float64
    w_hi::Float64
    anchor_lo::Float64        # anchor band, as a frequency (fraction of the 2n chromosomes)
    anchor_hi::Float64
    ps_obs::Vector{Float64}   # observed Ps per interval (empty until _prepare_sampling)
    band::Vector{Bool}        # intervals inside the anchor band
    tilt::Float64             # τ ~ U(-tilt, tilt); 0 = no tilt
    lnf::Vector{Float64}      # ln geometric-mid frequency per interval (tilt coordinate)
end

RjSampling(w_lo, w_hi, anchor_lo, anchor_hi, ps_obs, band) =
    RjSampling(w_lo, w_hi, anchor_lo, anchor_hi, ps_obs, band, 0.0, Float64[])

function RjSampling(; w = (0.0, 2.5), anchor = (0.15, 0.70), tilt = 0.05)
    lo, hi = Float64(first(w)), Float64(last(w))
    (isfinite(lo) && isfinite(hi) && 0 <= lo <= hi) ||
        throw(ArgumentError("w must satisfy 0 <= w[1] <= w[2] < Inf, got $w"))
    a_lo, a_hi = Float64(first(anchor)), Float64(last(anchor))
    0 < a_lo < a_hi <= 1 || throw(ArgumentError(
        "anchor must be a frequency band (lo, hi) with 0 < lo < hi <= 1, got $anchor"))
    (isfinite(tilt) && tilt >= 0) || throw(ArgumentError("tilt must be finite and >= 0, got $tilt"))
    return RjSampling(lo, hi, a_lo, a_hi, Float64[], Bool[], Float64(tilt), Float64[])
end

function _prepare_sampling(s::RjSampling, sfs::AbstractMatrix,
                           dac::Vector{Int64}, nn::Integer)
    issorted(dac) || throw(ArgumentError("RjSampling needs param.dac ascending"))
    ps_obs = _interval_ps(sfs, dac)
    n = length(dac)
    stops = [i < n ? dac[i+1] - 1 : nn - 1 for i in 1:n]
    # the band is a frequency; resolve it with this dataset's sample size
    band = [dac[i] >= s.anchor_lo * nn && stops[i] <= s.anchor_hi * nn for i in 1:n]
    any(band) || throw(ArgumentError(
        "no interval of dac = $dac lies inside the anchor band $(s.anchor_lo)-$(s.anchor_hi) " *
        "(DAC $(ceil(Int, s.anchor_lo * nn))-$(floor(Int, s.anchor_hi * nn)) at 2n = $nn): " *
        "add DACs so that one interval [dac[i], dac[i+1] - 1] falls inside it"))
    all(>(0), ps_obs) || throw(ArgumentError(
        "RjSampling needs Ps > 0 in every DAC interval; got $ps_obs for dac = $dac"))
    lnf = [log(sqrt(dac[i] * stops[i]) / nn) for i in 1:n]
    return RjSampling(s.w_lo, s.w_hi, s.anchor_lo, s.anchor_hi, ps_obs, band, s.tilt, lnf)
end

"Observed Ps per DAC interval `[dac_i, dac_(i+1) - 1]` (the last runs to the top) of a RAW SFS."
function _interval_ps(sfs::AbstractMatrix, dac::Vector{Int64})
    cum = cumulative_sfs(sfs)
    tail_ps = map(z -> sum(cum[cum[:, 1] .== z, 3]), dac)
    n = length(dac)
    return [i < n ? tail_ps[i] - tail_ps[i+1] : tail_ps[n] for i in 1:n]
end

"""
    valid_data(param, sfs, divergence; sampling = RjSampling()) -> Vector{Bool}

Which datasets `summary_statistics` and [`abcmk`](@ref) can use with `sampling`, one `Bool`
per dataset, so that datasets that would stop the run are discarded beforehand:

    ok = MKtest.valid_data(param, sfs, divergence)
    sfs, divergence = sfs[ok], divergence[ok]

Every sampler needs `Dn > 0`, `Ds > 0` and, at every DAC of `param.dac`, Pn > 0 and Ps > 0
in the tail (DAC >= k), so that the observed α(x) is finite. [`RjSampling`](@ref) also
needs Ps > 0 in every interval between consecutive DACs, where it measures the
synonymous distortion. The DAC layout itself (ascending, one interval inside the anchor
band) is not checked here: it is the same for every dataset, and `RjSampling` still
throws for it.

`sfs` and `divergence` as for `summary_statistics`: RAW per-frequency SFS `[DAC Pn Ps]`
and divergence `[Dn, Ds, ...]`. Throws an `ArgumentError` if their lengths differ.
"""
function valid_data(param::parameters, sfs::Vector, divergence::Vector;
                    sampling::AbstractPoissonSampling = RjSampling())
    length(sfs) == length(divergence) || throw(ArgumentError(
        "got $(length(sfs)) SFS for $(length(divergence)) divergence datasets"))
    return [_valid_data(sampling, x, d, param.dac) for (x, d) in zip(sfs, divergence)]
end

function _valid_data(::AbstractPoissonSampling, sfs::AbstractMatrix, d, dac::Vector{Int64})
    d[1] > 0 && d[2] > 0 || return false
    cum = cumulative_sfs(sfs)
    tails = cum[in(dac).(view(cum, :, 1)), 2:3]
    return size(tails, 1) == length(dac) && all(>(0), tails)
end

_valid_data(s::RjSampling, sfs::AbstractMatrix, d, dac::Vector{Int64}) =
    invoke(_valid_data, Tuple{AbstractPoissonSampling,AbstractMatrix,Any,Vector{Int64}}, s, sfs, d, dac) &&
    all(>(0), _interval_ps(sfs, dac))

"""
    _rj_tails!(tps, tpn, r, λps_col, λpn_col, s, w, τ = 0.0)

Distorted tail sums of one table row for [`RjSampling`](@ref): fills
`tps[k] = sum_{i>=k} λps_i r_i`, `tpn[k] = sum_{i>=k} λpn_i r_i^w t_i` from the row's TAIL
sums `λps_col`, `λpn_col` (split into intervals by differencing), with `r` normalised on
the anchor band and `t_i = exp(τ ln f_i) / c_τ` (1 if `τ = 0`) normalised so the band's
Pn mass is unchanged. `r` is scratch.
"""
function _rj_tails!(tps, tpn, r, λps_col, λpn_col, s::RjSampling, w::Float64, τ::Float64 = 0.0)
    n = length(λps_col)
    tot_obs = sum(s.ps_obs)
    tot_mod = 0.0
    @inbounds for i in 1:n
        tot_mod += max(0.0, i < n ? λps_col[i] - λps_col[i+1] : λps_col[n])
    end
    num = 0.0; den = 0.0
    @inbounds for i in 1:n
        e = max(0.0, i < n ? λps_col[i] - λps_col[i+1] : λps_col[n])
        r[i] = e > 0.0 && tot_mod > 0.0 ? (s.ps_obs[i] / tot_obs) / (e / tot_mod) : 1.0
        if s.band[i]
            num += e * r[i]; den += e
        end
    end
    c = den > 0.0 && num > 0.0 ? num / den : 1.0
    # tilt: normalise exp(τ ln f) on the band's (r-distorted) Pn mass
    cτ = 1.0
    if τ != 0.0
        tn = 0.0; td = 0.0
        @inbounds for i in 1:n
            s.band[i] || continue
            x = max(0.0, i < n ? λpn_col[i] - λpn_col[i+1] : λpn_col[n]) * (r[i] / c)^w
            tn += x * exp(τ * s.lnf[i]); td += x
        end
        cτ = td > 0.0 && tn > 0.0 ? tn / td : 1.0
    end
    acc_s = 0.0; acc_n = 0.0
    @inbounds for i in n:-1:1
        e_s = max(0.0, i < n ? λps_col[i] - λps_col[i+1] : λps_col[n])
        e_n = max(0.0, i < n ? λpn_col[i] - λpn_col[i+1] : λpn_col[n])
        ri = r[i] / c
        acc_s += e_s * ri
        acc_n += τ == 0.0 ? e_n * ri^w : e_n * ri^w * exp(τ * s.lnf[i]) / cτ
        tps[i] = acc_s
        tpn[i] = acc_n
    end
    return nothing
end

function poisson_fixation!(
    α_mat::AbstractMatrix,
    ω_mat::AbstractMatrix,
    sampled_dn::AbstractVector,
    sampled_ds::AbstractVector,
    observed_values::AbstractVector,
    λds::AbstractVector,
    λdn::AbstractVector,
    λweak::AbstractVector,
    λstrong::AbstractVector,
    s::RjSampling,
)
    return poisson_fixation!(α_mat, ω_mat, sampled_dn, sampled_ds, observed_values,
                             λds, λdn, λweak, λstrong, IndependentSampling())
end

function poisson_polymorphism!(
    sampled_pn::AbstractMatrix,
    sampled_ps::AbstractMatrix,
    observed_values::AbstractVector,
    λps::AbstractMatrix,
    λpn::AbstractMatrix,
    s::RjSampling,
)
    n_dac, N = size(λps)
    length(s.ps_obs) == n_dac || throw(ArgumentError(
        "RjSampling was not prepared for $n_dac DACs (call it through " *
        "summary_statistics, which fills the observed Ps per dataset)"))
    tps = zeros(n_dac); tpn = zeros(n_dac); r = zeros(n_dac)
    fixed_w = s.w_lo == s.w_hi
    u = fixed_w ? nothing : Distributions.Uniform(s.w_lo, s.w_hi)
    tilt = s.tilt > 0.0
    tilt && length(s.lnf) != n_dac && throw(ArgumentError(
        "RjSampling with tilt was not prepared for $n_dac DACs"))
    @inbounds for j in 1:N
        w = fixed_w ? s.w_lo : rand(u)
        τ = tilt ? s.tilt * (2 * rand() - 1) : 0.0
        _rj_tails!(tps, tpn, r, view(λps, :, j), view(λpn, :, j), s, w, τ)
        for i in 1:n_dac
            denom = tps[i] + tpn[i]
            if denom == 0.0
                sampled_pn[i, j] = 0.0
                sampled_ps[i, j] = 0.0
            else
                ov = observed_values[i]
                inv_denom = 1.0 / denom
                sampled_ps[i, j] = pois_rand(tps[i] * inv_denom * ov)
                sampled_pn[i, j] = pois_rand(tpn[i] * inv_denom * ov)
            end
        end
    end
    return nothing
end

"Sampler description for `provenance.txt`."
_sampling_desc(::AbstractPoissonSampling, param) = ""
_sampling_desc(s::RjSampling, param) =
    "  w ~ U($(s.w_lo), $(s.w_hi)), anchor freq $(s.anchor_lo)-$(s.anchor_hi)" *
    " (DAC $(ceil(Int, s.anchor_lo * param.nn))-$(floor(Int, s.anchor_hi * param.nn)))" *
    (s.tilt > 0 ? ", tilt τ ~ U(-$(s.tilt), $(s.tilt))" : "")

################################
######    Row priors      ######
################################

"""
    AbstractRowPrior

Which table rows `summary_statistics` draws, i.e. the ABC prior over the model parameters:
[`FlatAlphaPrior`](@ref) (the default) or [`TableRowPrior`](@ref) (uniform over the rows,
ABC-MK v1.0.0).
"""
abstract type AbstractRowPrior end

"""
    TableRowPrior()

Row prior of ABC-MK v1.0.0: table rows are drawn uniformly with replacement, so the ABC
prior is the design of the table built by [`rates`](@ref). Pass it as `row_prior` to
[`summary_statistics`](@ref) or [`abcmk`](@ref).
"""
struct TableRowPrior <: AbstractRowPrior end

"""
    FlatAlphaPrior(; bin = 0.01, floor = 0.001, by_B = true)

Default row prior of [`summary_statistics`](@ref) and [`abcmk`](@ref): table rows drawn
with replacement, weighted so that the prior on each row's realised α,
`(λW + λS) / λdn`, is flat within each `B` value of the table.

A row's weight is `π_B / max(n_bin, floor * n_B)`, where `n_bin` counts the rows of its
α bin (width `bin`) at its `B`, `n_B` the rows at that `B`, and
`π_B = n_B / (rows * occupied α bins at B)`. Every α bin of a `B` value then carries the
same weight and every `B` keeps its share of rows. `floor` stops the sparse edges of the
α range from taking huge weights. With `by_B = false`, α is flattened over the whole table
(weight `1 / max(n_bin, floor * rows)`), which also moves prior mass toward high `B`.
Rows with a non-finite α get weight 0.

Why: a table drawn with uniform α labels has a realised-α prior that falls with α,
because the realised α shrinks with `B`; rejection ABC then pulls α toward 0 wherever the
likelihood is broad.

Throws an `ArgumentError` unless `bin > 0` and `0 <= floor < 1`.
"""
struct FlatAlphaPrior <: AbstractRowPrior
    bin::Float64
    floor::Float64
    by_B::Bool
end

function FlatAlphaPrior(; bin = 0.01, floor = 0.001, by_B = true)
    (isfinite(bin) && bin > 0) || throw(ArgumentError("bin must be > 0, got $bin"))
    (isfinite(floor) && 0 <= floor < 1) || throw(ArgumentError("floor must be in [0, 1), got $floor"))
    return FlatAlphaPrior(Float64(bin), Float64(floor), Bool(by_B))
end

"Row weights over the candidate rows `idx` (nothing = uniform). `B`: the table's B column."
_row_weights(::TableRowPrior, dsdn::AbstractMatrix, idx; B = nothing) = nothing

function _row_weights(p::FlatAlphaPrior, dsdn::AbstractMatrix, idx; B = nothing)
    p.by_B && B === nothing && throw(ArgumentError("FlatAlphaPrior(by_B = true) needs the B column"))
    key = Vector{Tuple{Float64,Int}}(undef, length(idx))
    cnt = Dict{Tuple{Float64,Int},Int}()
    nB  = Dict{Float64,Int}()
    @inbounds for (k, i) in enumerate(idx)
        α = (dsdn[i, 3] + dsdn[i, 4]) / dsdn[i, 2]
        g = p.by_B ? Float64(B[i]) : 0.0
        key[k] = (g, isfinite(α) ? floor(Int, α / p.bin) : -1)
        nB[g] = get(nB, g, 0) + 1
        key[k][2] >= 0 && (cnt[key[k]] = get(cnt, key[k], 0) + 1)
    end
    # by_B: each B value carries its share of rows, split equally over its occupied α bins
    nbin = Dict{Float64,Int}()
    for k in keys(cnt); nbin[k[1]] = get(nbin, k[1], 0) + 1; end
    share(g) = p.by_B ? nB[g] / (length(idx) * nbin[g]) : 1.0
    return [x[2] >= 0 ? share(x[1]) / max(cnt[x], p.floor * nB[x[1]]) : 0.0 for x in key]
end

"Draw `n` rows of `idx` under `prior`."
function _sample_rows(prior::AbstractRowPrior, dsdn::AbstractMatrix, idx, n::Int; B = nothing)
    w = _row_weights(prior, dsdn, idx; B)
    return w === nothing ? sample(idx, n; replace = true) : sample(idx, Weights(w), n; replace = true)
end

_row_prior_desc(::TableRowPrior) = "TableRowPrior (uniform rows)"
_row_prior_desc(p::FlatAlphaPrior) = "FlatAlphaPrior (bin $(p.bin), floor $(p.floor), by_B $(p.by_B))"

"""
    sampling_summaries!(buf, fs, d, neut, sel, dsdn, gn, sh, gL, gH, B, output, summstat_size, sampling)

Write one dataset's summaries (see `_assemble_summaries!`) to `output`; returns the number
of rows written.
"""
function sampling_summaries!(
    buf::SamplingBuffers,
    fs::AbstractVector,
    d::AbstractVector,
    neut::AbstractMatrix,
    sel::AbstractMatrix,
    dsdn::AbstractMatrix,
    gn::AbstractVector,
    sh::AbstractVector,
    gL::AbstractVector,
    gH::AbstractVector,
    B::AbstractVector,
    output::String,
    summstat_size::Int64,
    sampling::AbstractPoissonSampling=IndependentSampling(),
)
    out = _assemble_summaries!(buf, fs, d, neut, sel, dsdn, gn, sh, gL, gH, B, summstat_size, sampling)
    write_files(out, output)
    return size(out, 1)
end

"""
    _assemble_summaries!(buf, fs, d, neut, sel, dsdn, gn, sh, gL, gH, B, summstat_size,
                         sampling = IndependentSampling()) -> Matrix{Float64}

Draw one dataset's Poisson realisations for every sampled table row, assemble
`[α(3), ω(4, dropped if any is NaN), gn, sh, gL, gH, B, α(x)(n_dac)]`, drop the
non-finite rows and keep at most `summstat_size`: the matrix `sampling_summaries!`
writes. The in-memory ABC route (`abcmk(...; in_memory = true)`) uses it directly.
"""
function _assemble_summaries!(
    buf::SamplingBuffers,
    fs::AbstractVector,
    d::AbstractVector,
    neut::AbstractMatrix,
    sel::AbstractMatrix,
    dsdn::AbstractMatrix,
    gn::AbstractVector,
    sh::AbstractVector,
    gL::AbstractVector,
    gH::AbstractVector,
    B::AbstractVector,
    summstat_size::Int64,
    sampling::AbstractPoissonSampling=IndependentSampling(),
)
    N = buf.N
    n_dac = buf.n_dac

    λds      = view(dsdn, :, 1)
    λdn      = view(dsdn, :, 2)
    λweak    = view(dsdn, :, 3)
    λstrong  = view(dsdn, :, 4)

    poisson_fixation!(buf.α_mat, buf.ω_mat, buf.sampled_dn, buf.sampled_ds,
                      d, λds, λdn, λweak, λstrong, sampling)
    poisson_polymorphism!(buf.sampled_pn, buf.sampled_ps, fs, neut, sel, sampling)

    # α_summaries[i,k] = round(1 - (sampled_ds[i]/sampled_dn[i]) * (sampled_pn[k,i]/sampled_ps[k,i]), 5)
    α_summaries = buf.α_summaries
    @inbounds for k in 1:n_dac
        for i in 1:N
            ratio_d = buf.sampled_ds[i] / buf.sampled_dn[i]
            ratio_p = buf.sampled_pn[k, i] / buf.sampled_ps[k, i]
            α_summaries[i, k] = round(1.0 - ratio_d * ratio_p, digits=5)
        end
    end

    # Drop the ω block entirely if any ω value is NaN (ABCreg then reads P = 8).
    has_ω = !any(isnan, buf.ω_mat)
    n_cols = has_ω ? 12 + n_dac : 8 + n_dac

    # Assemble expected_values column-by-column into preallocated buffer.
    #   [α(3), ω(4 if has_ω), gn, sh, gL, gH, B, α_summaries(n_dac)]
    ev = view(buf.expected_values, 1:N, 1:n_cols)
    col = 1
    @inbounds for k in 1:3
        for i in 1:N; ev[i, col] = round(buf.α_mat[i, k], digits=5); end
        col += 1
    end
    if has_ω
        @inbounds for k in 1:4
            for i in 1:N; ev[i, col] = round(buf.ω_mat[i, k], digits=5); end
            col += 1
        end
    end
    @inbounds for i in 1:N; ev[i, col] = gn[i]; end; col += 1
    @inbounds for i in 1:N; ev[i, col] = sh[i]; end; col += 1
    @inbounds for i in 1:N; ev[i, col] = gL[i]; end; col += 1
    @inbounds for i in 1:N; ev[i, col] = gH[i]; end; col += 1
    @inbounds for i in 1:N; ev[i, col] = B[i];  end; col += 1
    @inbounds for k in 1:n_dac
        for i in 1:N; ev[i, col] = α_summaries[i, k]; end
        col += 1
    end

    # Drop only rows with non-finite values (ABCreg fails on them). Rows are not filtered
    # on the sampled α: that would bias the regression adjustment toward (0, 1).
    n_kept = 0
    @inbounds for i in 1:N
        keep = true
        for c in 1:n_cols
            if !isfinite(ev[i, c])
                keep = false; break
            end
        end
        buf.mask[i] = keep
        n_kept += keep
    end

    # Subsample if more rows survive than requested: keep the first summstat_size kept rows.
    if n_kept > summstat_size
        target = summstat_size
        kept_seen = 0
        @inbounds for i in 1:N
            if buf.mask[i]
                kept_seen += 1
                if kept_seen > target
                    buf.mask[i] = false
                end
            end
        end
        n_kept = target
    end

    # Materialise the kept rows (single allocation).
    out = Matrix{Float64}(undef, n_kept, n_cols)
    row = 0
    @inbounds for i in 1:N
        if buf.mask[i]
            row += 1
            for c in 1:n_cols
                out[row, c] = ev[i, c]
            end
        end
    end

    return out
end

"""
    summary_statistics(param, sfs, divergence; h5_file, summstat_size, output_folder,
                       alpha = nothing, B_bins = nothing, sampling = RjSampling(),
                       row_prior = FlatAlphaPrior())

Write the ABC summary statistics of each dataset: draw `summstat_size` rows of the
reference table and turn each into the Poisson realisation of the data under that model.

`sfs` and `divergence` are vectors with one dataset each, as returned by `parse_sfs`
(a single matrix and vector are also accepted). The SFS must be RAW, per frequency
`[DAC Pn Ps]`: the tails are built here, and a cumulative SFS would be cumulated twice.
Use [`valid_data`](@ref) first to drop the datasets that would stop the run.

For dataset `i`, `output_folder` receives `summstat_i.txt.gz` (one row per sampled
table row: `α_weak, α_strong, α`, then `ωₐ_weak, ωₐ_strong, ωₐ, ωₙₐ` unless any is
NaN, then `|gam_dfe|, shape, gL, gH, B x 1000`, then α(x) at each `param.dac`) and
`alphas_i.txt.gz` (the observed α(x)), plus a `provenance.txt`. Earlier `summstat_*`,
`alphas_*` and `out_*` files in `output_folder` are deleted.

# Keywords
 - `h5_file::String`: reference table written by [`rates`](@ref) (JLD2, group
   `"N/n/cutoff"` of `param`).
 - `summstat_size::Real`: rows per dataset, a positive whole number (`1e5` works).
 - `output_folder::String`: created if missing; use one folder per analysis.
 - `alpha::Union{Nothing,Vector{Float64}} = nothing`: keep only rows with `al_tot` in
   `[alpha[1], alpha[end]]`. `al_tot` is the model's α label, which equals the realised
   fixation α only at `B = 0.999`.
 - `B_bins::Union{Nothing,Vector{Float64}} = nothing`: keep only rows with `B` in
   `[B_bins[1], B_bins[end]]`.
 - `sampling::AbstractPoissonSampling = RjSampling()`: [`RjSampling`](@ref) or
   [`IndependentSampling`](@ref) (v1.0.0).
 - `row_prior::AbstractRowPrior = FlatAlphaPrior()`: [`FlatAlphaPrior`](@ref) or
   [`TableRowPrior`](@ref) (v1.0.0).

# Returns
The paths of the `alphas_*` files. Throws if `param.dac` contains 1, if a dataset has
no polymorphism at some DAC or no divergence, or if fewer than `summstat_size` rows pass
the filters.
"""
function summary_statistics(
    param::parameters,
    sfs::Vector,
    divergence::Vector;
    h5_file::String,
    summstat_size::Real,
    output_folder::String,
    alpha::Union{Nothing,Vector{Float64}} = nothing,
    B_bins::Union{Nothing,Vector{Float64}} = nothing,
    sampling::AbstractPoissonSampling = RjSampling(),
    row_prior::AbstractRowPrior = FlatAlphaPrior(),
)

    @assert 1 ∉ param.dac "Please exclude singletons from the inference"
    summstat_size = _as_count(summstat_size, "summstat_size")

    ## Opening files
    assertion_params(param)
    mkpath(output_folder)

    # Remove any stale summstat_* / alphas_* / out_* left in output_folder by a previous
    # call. Without this, the ABC step globs the folder and pairs stale files (potentially
    # with a different DAC count) with the new ones.
    for fn in readdir(output_folder; join=true)
        bn = basename(fn)
        if startswith(bn, "summstat_") || startswith(bn, "alphas_") || startswith(bn, "out_")
            rm(fn)
        end
    end

    st = _summaries_setup(param, sfs, divergence; h5_file, summstat_size, alpha, B_bins,
                          sampling, row_prior, provenance_folder = output_folder)

    # Making summaries
    summ_output = output_folder .* "/summstat_" .* string.(1:length(st.sfs_p)) .* ".txt.gz"

    @info "Estimating and writting summary statistics at $output_folder"

    # One SamplingBuffers per dataset; the threaded loop only mutates its own buffer.
    ThreadsX.map(
        (i, z) -> sampling_summaries!(_dataset_buffer(st), st.sfs_p[i], st.divergence_p[i],
                                      st.neut, st.sel, st.dsdn, st.gn, st.sh, st.gL, st.gH,
                                      st.B, z, st.summstat_size, st.samplers[i]),
        eachindex(st.sfs_p),
        summ_output,
    );

    α_output = output_folder * "/alphas_" .* string.(1:length(st.sfs_p)) .* ".txt.gz"

    @. write_files(st.α, α_output)
end

"""
    _summaries_setup(param, sfs, divergence; h5_file, summstat_size, alpha, B_bins,
                     sampling, row_prior, provenance_folder = nothing)

Everything `summary_statistics` does before the per-dataset Poisson draws, shared with the
in-memory ABC route: read the rates table, filter and draw its rows under `row_prior`, lay
out the sampled spectra, and prepare one sampler and the observed α(x) per dataset. The
row draw is the only random step here, so both routes consume the global RNG identically.
Writes `provenance.txt` into `provenance_folder` when one is given.
"""
function _summaries_setup(
    param::parameters,
    sfs::Vector,
    divergence::Vector;
    h5_file::String,
    summstat_size::Int,
    alpha::Union{Nothing,Vector{Float64}},
    B_bins::Union{Nothing,Vector{Float64}},
    sampling::AbstractPoissonSampling,
    row_prior::AbstractRowPrior,
    provenance_folder::Union{Nothing,String} = nothing,
)
    @assert 1 ∉ param.dac "Please exclude singletons from the inference"
    assertion_params(param)

    @info "Opening and random sampling expected SFS and fixation rates"
    tmp = _read_rates(h5_file, param)
    @info "rates built with: " * tmp.options

    # Record how the summaries were made. `Base.open`, since this module imports
    # `GZip.open`.
    isnothing(provenance_folder) || Base.open(joinpath(provenance_folder, "provenance.txt"), "w") do io
        println(io, "rates_file  = ", h5_file)
        println(io, "rates_opts  = ", tmp.options)
        println(io, "sampling    = ", nameof(typeof(sampling)), _sampling_desc(sampling, param))
        println(io, "row_prior   = ", _row_prior_desc(row_prior))
        println(io, "S           = ", length(param.dac))
        println(io, "N, n        = ", param.N, ", ", param.n)
        println(io, "cutoff      = ", param.cutoff)
        println(io, "dac         = ", param.dac)
        println(io, "summstat    = ", summstat_size)
        println(io, "alpha_filt  = ", alpha, "  (on al_tot)")
        println(io, "B_bins_filt = ", B_bins)
        println(io, "written     = ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()))
    end

    α, sfs_p, divergence_p = data_to_poisson(sfs, divergence, param.dac)

    # Stop execution if any of the SFS contains 0 values
    # Low polymorphism distribution will cause ABC errors due to Inf values on α estimations
    if any(0 .∈ sfs_p) | any(0 .∈ getindex.(divergence_p,1))

        throw(
            ArgumentError(
                "Your SFS contains 0 values at the selected DACs or the divergence is 0. Please consider to bin the SFS and re-estimate the rates using the selected bin as sample the new sample size.",
            ),
        )
    end

    # Sample index to subset models
    idx = trues(size(tmp.models, 1))
    idx_alpha = trues(length(idx))
    idx_B = trues(length(idx))

    if !isnothing(alpha)
        # `al_tot` is the alpha LABEL: the realised fixation alpha only at B = 0.999.
        a_col = tmp.models.al_tot
        idx_alpha = a_col .>= alpha[1] .&& a_col .<= alpha[end]
    end

    if !isnothing(B_bins)
        idx_B = tmp.models.B .>= B_bins[1].*1000 .&& tmp.models.B .<= B_bins[end].*1000
    end

    idx = findall(idx .&& idx_alpha .&& idx_B)

    @assert length(idx) >= summstat_size "Check the filters or increase the number of rates solutions"

    # Increasing manually the summstat_size to ensure we retrieved the input summaries size
    increased_summstat = 1 * summstat_size

    sampled_idx = _sample_rows(row_prior, tmp.dsdn, idx, increased_summstat; B = tmp.models.B)

    @assert length(sampled_idx) >=  increased_summstat "Check the filters or increase the number of rates solutions"

    models_full = Matrix{Float64}(tmp.models)
    models = models_full[sampled_idx, :]
    dsdn   = tmp.dsdn[sampled_idx, :]

    # Pre-transposed SFS layout: (n_dac × N).
    n_dac = length(param.dac)
    N = increased_summstat
    neut = Matrix{Float64}(undef, n_dac, N)
    sel  = Matrix{Float64}(undef, n_dac, N)
    @inbounds for (i, v) in enumerate(param.dac)
        neut_col = tmp.neut[v]
        sel_col  = tmp.sel[v]
        for j in 1:N
            neut[i, j] = neut_col[sampled_idx[j]]
            sel[i,  j] = sel_col[sampled_idx[j]]
        end
    end

    # Model columns by position, the same in old (8 columns) and new (10) rates files:
    # B, al_low, al_tot, gam_flanking, gL, gH, shape, gam_dfe[, rho, ne_ratio].
    gn = abs.(view(models, :, 8))                  # |gam_dfe|
    sh = round.(view(models, :, 7), digits=5)      # shape (rounded)
    gL = ceil.(abs.(view(models, :, 5)))
    gH = ceil.(abs.(view(models, :, 6)))
    B  = collect(view(models, :, 1))               # B * 1000

    # One sampler per dataset: RjSampling needs each dataset's observed Ps.
    samplers = [_prepare_sampling(sampling, sfs[i], param.dac, param.nn) for i in eachindex(sfs)]

    return (; α, sfs_p, divergence_p, samplers, neut, sel, dsdn, gn, sh, gL, gH, B,
            N, n_dac, summstat_size)
end

"Fresh `SamplingBuffers` for one dataset of a `_summaries_setup` result."
_dataset_buffer(st) = SamplingBuffers(st.N, st.n_dac)

# One set, as `parse_sfs` returns it for a single gene list.
summary_statistics(param::parameters, sfs::AbstractMatrix{<:Real}, divergence::AbstractVector{<:Real}; kwargs...) =
    summary_statistics(param, [Matrix{Float64}(sfs)], [Vector{Float64}(divergence)]; kwargs...)

function pol_correction!(
    sfs_all::Vector{Matrix{Float64}},
    sfs_in::Vector{Matrix{Float64}};
    column::Vector{Int} = [2],
)
    pn_all = map(x -> sum(x[:, column]), sfs_all)
    sfs_pn = map(x -> sum(x[:, column]), sfs_in)
    ratio_all_sub = map((x, y) -> x / y, pn_all, sfs_pn)

    map((x, y) -> y[:, column] .= x * y[:, column], ratio_all_sub, sfs_in)
end

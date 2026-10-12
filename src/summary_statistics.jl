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
statistics: [`IndependentSampling`](@ref) (the default) or [`RjSampling`](@ref), which
accounts for demography in the polymorphism.
"""
abstract type AbstractPoissonSampling end

"""
    IndependentSampling()

Every class drawn as an independent Poisson. The default.

Let `O_k` be the observed polymorphisms at DAC >= k (Pn + Ps) and `D` the observed
divergence. Per table row, with tail rates `λps_k`, `λpn_k` and fixation rates `λds`,
`λdn`, `λW`, `λS`:

    ps_k ~ Poisson(λps_k / (λps_k + λpn_k) · O_k)    pn_k ~ Poisson(λpn_k / (λps_k + λpn_k) · O_k)
    ds   ~ Poisson(λds / (λds + λdn) · D)            dn   ~ Poisson(λdn / (λds + λdn) · D)

`sw + ss > dn` is possible, so `α_total > 1` occurs as a sampling artifact; such rows are
kept (only non-finite rows are dropped).
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
    RjSampling(; w = (0.0, 1.5), anchor = (0.15, 0.70))

`IndependentSampling` with a demographic nuisance on the polymorphism composition:
GRAPES' r_j, normalised on an anchor band, with the nonsynonymous class distorted by
r_j^w.

Per observed dataset, the synonymous distortion per interval `i` between consecutive DACs
of `param.dac` (the last interval runs to the top of the spectrum) is measured against each
table row's own neutral expectation:

    r_i = (Ps_obs,i / sum Ps_obs) / (λps_i / sum λps)

normalised so that its λps-weighted mean over the ANCHOR band is 1: the intervals lying
inside `anchor`, a band of derived-allele FREQUENCY, `anchor[1] <= k/nn` and
`stop/nn <= anchor[2]` with `nn = 2n` sampled chromosomes -- so the same band works at any
sample size. Each row then draws `w ~ Uniform(w...)` and its tail composition becomes

    P_k = sum_{i>=k} λps_i r_i / (sum_{i>=k} λps_i r_i + sum_{i>=k} λpn_i r_i^w)

after which `ps_k ~ Poisson(P_k O_k)`, `pn_k ~ Poisson((1 - P_k) O_k)` as in
`IndependentSampling`; fixations are `IndependentSampling`'s. `w` is marginalised: it is
not written to the summaries, so the file layout and ABCreg's `P` are unchanged.

- `r ≡ 1` (the data's synonymous shape equals the model's) gives `IndependentSampling`'s
  composition for any `w`; `w = 1` gives GRAPES-style conditioning on the observed Ps.
- Why an anchor: rescaling `r` by `c` multiplies Pn/Ps by `c^(w-1)`, a pure level that only
  α may set. The band is where `r = 1` means "not distorted"; the default 15-70% of the
  sample (DAC 199-925 at n = 661; DAC 15-70 at n = 50) assumes the distortion is RECENT and
  stays below the 70% cut that empirical data need against mispolarization. `param.dac`
  must place at least one whole interval inside it (e.g. at n = 50, DACs 20 and 50 give
  the intervals 20-49 and 50-69). Plain r_j (no band) was measured worse.
- Requirements: `param.dac` ascending, Ps > 0 in every DAC interval, at least one whole
  interval inside the anchor band. Each violation throws an `ArgumentError`.
- Measured against SLiM on the forward model `rates()` builds (ρ drawn per model; Tennessen
  history, 661 diploids, 6 ρ x B cells; two-epoch histories, 50 diploids, 8 cells; ABCreg,
  tol 1000 of 1e5): α RMSE 0.028 vs 0.040 for `IndependentSampling` on the ρ grid, 0.047
  vs 0.068 on the two-epoch grid, α coverage 0.95-0.96. One α and one DFE shape per grid,
  one seed per run, bootstraps of one simulated pool, no 70% cut applied.
- WARNING: pass the RAW per-frequency SFS. `summary_statistics` builds the tails itself;
  a spectrum cumulated twice flattens α(x) and biases α.
- WARNING: not for contractions. A bottleneck distorts the anchor band itself and leaves a
  flat deleterious Pn/Ps excess that no `r^w` produces; the anchor assigns it to α (α
  biased by about -0.10 under a 20x bottleneck at B = 0.999).
"""
struct RjSampling <: AbstractPoissonSampling
    w_lo::Float64
    w_hi::Float64
    anchor_lo::Float64        # anchor band, as a frequency (fraction of the 2n chromosomes)
    anchor_hi::Float64
    ps_obs::Vector{Float64}   # observed Ps per interval (empty until _prepare_sampling)
    band::Vector{Bool}        # intervals inside the anchor band
end

function RjSampling(; w = (0.0, 1.5), anchor = (0.15, 0.70))
    lo, hi = Float64(first(w)), Float64(last(w))
    (isfinite(lo) && isfinite(hi) && 0 <= lo <= hi) ||
        throw(ArgumentError("w must satisfy 0 <= w[1] <= w[2] < Inf, got $w"))
    a_lo, a_hi = Float64(first(anchor)), Float64(last(anchor))
    0 < a_lo < a_hi <= 1 || throw(ArgumentError(
        "anchor must be a frequency band (lo, hi) with 0 < lo < hi <= 1, got $anchor"))
    return RjSampling(lo, hi, a_lo, a_hi, Float64[], Bool[])
end

function _prepare_sampling(s::RjSampling, sfs::AbstractMatrix,
                           dac::Vector{Int64}, nn::Integer)
    issorted(dac) || throw(ArgumentError("RjSampling needs param.dac ascending"))
    cum = cumulative_sfs(sfs)
    tail_ps = map(z -> sum(cum[cum[:, 1] .== z, 3]), dac)
    n = length(dac)
    ps_obs = [i < n ? tail_ps[i] - tail_ps[i+1] : tail_ps[n] for i in 1:n]
    stops = [i < n ? dac[i+1] - 1 : nn - 1 for i in 1:n]
    # the band is a frequency; resolve it with this dataset's sample size
    band = [dac[i] >= s.anchor_lo * nn && stops[i] <= s.anchor_hi * nn for i in 1:n]
    any(band) || throw(ArgumentError(
        "no interval of dac = $dac lies inside the anchor band $(s.anchor_lo)-$(s.anchor_hi) " *
        "(DAC $(ceil(Int, s.anchor_lo * nn))-$(floor(Int, s.anchor_hi * nn)) at 2n = $nn): " *
        "add DACs so that one interval [dac[i], dac[i+1] - 1] falls inside it"))
    all(>(0), ps_obs) || throw(ArgumentError(
        "RjSampling needs Ps > 0 in every DAC interval; got $ps_obs for dac = $dac"))
    return RjSampling(s.w_lo, s.w_hi, s.anchor_lo, s.anchor_hi, ps_obs, band)
end

"""
    _rj_tails!(tps, tpn, r, λps_col, λpn_col, s, w)

Distorted tail sums of one table row for [`RjSampling`](@ref): fills
`tps[k] = sum_{i>=k} λps_i r_i`, `tpn[k] = sum_{i>=k} λpn_i r_i^w` from the row's TAIL sums
`λps_col`, `λpn_col` (split into intervals by differencing), with `r` normalised on the
anchor band. `r` is scratch.
"""
function _rj_tails!(tps, tpn, r, λps_col, λpn_col, s::RjSampling, w::Float64)
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
    acc_s = 0.0; acc_n = 0.0
    @inbounds for i in n:-1:1
        e_s = max(0.0, i < n ? λps_col[i] - λps_col[i+1] : λps_col[n])
        e_n = max(0.0, i < n ? λpn_col[i] - λpn_col[i+1] : λpn_col[n])
        ri = r[i] / c
        acc_s += e_s * ri
        acc_n += e_n * ri^w
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
    @inbounds for j in 1:N
        w = fixed_w ? s.w_lo : rand(u)
        _rj_tails!(tps, tpn, r, view(λps, :, j), view(λpn, :, j), s, w)
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
    " (DAC $(ceil(Int, s.anchor_lo * param.nn))-$(floor(Int, s.anchor_hi * param.nn)))"

################################
######    Row priors      ######
################################

"""
    AbstractRowPrior

Which table rows `summary_statistics` draws, i.e. the ABC prior over the model parameters:
[`TableRowPrior`](@ref) (the default, uniform over the rows) or [`FlatAlphaPrior`](@ref).
"""
abstract type AbstractRowPrior end

"""
    TableRowPrior()

Rows drawn uniformly (with replacement): the prior is the table's own design. The default.
"""
struct TableRowPrior <: AbstractRowPrior end

"""
    FlatAlphaPrior(; bin = 0.01, floor = 0.001, by_B = false)

Rows drawn with replacement, weighted so that the prior on each row's realised α
`(λW + λS) / λdn` is flat: weight `1 / max(count of its α bin, floor · rows)`, bins of
width `bin`. The floor keeps the sparse edges of the α range from taking huge weights.
With `by_B = true` α is flattened within each B value of the table instead (counts and
floor per B), so the prior on B keeps the table's uniform design; otherwise flattening α
also moves prior mass toward high B, where the high α values are.

Why: a table drawn with uniform α LABELS has a realised-α prior that falls with α, because
the realised α shrinks with B (two-epoch table: 18.5% of rows below 0.05, 1.7% in
[0.45, 0.50)). Rejection ABC then pulls α toward 0 wherever the likelihood is broad.
Measured in the offline lab (`results/rj/lab*`, 2026-10-11): bottleneck α RMSE 0.105 -> 0.039
(bt_strong, B .999) and 0.038 -> 0.019 (bt_moderate, B 0.4), ρ ladder unchanged (|Δ| ≤ 0.001).
"""
struct FlatAlphaPrior <: AbstractRowPrior
    bin::Float64
    floor::Float64
    by_B::Bool
end

function FlatAlphaPrior(; bin = 0.01, floor = 0.001, by_B = false)
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

Draw one dataset's Poisson realisations for every sampled table row, assemble
`[α(3), ω(4, dropped if any is NaN), gn, sh, gL, gH, B, α(x)(n_dac)]`, drop the
non-finite rows, keep at most `summstat_size`, and write them to `output`.
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

    # Materialise the kept rows (single allocation), then write.
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

    write_files(out, output)
    return n_kept
end

"""
    summary_statistics(param,sfs,divergence;h5_file,summstat_size,output_folder,...)

Estimate summary statistics using observed data and analytical rates. *output_folder* will be used to output summary statistics

# Arguments
 - `param::parameters` : Mutable structure containing the models
 - `h5_file::String` : JLD2 file containing solved models, fixation and polymorphic rates
 - `sfs::Vector`: SFS data parsed from parse_sfs(). The RAW per-frequency SFS: the tails
   are built here, and a spectrum cumulated twice flattens α(x) and biases α.
 - `divergence::Vector`: divergence data from parse_sfs().
 - `output_folder::String` : path to save summary statistics and observed data.
 - `summstat_size::Int64` : number of summary statistics to sample.
 - `sampling::AbstractPoissonSampling` : how the Poisson realisations are drawn.
   `IndependentSampling()` (default) or `RjSampling()`.
 - `alpha::Union{Nothing,Vector{Float64}}` : keep only models with `al_tot` in
   `[alpha[1], alpha[end]]`. `al_tot` is the model's alpha LABEL, which equals the realised
   fixation alpha only in the B = 0.999 bin.
 - `B_bins::Union{Nothing,Vector{Float64}}` : keep only models with B in this range.
 - `row_prior::AbstractRowPrior` : how table rows are drawn. `TableRowPrior()` (default)
   or `FlatAlphaPrior()`.
"""
function summary_statistics(
    param::parameters,
    sfs::Vector,
    divergence::Vector;
    h5_file::String,
    summstat_size::Int64,
    output_folder::String,
    alpha::Union{Nothing,Vector{Float64}} = nothing,
    B_bins::Union{Nothing,Vector{Float64}} = nothing,
    sampling::AbstractPoissonSampling = IndependentSampling(),
    row_prior::AbstractRowPrior = TableRowPrior(),
)

    @assert 1 ∉ param.dac "Please exclude singletons from the inference"

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

    @info "Opening and random sampling expected SFS and fixation rates"
    tmp = _read_rates(h5_file, param)
    @info "rates built with: " * tmp.options

    # Record how the summaries were made. `Base.open`, since this module imports
    # `GZip.open`.
    Base.open(joinpath(output_folder, "provenance.txt"), "w") do io
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

    # Making summaries
    summ_output = output_folder .* "/summstat_" .* string.(1:size(sfs_p, 1)) .* ".txt"

    @info "Estimating and writting summary statistics at $output_folder"

    # One sampler per dataset: RjSampling needs each dataset's observed Ps.
    samplers = [_prepare_sampling(sampling, sfs[i], param.dac, param.nn) for i in eachindex(sfs)]

    # One SamplingBuffers per dataset; the threaded loop only mutates its own buffer.
    ThreadsX.map(
        (x, y, z, smp) -> begin
            buf = SamplingBuffers(N, n_dac)
            sampling_summaries!(buf, x, y, neut, sel, dsdn,
                                gn, sh, gL, gH, B, z, summstat_size, smp)
        end,
        sfs_p,
        divergence_p,
        summ_output,
        samplers
    );

    α_output = output_folder * "/alphas_" .* string.(1:size(sfs_p, 1)) .* ".txt"

    @. write_files(α, α_output)
end

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

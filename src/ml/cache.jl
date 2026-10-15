# =============================================================================
# Precomputed tables and the cache that owns them
#
# There is no global cache. `DFECache` is an ordinary object: `dfe_ml` makes a
# fresh one unless the caller passes one in, and `dfe_model_comparison` shares one
# across its models. Every table that is filled lazily lives inside the
# `PrecomputedData` that needs it, behind that object's lock.
# =============================================================================

function precalculate_expected_counts(s_grid::Vector{Float64}, n::Int,
                                      integ::AbstractIntegrator)::Matrix{Float64}
    n_s = length(s_grid)
    n_bins = n - 1
    counts = Matrix{Float64}(undef, n_s, n_bins)

    logfact = zeros(Float64, n + 1)
    for i in 2:(n+1)
        logfact[i] = logfact[i-1] + log(i - 1)
    end
    logc = Vector{Float64}(undef, n_bins)
    for j in 1:n_bins
        logc[j] = logfact[n+1] - logfact[j+1] - logfact[n-j+1]
    end

    # Each iteration writes one row and reads nothing another writes, so the values
    # do not depend on the thread count.
    Threads.@threads for iS in eachindex(s_grid)
        S = s_grid[iS]
        @inbounds for j in 1:n_bins
            if abs(S) < QUASI_ZERO
                counts[iS, j] = 1.0 / j
            else
                counts[iS, j] = expected_count(integ, S, j, n, logc[j])
            end
        end
    end
    return counts
end

expectation_counts(::PRFSojourn, s_eff::Vector{Float64}, n::Int, integ::AbstractIntegrator) =
    precalculate_expected_counts(s_eff, n, integ)

"""
    DnCorrectionTable(n, integ)

GRAPES' dN correction for sample size `n` at integer |S| (0..500 below zero,
0..10000 above), computed with `integ` on first use and kept. Thread safe.
"""
struct DnCorrectionTable{I<:AbstractIntegrator}
    n::Int
    integ::I
    neg::Vector{Float64}   # index k + 1, NaN until computed
    pos::Vector{Float64}
    lock::ReentrantLock
end
DnCorrectionTable(n::Int, integ::AbstractIntegrator) =
    DnCorrectionTable(n, integ, fill(NaN, Int(LARGE_S_VALUE) + 1),
                      fill(NaN, Int(VERY_LARGE_S_VALUE) + 1), ReentrantLock())

"""
    PrecomputedData

Everything that depends on the sample size and on the integrator, S grid, dominance
and SFS expectation, but not on the data or the DFE: the expected-counts table
(`length(s_grid)` x `n-1`), fixation probabilities, cell bounds, and two lazily
filled tables (GRAPES' dN correction by integer |S|, and the expected-counts row of
each discrete DFE component). Built by `get_precomputed`; safe to share across
threads (the lazy tables are filled under `lock`).
"""
struct PrecomputedData{I<:AbstractIntegrator,E<:AbstractSFSExpectation,D<:AbstractDominance}
    n::Int
    s_grid::Vector{Float64}          # S of the DFE (real coefficients)
    s_widths::Vector{Float64}
    s_eff::Vector{Float64}           # effective S under dominance (== s_grid when additive)
    expected_counts::Matrix{Float64} # (n_s, n-1)
    fixprob::Vector{Float64}
    fixprob_weighted::Vector{Float64}
    neutral_shape::Vector{Float64}
    integ::I
    expectation::E
    dominance::D
    s_lo::Vector{Float64}
    s_hi::Vector{Float64}
    s_ext_lo::Float64
    s_ext_hi::Float64
    s_edges::Vector{Float64}
    lo_idx::Vector{Int}
    hi_idx::Vector{Int}
    dn::DnCorrectionTable{I}
    discrete_rows::Dict{Float64,Vector{Float64}}
    lock::ReentrantLock
end

"""
    DFECache()

Owner of the `PrecomputedData` tables, keyed by value on
`(n, integrator, s_grid, dominance, expectation)`. Pass one to several `dfe_ml`
calls (`cache = c`) to reuse the tables across datasets of the same sample size;
safe to share across threads.
"""
struct DFECache
    tables::Dict{Any,PrecomputedData}
    lock::ReentrantLock
end
DFECache() = DFECache(Dict{Any,PrecomputedData}(), ReentrantLock())

Base.length(c::DFECache) = length(c.tables)

"""
    get_precomputed(cache, n; integrator, grid, dominance, expectation) -> PrecomputedData

The tables for sample size `n` (chromosomes), built on first use.
"""
function get_precomputed(cache::DFECache, n::Int;
                         integrator::AbstractIntegrator = ExactSojourn(),
                         grid::Union{Nothing,AbstractSGrid} = nothing,
                         dominance::AbstractDominance = Additive(),
                         expectation::AbstractSFSExpectation = PRFSojourn())
    g = grid === nothing ? default_s_grid(expectation) : grid
    key = (n, integrator, g, dominance, expectation)
    lock(cache.lock) do
        get!(cache.tables, key) do
            _build_precomputed(n, integrator, g, dominance, expectation)
        end
    end
end

function _build_precomputed(n::Int, integ::AbstractIntegrator, g::AbstractSGrid,
                            dominance::AbstractDominance, expectation::AbstractSFSExpectation)
    n >= 3 || throw(ArgumentError("dfe_ml needs at least 3 chromosomes, got n = $n"))
    s_grid = s_grid_points(g)
    s_widths = s_grid_widths(g, s_grid)

    # Romberg cannot converge below S = -1500 (`exp(-S(1-x))` overflows and the
    # density is NaN), so it is refused up front instead of spinning to JMAX.
    if integ isa RombergIntegrator && minimum(s_grid) < -1500.0
        error("RombergIntegrator does not converge below S = -1500 " *
              "(got S_min = $(minimum(s_grid))); use ExactSojourn(), " *
              "GaussLegendreSojourn() or a grid with |S| <= 1500")
    end

    s_eff = [effective_s(dominance, S) for S in s_grid]
    expected_counts = expectation_counts(expectation, s_eff, n, integ)

    fixprob = similar(s_grid)
    fixprob_weighted = similar(s_grid)
    @inbounds for i in eachindex(s_grid)
        fp = fixation_probability(s_eff[i])
        fixprob[i] = fp
        fixprob_weighted[i] = fp * s_widths[i]
    end

    s_lo, s_hi = s_cell_bounds(g, s_grid, s_widths)
    ext_lo, ext_hi = _exterior_edges(s_grid, s_widths)
    edges, lo_idx, hi_idx = _cell_edge_index(s_lo, s_hi)
    return PrecomputedData(n, s_grid, s_widths, s_eff, expected_counts, fixprob,
                           fixprob_weighted, neutral_shape(expectation, n), integ,
                           expectation, dominance, s_lo, s_hi, ext_lo, ext_hi,
                           edges, lo_idx, hi_idx, DnCorrectionTable(n, integ),
                           Dict{Float64,Vector{Float64}}(), ReentrantLock())
end

function _dn_correction_integral(S::Float64, n::Int, integ::AbstractIntegrator)
    abs(S) <= QUASI_ZERO && return 1.0 / n
    return dn_correction_term(integ, S, n)
end

"""
    dn_correction(t::DnCorrectionTable, S) -> Float64
    dn_correction(pre::PrecomputedData, S) -> Float64

GRAPES' dN correction at `S`, looked up at the nearest integer |S| (clamped to 500
below zero and 10000 above, as GRAPES does) and computed on first use.
"""
function dn_correction(t::DnCorrectionTable, S::Float64)
    abs(S) <= QUASI_ZERO && return 1.0 / t.n
    k = Int(round(abs(S)))
    if S > 0
        k = min(k, Int(VERY_LARGE_S_VALUE))
        tab = t.pos
        Sk = Float64(k)
    else
        k = min(k, Int(LARGE_S_VALUE))
        tab = t.neg
        Sk = -Float64(k)
    end
    v = @inbounds tab[k + 1]
    isnan(v) || return v
    return _fill_dn!(t, tab, k, Sk)
end

dn_correction(pre::PrecomputedData, S::Float64) = dn_correction(pre.dn, S)

# Computed under the lock; a reader that raced a writer re-checks here.
function _fill_dn!(t::DnCorrectionTable, tab::Vector{Float64}, k::Int, Sk::Float64)
    lock(t.lock)
    try
        v = tab[k + 1]
        if isnan(v)
            v = _dn_correction_integral(Sk, t.n, t.integ)
            tab[k + 1] = v
        end
        return v
    finally
        unlock(t.lock)
    end
end

"""
    discrete_row(pre, S) -> Vector{Float64}

Expected counts (unfolded, length `n-1`) of a point mass at `S`, from the same
producer, integrator and dominance as `pre.expected_counts`. Computed once per S.
"""
function discrete_row(pre::PrecomputedData, S::Float64)
    lock(pre.lock)
    try
        return get!(pre.discrete_rows, S) do
            vec(expectation_counts(pre.expectation, [effective_s(pre.dominance, S)], pre.n, pre.integ))
        end
    finally
        unlock(pre.lock)
    end
end

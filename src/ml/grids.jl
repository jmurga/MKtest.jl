# =============================================================================
# S grids, cell bounds, and the DFE quadratures
# =============================================================================

function precalculate_s_grid()::Vector{Float64}
    n_total = 4 * S_GRID_K1
    s_grid = Vector{Float64}(undef, n_total)
    idx = 1

    # large negative: [-LARGE_S_VALUE, -GAMMA_INTEGRAL_BREAK]
    a1 = (LARGE_S_VALUE - GAMMA_INTEGRAL_BREAK) / S_GRID_K1
    for i in 0:(S_GRID_K1-1)
        s_grid[idx] = -LARGE_S_VALUE + (i + 0.5) * a1
        idx += 1
    end

    # medium: [-GAMMA_INTEGRAL_BREAK, +GAMMA_INTEGRAL_BREAK]
    a2 = 2.0 * GAMMA_INTEGRAL_BREAK / S_GRID_K2
    for i in 0:(S_GRID_K2-1)
        s_grid[idx] = -GAMMA_INTEGRAL_BREAK + (i + 0.5) * a2
        idx += 1
    end

    # large positive: [GAMMA_INTEGRAL_BREAK, LARGE_S_VALUE]
    for i in 0:(S_GRID_K1-1)
        s_grid[idx] = GAMMA_INTEGRAL_BREAK + (i + 0.5) * a1
        idx += 1
    end

    # very large positive: [LARGE_S_VALUE, VERY_LARGE_S_VALUE]
    a3 = (VERY_LARGE_S_VALUE - LARGE_S_VALUE) / S_GRID_K1
    for i in 0:(S_GRID_K1-1)
        s_grid[idx] = LARGE_S_VALUE + (i + 0.5) * a3
        idx += 1
    end

    return s_grid
end

function precalculate_s_widths(s_grid::Vector{Float64})::Vector{Float64}
    n = length(s_grid)
    widths = Vector{Float64}(undef, n)
    for i in 1:n
        if i == 1
            widths[i] = s_grid[i+1] - s_grid[i]
        elseif i == n
            widths[i] = s_grid[i] - s_grid[i-1]
        else
            widths[i] = (s_grid[i+1] - s_grid[i-1]) / 2.0
        end
    end
    return widths
end

s_grid_points(::GrapesGrid) = precalculate_s_grid()
s_grid_widths(::GrapesGrid, grid::Vector{Float64}) = precalculate_s_widths(grid)

function s_grid_points(g::LogGammaGrid)::Vector{Float64}
    lo, hi = g.int_bounds
    k = max(2, g.Npts ÷ 2)
    neg = -2 .* exp.(range(log(hi), log(lo); length = k))
    pos =  2 .* exp.(range(log(lo), log(VERY_LARGE_S_VALUE / 2); length = k))
    return vcat(neg, pos)
end

# Widths from geometric midpoints on each side of zero; the two ends are half-cells
# (the exterior terms carry the mass beyond them).
function s_grid_widths(::LogGammaGrid, grid::Vector{Float64})::Vector{Float64}
    n = length(grid)
    w = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        si = grid[i]
        same(j) = 1 <= j <= n && sign(grid[j]) == sign(si)
        lo = same(i - 1) ? sqrt(abs(si * grid[i-1])) : abs(si)
        hi = same(i + 1) ? sqrt(abs(si * grid[i+1])) : abs(si)
        w[i] = abs(hi - lo)
        w[i] == 0 && (w[i] = same(i - 1) ? abs(si - grid[i-1]) / 2 :
                             same(i + 1) ? abs(grid[i+1] - si) / 2 : 0.0)
    end
    return w
end

"""
    s_cell_bounds(g, grid, widths) -> (lo, hi)

Per-cell integration bounds, used by `IntervalCDF`. Pairs rather than one edge
vector, because cells need not tile (`LogGammaGrid` leaves a gap at S = 0, which the
exterior term covers). Default: `S_i -+ w_i/2`.
"""
function s_cell_bounds(::AbstractSGrid, s_grid::Vector{Float64}, s_widths::Vector{Float64})
    lo = similar(s_grid); hi = similar(s_grid)
    @inbounds for i in eachindex(s_grid)
        lo[i] = s_grid[i] - s_widths[i] / 2
        hi[i] = s_grid[i] + s_widths[i] / 2
    end
    return lo, hi
end

# GrapesGrid is four uniform blocks: its exact cells tile. (The central-difference
# widths overlap at the block boundaries, one of which is |S| = 0.5.)
function s_cell_bounds(::GrapesGrid, grid::Vector{Float64}, ::Vector{Float64})
    n = length(grid)
    lo = Vector{Float64}(undef, n); hi = Vector{Float64}(undef, n)
    a1 = (LARGE_S_VALUE - GAMMA_INTEGRAL_BREAK) / S_GRID_K1
    a2 = 2.0 * GAMMA_INTEGRAL_BREAK / S_GRID_K2
    a3 = (VERY_LARGE_S_VALUE - LARGE_S_VALUE) / S_GRID_K1
    blocks = ((-LARGE_S_VALUE, a1, S_GRID_K1),          # [-500, -0.5]
              (-GAMMA_INTEGRAL_BREAK, a2, S_GRID_K2),   # [-0.5, 0.5]
              (GAMMA_INTEGRAL_BREAK, a1, S_GRID_K1),    # [0.5, 500]
              (LARGE_S_VALUE, a3, S_GRID_K1))           # [500, 10000]
    idx = 1
    for (start, a, k) in blocks
        for i in 0:(k - 1)
            lo[idx] = start + i * a
            hi[idx] = start + (i + 1) * a
            idx += 1
        end
    end
    idx == n + 1 || error("GrapesGrid cell bounds: built $(idx-1) cells for $n points")
    return lo, hi
end

function s_cell_bounds(::LogGammaGrid, grid::Vector{Float64}, widths::Vector{Float64})
    n = length(grid)
    lo = similar(grid); hi = similar(grid)
    @inbounds for i in 1:n
        si = grid[i]
        same(j) = 1 <= j <= n && sign(grid[j]) == sign(si)
        a = same(i - 1) ? sign(si) * sqrt(abs(si * grid[i-1])) : si
        b = same(i + 1) ? sign(si) * sqrt(abs(si * grid[i+1])) : si
        lo[i], hi[i] = minmax(a, b)
        if hi[i] == lo[i]
            h = same(i - 1) ? abs(si - grid[i-1]) / 2 :
                same(i + 1) ? abs(grid[i+1] - si) / 2 : 0.0
            lo[i], hi[i] = si - h, si + h
        end
    end
    return lo, hi
end

# Distinct cell bounds, ascending, with `edges[lo_idx[i]] == s_lo[i]` exactly, so
# `IntervalCDF` evaluates the CDF once per edge and shared edges are bit-identical.
function _cell_edge_index(s_lo::Vector{Float64}, s_hi::Vector{Float64})
    edges = sort!(unique!(vcat(s_lo, s_hi)))
    lo_idx = [searchsortedfirst(edges, x) for x in s_lo]
    hi_idx = [searchsortedfirst(edges, x) for x in s_hi]
    return edges, lo_idx, hi_idx
end

# The two |S| edges outside which the deleterious DFE is handled analytically when
# `uses_exterior`: below `lo` effectively neutral, beyond `hi` effectively lethal.
# `lo` is 0 for a grid whose cells tile down to S = 0 (GrapesGrid).
function _exterior_edges(s_grid::Vector{Float64}, s_widths::Vector{Float64})
    inner = 0; outer = 0
    @inbounds for i in eachindex(s_grid)
        s_grid[i] < 0.0 || continue
        (inner == 0 || s_grid[i] > s_grid[inner]) && (inner = i)
        (outer == 0 || s_grid[i] < s_grid[outer]) && (outer = i)
    end
    inner == 0 && return (0.0, 0.0)
    lo = max(abs(s_grid[inner]) - s_widths[inner] / 2, 0.0)
    hi = abs(s_grid[outer]) + s_widths[outer] / 2
    return (lo, hi)
end

@inline function _last_index_leq(s_grid::Vector{Float64}, x::Float64)
    x == Inf && return length(s_grid)
    x == -Inf && return 0
    return searchsortedlast(s_grid, x)
end

# -----------------------------------------------------------------------------
# DFE quadratures
# -----------------------------------------------------------------------------

"""
    fill_dfe_vals!(out, quad, m, p, pre) -> out

The DFE's average density in each S cell of `pre`, which the engine multiplies by
`pre.s_widths`. The engine calls the 6-argument form with a scratch vector of length
`length(pre.s_edges)`; a user quadrature only needs the 5-argument one.
"""
fill_dfe_vals!(out, ::MidpointDensity, m::AbstractDFEModel, p, pre) =
    dfe_density!(out, m, p, pre.s_grid)

fill_dfe_vals!(out, q::AbstractDFEQuadrature, m::AbstractDFEModel, p, pre, scratch) =
    fill_dfe_vals!(out, q, m, p, pre)

fill_dfe_vals!(out, q::IntervalCDF, m::AbstractDFEModel, p, pre) =
    fill_dfe_vals!(out, q, m, p, pre, Vector{Float64}(undef, length(pre.s_edges)))

function fill_dfe_vals!(out, ::IntervalCDF, m::AbstractDFEModel, p, pre, phi)
    dfe_cdf_continuous!(phi, m, p, pre.s_edges)
    @inbounds for i in eachindex(out)
        w = pre.s_widths[i]
        # mass / width, so that the caller's `* s_widths[i]` recovers the mass
        out[i] = w > 0 ? max(phi[pre.hi_idx[i]] - phi[pre.lo_idx[i]], 0.0) / w : 0.0
    end
    return out
end

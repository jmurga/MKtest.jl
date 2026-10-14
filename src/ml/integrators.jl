# =============================================================================
# Integrators: the expected-SFS and dN-correction integrals over the PRF sojourn
#
#   expected_count(S, j, n) = 1/2 * int psi(x, S) C(n, j) x^j (1 - x)^(n-j) dx
#   dn_correction(S, n)     =       int psi(x, S)/2 * x^n dx
#
# with psi(x, S) = 2 (1 - e^{-S(1-x)}) / (x (1 - x) (1 - e^{-S})).
# =============================================================================

# Romberg (GRAPES `qromb`)
const EPS_R = 1.0e-6
const JMAX = 100
const K_ROMB = 5

function _polint(xa::Vector{Float64}, ya::Vector{Float64}, x::Float64)
    n = length(xa)
    c = copy(ya)
    d = copy(ya)
    ns = 1
    dif = abs(x - xa[1])

    for i in 1:n
        dift = abs(x - xa[i])
        if dift < dif
            ns = i
            dif = dift
        end
    end

    y = ya[ns]
    ns -= 1

    dy = 0.0
    for m in 1:(n-1)
        for i in 1:(n-m)
            ho = xa[i] - x
            hp = xa[i+m] - x
            w = c[i+1] - d[i]
            den = ho - hp
            if den == 0.0
                error("polint failure")
            end
            den = w / den
            d[i] = hp * den
            c[i] = ho * den
        end
        if 2 * ns < (n - m)
            dy = c[ns+1]
        else
            dy = d[ns]
            ns -= 1
        end
        y += dy
    end
    return y, dy
end

function _qromb(f, a::Float64, b::Float64)
    s = 0.0
    h = zeros(Float64, JMAX + 1)
    sarr = zeros(Float64, JMAX + 1)
    h[1] = 1.0

    for j in 1:JMAX
        if j == 1
            s = 0.5 * (b - a) * (f(a) + f(b))
        else
            it = 1 << (j - 2)
            tnm = Float64(it)
            del = (b - a) / tnm
            x = a + 0.5 * del
            acc = 0.0
            for _ in 1:it
                acc += f(x)
                x += del
            end
            s = 0.5 * (s + (b - a) * acc / tnm)
        end
        sarr[j] = s
        if j >= K_ROMB
            y, dy = _polint(h[j-K_ROMB+1:j], sarr[j-K_ROMB+1:j], 0.0)
            if abs(dy) < EPS_R * abs(y)
                return y
            end
        end
        h[j+1] = 0.25 * h[j]
    end
    error("Too many steps in qromb")
end

@inline function _sojourn_density(x::Float64, S::Float64)
    return 2.0 * (1.0 - exp(-S * (1.0 - x))) / (x * (1.0 - x) * (1.0 - exp(-S)))
end

"""
    fixation_probability(S) -> Float64

Fixation probability of a new mutation with coefficient `S`, relative to a neutral
one: `S / (1 - e^{-S})`, 1 at S = 0, 0 below S = -500, `S` above S = 500.
"""
@inline function fixation_probability(S::Float64)
    if abs(S) < 1e-10
        return 1.0
    end
    if S < -500.0
        return 0.0
    end
    if S > 500.0
        return S
    end
    exp_neg = exp(-S)
    return S / (1.0 - exp_neg)
end

"""
    RombergIntegrator()

Romberg quadrature over `[QUASI_ZERO, 1 - QUASI_ZERO]` with `EPS_R = 1e-6`: GRAPES'
own rule (`Grapes.cpp` `qromb`), and the integrator the GRAPES parity tests use.
The truncated interval drops mass from the rare bins: against the exact closed form
bin 1 is off by 0.19% at n = 20, 0.99% at n = 100 and 6.4% at n = 661, at every S.
It also overflows below S = -709, so grids reaching below S = -1500 are refused.
Use it to reproduce GRAPES, not for new analyses.
"""
struct RombergIntegrator <: AbstractIntegrator end

integrate(::RombergIntegrator, f, a::Float64, b::Float64) = _qromb(f, a, b)

"""
    GaussLegendreSojourn(; panel = 0.5, tail = 45.0)

Composite 10-point Gauss-Legendre under the logit substitution `x = 1/(1+e^{-t})`,
which cancels the `x(1-x)` denominator and removes the `QUASI_ZERO` truncation.
Panels follow the sampling kernel analytically from `(n, j, S)`. Max relative error
1.8e-11 against a BigFloat oracle over n <= 1322. The general rule for any `f`, and
the dN-correction producer.
"""
struct GaussLegendreSojourn <: AbstractIntegrator
    panel::Float64
    tail::Float64
end
GaussLegendreSojourn(; panel::Real = 0.5, tail::Real = 45.0) =
    GaussLegendreSojourn(Float64(panel), Float64(tail))

"""
    ExactSojourn(fallback = GaussLegendreSojourn(); max_series = 500.0)

The expected-SFS integral in closed form, the default:
`counts(S, j, n) = n / (j (n - j)) * (1 - M(n - j, n, -S)) / (1 - e^{-S})` with `M`
Kummer's 1F1, summed as an all-positive series with a running scale (no `e^{|S|}`
is formed). Agrees with a BigFloat quadrature of the original integrand to 1.9e-16.
Above `|S| = max_series` the series is slower than quadrature, so `fallback` takes
over; the dN correction and `integrate(::ExactSojourn, f, a, b)` always use it.
"""
struct ExactSojourn{F<:AbstractIntegrator} <: AbstractIntegrator
    fallback::F
    max_series::Float64
end
ExactSojourn(fallback::AbstractIntegrator = GaussLegendreSojourn(); max_series::Real = LARGE_S_VALUE) =
    ExactSojourn(fallback, Float64(max_series))

integrate(e::ExactSojourn, f, a::Float64, b::Float64) = integrate(e.fallback, f, a, b)

"""
    expected_count(integ, S, j, n, logc) -> Float64

Expected SFS bin `j` of `n` at selection coefficient `S` and theta = 1, with
`logc = log C(n, j)`. Optional hook; the fallback is GRAPES' expression through
`integrate`.
"""
function expected_count(integ::AbstractIntegrator, S::Float64, j::Int, n::Int,
                        logc::Float64)::Float64
    f = x -> _sojourn_density(x, S) * exp(logc + j * log(x) + (n - j) * log(1.0 - x))
    return integrate(integ, f, QUASI_ZERO, 1.0 - QUASI_ZERO) / 2.0
end

"""
    dn_correction_term(integ, S, n) -> Float64

GRAPES' correction for fixed differences that are still segregating, at `S` and
sample size `n`. Optional hook; the fallback is GRAPES' expression.
"""
function dn_correction_term(integ::AbstractIntegrator, S::Float64, n::Int)::Float64
    f = z -> (1.0 - exp(-S * (1.0 - z))) / (z * (1.0 - z) * (1.0 - exp(-S))) * z^n
    return integrate(integ, f, QUASI_ZERO, 1.0 - QUASI_ZERO)
end

# --- logit-substitution primitives -------------------------------------------
@inline _log1pexp(z::Float64) = z > 0 ? z + log1p(exp(-z)) : log1p(exp(z))
@inline _logistic(t::Float64) = t >= 0 ? 1.0 / (1.0 + exp(-t)) :
                                         (e = exp(t); e / (1.0 + e))
@inline _logexpm1(z::Float64) = z + log(-expm1(-z))

# w(x) = (1 - e^{-S(1-x)}) / (1 - e^{-S}) at x = logistic(t), with every exponential
# argument non-positive, so nothing overflows or cancels at |S| = 10000.
@inline function _sojourn_weight(t::Float64, S::Float64)
    if S > 0
        return (-expm1(-S * _logistic(-t))) / (-expm1(-S))
    else
        return exp(S * _logistic(t)) * (-expm1(S * _logistic(-t))) / (-expm1(S))
    end
end

# Window in t for the kernel C(n,j) x^{j-1}(1-x)^{n-j-1}; for S < 0 extended to
# t ~ -log|S|, where e^{S x} moves the mass.
@inline function _logit_window(j::Int, n::Int, S::Float64, tail::Float64)
    tstar = log(j / (n - j))
    sd    = sqrt(n / (j * (n - j)))
    tpk   = S < -1.0 ? min(tstar, -log(-S)) : tstar
    lo    = tpk   - max(8.0 * sd, tail / j)
    hi    = tstar + max(8.0 * sd, tail / (n - j))
    return lo, hi, sd
end

const _GL_MAX_PANELS = 8192

# 10-point Gauss-Legendre nodes and weights on [-1, 1].
const _GL10_X = (-0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
                 -0.4333953941292472, -0.1488743389816312,
                  0.1488743389816312,  0.4333953941292472,  0.6794095682990244,
                  0.8650633666889845,  0.9739065285171717)
const _GL10_W = (0.0666713443086881, 0.1494513491505806, 0.2190863625159820,
                 0.2692667193099963, 0.2955242247147529, 0.2955242247147529,
                 0.2692667193099963, 0.2190863625159820, 0.1494513491505806,
                 0.0666713443086881)

@inline function _gl_composite(f, lo::Float64, hi::Float64, npan::Int)
    h = (hi - lo) / npan
    hw = 0.5 * h
    acc = 0.0
    @inbounds for p in 0:(npan-1)
        mid = lo + (p + 0.5) * h
        s = 0.0
        for m in 1:10
            s += _GL10_W[m] * f(mid + hw * _GL10_X[m])
        end
        acc += hw * s
    end
    return acc
end

# Generic rule for any integrand: panel doubling until the estimate stops moving.
function integrate(g::GaussLegendreSojourn, f, a::Float64, b::Float64)
    prev = _gl_composite(f, a, b, 8)
    npan = 16
    while npan <= _GL_MAX_PANELS
        cur = _gl_composite(f, a, b, npan)
        abs(cur - prev) <= 1e-12 * abs(cur) && return cur
        prev = cur
        npan *= 2
    end
    return prev
end

function expected_count(g::GaussLegendreSojourn, S::Float64, j::Int, n::Int,
                        logc::Float64)::Float64
    lo, hi, sd = _logit_window(j, n, S, g.tail)
    step = min(sd, g.panel)
    npan = clamp(ceil(Int, (hi - lo) / step), 1, _GL_MAX_PANELS)
    # dx = x(1-x) dt cancels psi's denominator, and psi's 2 cancels the 1/2.
    return _gl_composite(lo, hi, npan) do t
        _sojourn_weight(t, S) *
            exp(logc - j * _log1pexp(-t) - (n - j) * _log1pexp(t))
    end
end

function dn_correction_term(g::GaussLegendreSojourn, S::Float64, n::Int)::Float64
    aS = abs(S)
    lo = min(0.0, aS > 1.0 ? -log(aS) : 0.0) - max(8.0, g.tail / n)
    hi = g.tail + 8.0
    step = min(g.panel, 0.25)
    npan = clamp(ceil(Int, (hi - lo) / step), 1, _GL_MAX_PANELS)
    return _gl_composite(lo, hi, npan) do t
        _sojourn_weight(t, S) * exp(-n * _log1pexp(-t))
    end
end

# The closed form, as all-positive series (Kummer's transform for S > 0), with a
# running scale exponent; terms peak at k ~ |S|.
function _exact_expected_count(S::Float64, j::Int, n::Int)::Float64
    pre = n / (j * (n - j))
    z = abs(S)
    total = 0.0
    logscale = 0.0
    k = 0
    if S > 0
        c = 1.0      # S^k / k!
        r = 1.0      # (j)_k / (n)_k
        while true
            k += 1
            c *= S / k
            r *= (j + k - 1) / (n + k - 1)
            term = c * (1.0 - r)
            total += term
            if total > 1e250
                total *= 1e-250; c *= 1e-250; term *= 1e-250
                logscale += 250 * log(10.0)
            end
            (k > z && term < total * 1e-17) && break
            k > 4_000_000 && error("ExactSojourn: series did not converge (S=$S, j=$j, n=$n)")
        end
    else
        t = 1.0      # (n-j)_k z^k / ((n)_k k!)
        while true
            k += 1
            t *= (n - j + k - 1) * z / ((n + k - 1) * k)
            total += t
            if total > 1e250
                total *= 1e-250; t *= 1e-250
                logscale += 250 * log(10.0)
            end
            (k > z && t < total * 1e-17) && break
            k > 4_000_000 && error("ExactSojourn: series did not converge (S=$S, j=$j, n=$n)")
        end
    end
    total <= 0.0 && return 0.0
    return pre * exp(log(total) + logscale - _logexpm1(z))
end

expected_count(e::ExactSojourn, S::Float64, j::Int, n::Int, logc::Float64)::Float64 =
    abs(S) <= e.max_series ? _exact_expected_count(S, j, n) :
                             expected_count(e.fallback, S, j, n, logc)

# No stable closed form for the dN term (its series alternates for S > 0).
dn_correction_term(e::ExactSojourn, S::Float64, n::Int)::Float64 =
    dn_correction_term(e.fallback, S, n)

# -----------------------------------------------------------------------------
# Numerical tail masses of a density on (0, inf), for models without a closed form
# -----------------------------------------------------------------------------

# Simpson on [0, Y], n panels (n even).
function _simpson(f, Y::Float64, n::Int)
    h = Y / n
    acc = f(0.0) + f(Y)
    @inbounds for i in 1:(n - 1)
        acc += (iseven(i) ? 2.0 : 4.0) * f(i * h)
    end
    return acc * h / 3
end

# Integrate `g` on [0, Ymax] in slabs of fixed width, refining each until it stops
# moving, and stopping once a slab no longer contributes.
function _log_quad(g, Ymax::Float64; rtol::Float64 = 1e-12, width::Float64 = 8.0)
    total = 0.0
    y0 = 0.0
    while y0 < Ymax
        w = min(width, Ymax - y0)
        slab = _simpson(y -> g(y0 + y), w, 16)
        for p in (32, 64, 128, 256, 512, 1024)
            s = _simpson(y -> g(y0 + y), w, p)
            converged = abs(s - slab) <= rtol * max(abs(s), 1e-300)
            slab = s
            converged && break
        end
        total += slab
        y0 += w
        abs(slab) <= rtol * max(abs(total), eps()) && break
    end
    return total
end

# Representable range of `s`: past these `s^(shape-1)` overflows.
const _S_FLOOR = 1e-300
const _S_CEIL = 1e300

# Mass of `density` on (0, b]; handles an integrable singularity at 0.
function _mass_below(density, b::Float64)
    b <= 0 && return 0.0
    Ymax = log(b / _S_FLOOR)
    return _log_quad(Ymax) do y
        s = b * exp(-y)
        v = density(s) * s
        isfinite(v) ? v : 0.0
    end
end

# Mass of `density` on [a, inf).
function _mass_above(density, a::Float64)
    a <= 0 && return 0.0
    Ymax = log(_S_CEIL / a)
    return _log_quad(Ymax) do y
        s = a * exp(y)
        v = density(s) * s
        isfinite(v) ? v : 0.0
    end
end

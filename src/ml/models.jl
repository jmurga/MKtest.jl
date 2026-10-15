# =============================================================================
# DFE models: traits, the generic `distribution` path, the built-in shapes,
# tail masses and CDFs.
#
# The built-in shapes keep their hand-written GRAPES densities, so their numbers
# do not move. A user model may instead define only `dfe_parameters` and
# `distribution`; every method below then has a generic fallback.
# =============================================================================

# -----------------------------------------------------------------------------
# Traits and their defaults
# -----------------------------------------------------------------------------

"""
    dfe_log_scale(m::AbstractDFEModel) -> Tuple of Bool

Which parameters the optimizer moves on a log scale, in `dfe_parameters` order.
Default: those whose lower bound is positive.
"""
dfe_log_scale(m::AbstractDFEModel) = map(b -> first(b) > 0, values(dfe_parameters(m)))

"""
    dfe_start_points(m::AbstractDFEModel) -> Vector{Vector{Float64}}

Starting points of the optimizer, in natural space and `dfe_parameters` order. The
engine always adds the centre of the box in optimizer space, and `nb_rand_start`
random points. Default: the first 8 points of a Halton sequence (bases 2, 3, 5, ...)
over the box in optimizer space, which fill it rather than sit on its diagonal and are
the same at every call. Override with starts you trust when the box is wide.
"""
function dfe_start_points(m::AbstractDFEModel)
    b = collect(values(dfe_parameters(m)))
    ls = collect(dfe_log_scale(m))
    lo = [ls[i] ? log(Float64(b[i][1])) : Float64(b[i][1]) for i in eachindex(b)]
    hi = [ls[i] ? log(Float64(b[i][2])) : Float64(b[i][2]) for i in eachindex(b)]
    return [[_from_opt(lo[i] + _halton(k, _PRIMES[i]) * (hi[i] - lo[i]), ls[i]) for i in eachindex(b)]
            for k in 1:8]
end

const _PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53)

# k-th element (k >= 1) of the van der Corput sequence in base `p`.
function _halton(k::Int, p::Int)
    f, r = 1.0, 0.0
    while k > 0
        f /= p
        r += f * (k % p)
        k ÷= p
    end
    return r
end

"""
    dfe_name(m::AbstractDFEModel) -> String

Name of the model in output tables. Default: the type name.
"""
dfe_name(m::AbstractDFEModel) = string(nameof(typeof(m)))

"""
    dfe_param_names(m) -> Vector{Symbol}

Parameter names of `m`, in order (the keys of `dfe_parameters`).
"""
dfe_param_names(m::AbstractDFEModel) = collect(keys(dfe_parameters(m)))

_nparams(m::AbstractDFEModel) = length(dfe_parameters(m))

"""
    discrete_components(m, p) -> iterable of (S, weight)

Point masses of the DFE, carried exactly rather than by a grid cell. Empty by
default. The weights carry their own mixture proportion.
"""
discrete_components(::AbstractDFEModel, p) = ()

# -----------------------------------------------------------------------------
# The generic `distribution` path
# -----------------------------------------------------------------------------

# Parameters as a tuple of the model's length: `_nparams` folds to a constant for a
# model whose `dfe_parameters` is a literal, so the call is type stable.
@inline _params_tuple(m::AbstractDFEModel, p::AbstractVector{Float64}) =
    ntuple(i -> @inbounds(p[i]), Val(_nparams(m)))

@inline _dist(m::AbstractDFEModel, p::AbstractVector{Float64}) =
    distribution(m, _params_tuple(m, p)...)

_has_distribution(m::AbstractDFEModel) =
    applicable(distribution, m, ntuple(_ -> 1.0, _nparams(m))...)

# Untyped beyond `m`, so a user method such as `dfe_density!(out, m::MyWrapper, p, s)`
# is strictly more specific and never ambiguous with this fallback.
function dfe_density!(out, m::AbstractDFEModel, p, s_grid)
    d = _dist(m, p)
    @inbounds for i in eachindex(out, s_grid)
        s = s_grid[i]
        out[i] = s < 0 ? pdf(d, -s) : 0.0
    end
    return out
end

# -----------------------------------------------------------------------------
# Density helpers of the built-in shapes (GRAPES parameterisation: shape, MEAN)
# -----------------------------------------------------------------------------

@inline function _gamma_pdf_neg(s::Float64, shape_m1::Float64, inv_scale::Float64, log_norm::Float64)
    x = -s
    return exp(shape_m1 * log(x) - x * inv_scale - log_norm)
end

@inline function _gamma_pdf_pos(s::Float64, shape_m1::Float64, inv_scale::Float64, log_norm::Float64)
    return exp(shape_m1 * log(s) - s * inv_scale - log_norm)
end

@inline function _beta_pdf_unit(x::Float64, a_m1::Float64, b_m1::Float64, log_norm::Float64)
    return exp(a_m1 * log(x) + b_m1 * log1p(-x) - log_norm)
end

@inline function _gamma_consts(shape::Float64, mean::Float64)
    scale     = mean / shape
    shape_m1  = shape - 1.0
    inv_scale = 1.0 / scale
    log_norm  = shape * log(scale) + loggamma(shape)
    return shape_m1, inv_scale, log_norm
end

const _FVec = AbstractVector{Float64}

function density_gamma_zero!(out::_FVec, params::_FVec, s_grid::_FVec)
    shape, mean = params[1], params[2]
    sm1, iscale, lnorm = _gamma_consts(shape, mean)
    @inbounds for i in eachindex(s_grid)
        s = s_grid[i]
        out[i] = s < 0 ? _gamma_pdf_neg(s, sm1, iscale, lnorm) : 0.0
    end
    return out
end

function density_gamma_expo!(out::_FVec, params::_FVec, s_grid::_FVec)
    shape, mean, prop_ben, mean_ben = params[1], params[2], params[3], params[4]
    sm1, iscale, lnorm = _gamma_consts(shape, mean)
    w_del = 1.0 - prop_ben
    inv_mben = 1.0 / mean_ben
    @inbounds for i in eachindex(s_grid)
        s = s_grid[i]
        if s < 0
            out[i] = w_del * _gamma_pdf_neg(s, sm1, iscale, lnorm)
        elseif s > 0
            out[i] = prop_ben * exp(-s * inv_mben) * inv_mben
        else
            out[i] = 0.0
        end
    end
    return out
end

function density_gamma_gamma!(out::_FVec, params::_FVec, s_grid::_FVec)
    shape_del, mean_del, prop_ben, shape_ben, mean_ben =
        params[1], params[2], params[3], params[4], params[5]
    sm1d, iscaled, lnormd = _gamma_consts(shape_del, mean_del)
    sm1b, iscaleb, lnormb = _gamma_consts(shape_ben, mean_ben)
    w_del = 1.0 - prop_ben
    @inbounds for i in eachindex(s_grid)
        s = s_grid[i]
        if s < 0
            out[i] = w_del * _gamma_pdf_neg(s, sm1d, iscaled, lnormd)
        elseif s > 0
            out[i] = prop_ben * _gamma_pdf_pos(s, sm1b, iscaleb, lnormb)
        else
            out[i] = 0.0
        end
    end
    return out
end

function density_expo_zero!(out::_FVec, params::_FVec, s_grid::_FVec)
    mean = params[1]
    inv_mean = 1.0 / mean
    @inbounds for i in eachindex(s_grid)
        s = s_grid[i]
        out[i] = s < 0 ? exp(s * inv_mean) * inv_mean : 0.0
    end
    return out
end

function density_expo_expo!(out::_FVec, params::_FVec, s_grid::_FVec)
    mean_del, prop_ben, mean_ben = params[1], params[2], params[3]
    inv_mdel = 1.0 / mean_del
    inv_mben = 1.0 / mean_ben
    w_del = 1.0 - prop_ben
    @inbounds for i in eachindex(s_grid)
        s = s_grid[i]
        if s < 0
            out[i] = w_del * exp(s * inv_mdel) * inv_mdel
        elseif s > 0
            out[i] = prop_ben * exp(-s * inv_mben) * inv_mben
        else
            out[i] = 0.0
        end
    end
    return out
end

function density_expo_gamma!(out::_FVec, params::_FVec, s_grid::_FVec)
    mean_del, prop_ben, shape_ben, mean_ben = params[1], params[2], params[3], params[4]
    inv_mdel = 1.0 / mean_del
    sm1b, iscaleb, lnormb = _gamma_consts(shape_ben, mean_ben)
    w_del = 1.0 - prop_ben
    @inbounds for i in eachindex(s_grid)
        s = s_grid[i]
        if s < 0
            out[i] = w_del * exp(s * inv_mdel) * inv_mdel
        elseif s > 0
            out[i] = prop_ben * _gamma_pdf_pos(s, sm1b, iscaleb, lnormb)
        else
            out[i] = 0.0
        end
    end
    return out
end

function density_displaced_gamma!(out::_FVec, params::_FVec, s_grid::_FVec)
    shape, mean, displacement = params[1], params[2], params[3]
    sm1, iscale, lnorm = _gamma_consts(shape, mean)
    @inbounds for i in eachindex(s_grid)
        s = s_grid[i]
        x = -s - displacement   # > 0 only when s < -displacement
        out[i] = x > 0 ? exp(sm1 * log(x) - x * iscale - lnorm) : 0.0
    end
    return out
end

function density_scaled_beta!(out::_FVec, params::_FVec, s_grid::_FVec)
    Ba, Bb, med_prop = params[1], params[2], params[3]
    boundary  = SCALED_BETA_BOUNDARY
    scale     = 2.0 * boundary
    inv_scale = 1.0 / scale
    a_m1      = Ba - 1.0
    b_m1      = Bb - 1.0
    log_norm  = loggamma(Ba) + loggamma(Bb) - loggamma(Ba + Bb)
    factor    = med_prop * inv_scale
    @inbounds for i in eachindex(s_grid)
        s = s_grid[i]
        if s < -boundary || s > boundary
            out[i] = 0.0
        else
            x = s * inv_scale + 0.5
            out[i] = factor * _beta_pdf_unit(x, a_m1, b_m1, log_norm)
        end
    end
    return out
end

# |S| ~ LogNormal(mu, sigma) with mu = log(mean) - sigma^2/2.
function density_lognormal_zero!(out::_FVec, params::_FVec, s_grid::_FVec)
    mean_neg, sigma = params[1], params[2]
    mu                = log(mean_neg) - 0.5 * sigma^2
    log_sigma_sqrt2pi = log(sigma) + 0.5 * log(2.0 * π)
    half_inv_sigma2   = 0.5 / (sigma * sigma)
    @inbounds for i in eachindex(s_grid)
        s = s_grid[i]
        if s < 0.0
            lx = log(-s)
            out[i] = exp(-(lx - mu)^2 * half_inv_sigma2 - lx - log_sigma_sqrt2pi)
        else
            out[i] = 0.0
        end
    end
    return out
end

# -----------------------------------------------------------------------------
# Built-in shapes
# -----------------------------------------------------------------------------

"""
    GammaZero()

Gamma DFE of deleterious mutations (`del_shape`, `del_mean`), no beneficial class.
GRAPES' GammaZero.
"""
struct GammaZero <: AbstractDFEModel end

"""
    GammaExpo()

Gamma deleterious DFE plus an exponential beneficial class of weight `ben_prop` and
mean `ben_mean`. GRAPES' GammaExpo.
"""
struct GammaExpo <: AbstractDFEModel end

"""
    GammaGamma()

Gamma deleterious DFE plus a gamma beneficial class (`ben_prop`, `ben_shape`,
`ben_mean`). GRAPES' GammaGamma.
"""
struct GammaGamma <: AbstractDFEModel end

"""
    ExpoZero()

Exponential deleterious DFE (`del_mean`). GRAPES' ExpoZero.
"""
struct ExpoZero <: AbstractDFEModel end

"""
    ExpoExpo()

Exponential deleterious DFE plus an exponential beneficial class. GRAPES' ExpoExpo.
"""
struct ExpoExpo <: AbstractDFEModel end

"""
    ExpoGamma()

Exponential deleterious DFE plus a gamma beneficial class. GRAPES' ExpoGamma.
"""
struct ExpoGamma <: AbstractDFEModel end

"""
    DisplGamma()

Displaced gamma (`del_shape`, `del_mean`, displacement `del_shift <= 0`): the gamma variable
is `-S - s0`, so a little mass sits at positive S. GRAPES' DisplGamma.
"""
struct DisplGamma <: AbstractDFEModel end

"""
    ScaledBeta()

A beta density on S in [-25, 25] (`beta_a`, `beta_b`) scaled by `beta_mass`. GRAPES' ScaledBeta.
"""
struct ScaledBeta <: AbstractDFEModel end

"""
    LogNormalZero()

Lognormal |S| of deleterious mutations (`del_mean`, `del_sigma` = the sd of log |S|).
"""
struct LogNormalZero <: AbstractDFEModel end

# Parameter names follow `<component>_<quantity>` (component del / ben / neu / lethal, or a
# family prefix; quantity prop / mean / shape / sigma / shift); means are of |S|, S = 4 Ne s.
# `ben_prop`'s upper bound 0.1 is GRAPES' MAX_ORIENTATION_ERROR reused as a mixture
# weight bound. It is kept because GRAPES parity depends on it; widen per fit with
# `bounds = (ben_prop = (1e-6, 0.99),)`.
dfe_parameters(::GammaZero)  = (del_shape = (0.05, 100.0), del_mean = (1.0, 1e7))
dfe_parameters(::GammaExpo)  = (del_shape = (0.05, 100.0), del_mean = (1.0, 1e7),
                                ben_prop = (1e-6, 0.1), ben_mean = (1e-5, 1e4))
dfe_parameters(::GammaGamma) = (del_shape = (0.05, 100.0), del_mean = (1.0, 1e7),
                                ben_prop = (1e-6, 0.1), ben_shape = (0.05, 100.0),
                                ben_mean = (1e-5, 1e4))
dfe_parameters(::ExpoZero)   = (del_mean = (1.0, 1e7),)
dfe_parameters(::ExpoExpo)   = (del_mean = (1.0, 1e7), ben_prop = (1e-6, 0.1),
                                ben_mean = (1e-5, 1e4))
dfe_parameters(::ExpoGamma)  = (del_mean = (1.0, 1e7), ben_prop = (1e-6, 0.1),
                                ben_shape = (0.05, 100.0), ben_mean = (1e-5, 1e4))
dfe_parameters(::DisplGamma) = (del_shape = (0.05, 100.0), del_mean = (1.0, 1e7),
                                del_shift = (-100.0, 0.0))
dfe_parameters(::ScaledBeta) = (beta_a = (QUASI_ZERO, 10000.0), beta_b = (QUASI_ZERO, 10000.0),
                                beta_mass = (0.0, 1.0))
dfe_parameters(::LogNormalZero) = (del_mean = (1.0, 1e7), del_sigma = (0.1, 10.0))

# GRAPES' scales (shape and proportions linear, means log), not the lower > 0 rule.
dfe_log_scale(::GammaZero)  = (false, true)
dfe_log_scale(::GammaExpo)  = (false, true, false, true)
dfe_log_scale(::GammaGamma) = (false, true, false, true, true)
dfe_log_scale(::ExpoZero)   = (true,)
dfe_log_scale(::ExpoExpo)   = (true, false, true)
dfe_log_scale(::ExpoGamma)  = (true, false, false, true)
dfe_log_scale(::DisplGamma) = (false, true, false)
dfe_log_scale(::ScaledBeta) = (false, false, false)
dfe_log_scale(::LogNormalZero) = (true, false)

# GRAPES' effective starting points (all GammaExpo starts use del_mean = 10000, as
# GRAPES does through a typo in its source).
dfe_start_points(::GammaZero) = [[1.0, 1000.0], [0.5, 100000.0], [0.2, 10000000.0], [0.3, 20000.0]]
dfe_start_points(::GammaExpo) = [
    [1.0, 10000.0, 0.025, 1.0], [0.5, 10000.0, 0.025, 10.0],
    [0.2, 10000.0, 0.001, 100.0], [0.2, 10000.0, 0.025, 1000.0]]
dfe_start_points(::GammaGamma) = [
    [1.0, 10000.0, 0.025, 5.0, 1.0], [0.5, 10000.0, 0.025, 1.0, 10.0],
    [0.2, 10000.0, 0.025, 0.5, 100.0], [0.2, 10000.0, 0.025, 0.2, 1000.0]]
dfe_start_points(::ExpoZero) = [[100.0], [1000.0], [10000.0], [100000.0]]
dfe_start_points(::ExpoExpo) = [
    [100.0, 0.025, 1.0], [1000.0, 0.025, 10.0],
    [10000.0, 0.025, 100.0], [100000.0, 0.025, 1000.0]]
dfe_start_points(::ExpoGamma) = [
    [100.0, 0.025, 5.0, 1.0], [1000.0, 0.025, 1.0, 10.0],
    [10000.0, 0.025, 0.5, 100.0], [100000.0, 0.025, 0.2, 1000.0]]
dfe_start_points(::DisplGamma) = [
    [1.0, 10000.0, -2.0], [1.0, 10000.0, -0.2], [0.5, 100000.0, -0.5], [0.2, 1e7, -0.1]]
dfe_start_points(::ScaledBeta) = [
    [0.01, 3.0, 0.5], [0.1, 1.0, 0.5], [0.1, 5.0, 0.5], [1.0, 10.0, 0.5]]
dfe_start_points(::LogNormalZero) = [
    [1000.0, 1.0], [10000.0, 0.5], [100.0, 2.0], [100000.0, 1.5]]

dfe_density!(out::_FVec, ::GammaZero, p::_FVec, s::_FVec)  = density_gamma_zero!(out, p, s)
dfe_density!(out::_FVec, ::GammaExpo, p::_FVec, s::_FVec)  = density_gamma_expo!(out, p, s)
dfe_density!(out::_FVec, ::GammaGamma, p::_FVec, s::_FVec) = density_gamma_gamma!(out, p, s)
dfe_density!(out::_FVec, ::ExpoZero, p::_FVec, s::_FVec)   = density_expo_zero!(out, p, s)
dfe_density!(out::_FVec, ::ExpoExpo, p::_FVec, s::_FVec)   = density_expo_expo!(out, p, s)
dfe_density!(out::_FVec, ::ExpoGamma, p::_FVec, s::_FVec)  = density_expo_gamma!(out, p, s)
dfe_density!(out::_FVec, ::DisplGamma, p::_FVec, s::_FVec) = density_displaced_gamma!(out, p, s)
dfe_density!(out::_FVec, ::ScaledBeta, p::_FVec, s::_FVec) = density_scaled_beta!(out, p, s)
dfe_density!(out::_FVec, ::LogNormalZero, p::_FVec, s::_FVec) = density_lognormal_zero!(out, p, s)

"""
    NeuGamma()

A point mass `neu_prop` at S = 0 plus `1 - neu_prop` on a gamma deleterious DFE (fitdadi's
`neugamma`). The neutral class is a discrete component, not a grid cell.
"""
struct NeuGamma <: AbstractDFEModel end

dfe_parameters(::NeuGamma) = (neu_prop = (0.0, 1.0), del_shape = (0.05, 100.0), del_mean = (1.0, 1e7))
dfe_log_scale(::NeuGamma) = (false, false, true)
dfe_start_points(::NeuGamma) = [
    [0.1, 1.0, 1000.0], [0.25, 0.5, 100000.0], [0.5, 0.2, 10000000.0], [0.05, 0.3, 20000.0]]

function density_neugamma!(out::_FVec, params::_FVec, s_grid::_FVec)
    p_neu, shape, mean = params[1], params[2], params[3]
    sm1, iscale, lnorm = _gamma_consts(shape, mean)
    w = 1.0 - p_neu
    @inbounds for i in eachindex(s_grid)
        sv = s_grid[i]
        out[i] = sv < 0 ? w * _gamma_pdf_neg(sv, sm1, iscale, lnorm) : 0.0
    end
    return out
end

dfe_density!(out::_FVec, ::NeuGamma, p::_FVec, s::_FVec) = density_neugamma!(out, p, s)
discrete_components(::NeuGamma, p::AbstractVector{Float64}) = ((0.0, p[1]),)

const _LOG_SQRT_PI = log(sqrt(pi))

"""
    FGMBesselK()

The DFE implied by Fisher's geometric model (Lourenco, Galtier & Glemin 2011),
GRAPES' FGMBesselK, with GRAPES' bounds, starts and its `C > 50 -> 0` cutoff. No
closed-form CDF; masses use numerical quadrature, as GRAPES does.
"""
struct FGMBesselK <: AbstractDFEModel end

dfe_parameters(::FGMBesselK) = (fgm_m = (1.0, 10.0), fgm_nsz2 = (1.0, 10000.0),
                                fgm_sigma = (0.05, 100.0), fgm_scale = (QUASI_ZERO, 1e4))
dfe_log_scale(::FGMBesselK) = (false, false, false, true)
dfe_start_points(::FGMBesselK) = [
    [2.0, 100.0, 0.5, 1000.0], [2.0, 10.0, 0.3, 100.0],
    [4.0, 1000.0, 1.0, 1000.0], [1.5, 100.0, 0.54, 5000.0]]

# Scalar density, so S = 0 can be the mean of the two one-sided limits (GRAPES).
function _besselk_density(x0::Float64, m::Float64, nsz2::Float64, sig::Float64, scale::Float64)
    if x0 == 0.0
        return (_besselk_density(-QUASI_ZERO, m, nsz2, sig, scale) +
                _besselk_density(QUASI_ZERO, m, nsz2, sig, scale)) / 2.0
    end
    x = x0 / scale
    s2 = sig * sig
    x2 = x * x
    lp = -log(2.0) * (m + 1.0) / 2.0
    lp -= x * nsz2 / 4.0
    lp += log(nsz2) / 2.0
    A = (1.0 + 4.0 / (nsz2 * s2)) / x2
    lp += log(A) * (1.0 - m) / 4.0
    lp -= m * log(sig)
    lp -= _LOG_SQRT_PI
    lp -= loggamma(m / 2.0)
    B = 1.0 / (nsz2 * x2 * (nsz2 + 4.0 / s2))
    C = 1.0 / (4.0 * sqrt(B))
    C > 50.0 && return 0.0
    lp += log(besselk((m - 1.0) / 2.0, C))
    v = exp(lp) / scale
    return isfinite(v) ? v : 0.0
end

function density_fgm_besselk!(out::_FVec, params::_FVec, s_grid::_FVec)
    m, nsz2, sig, scale = params[1], params[2], params[3], params[4]
    @inbounds for i in eachindex(s_grid)
        out[i] = _besselk_density(s_grid[i], m, nsz2, sig, scale)
    end
    return out
end

dfe_density!(out::_FVec, ::FGMBesselK, p::_FVec, s::_FVec) = density_fgm_besselk!(out, p, s)

"""
    FixedSpikes()
    FixedSpikes(positions)

All mass on a fixed ladder of S values with free weights (Kousathanas & Keightley
2013); the default ladder is `(0, -1, -10, -100, -1000, -10000)`. The K weights are
parameterised by K-1 stick-breaking fractions `spike_w1 ... spike_w{K-1}` in (0, 1); use
`spike_weights(m, p)` for the weights.
"""
struct FixedSpikes{K} <: AbstractDFEModel
    positions::NTuple{K,Float64}
end

const KK_SPIKE_LADDER = (0.0, -1.0, -10.0, -100.0, -1000.0, -10000.0)
FixedSpikes() = FixedSpikes(KK_SPIKE_LADDER)
FixedSpikes(pos::AbstractVector{<:Real}) = FixedSpikes(Tuple(Float64.(pos)))

dfe_name(::FixedSpikes) = "FixedSpikes"
dfe_parameters(::FixedSpikes{K}) where {K} =
    NamedTuple{ntuple(i -> Symbol("spike_w", i), K - 1)}(ntuple(_ -> (1e-6, 1.0 - 1e-6), K - 1))
dfe_log_scale(::FixedSpikes{K}) where {K} = ntuple(_ -> false, K - 1)

# Mass on the neutral end, on the deleterious end, and the flat simplex point.
function dfe_start_points(::FixedSpikes{K}) where {K}
    flat  = [1.0 / (K - i + 1) for i in 1:(K - 1)]
    front = fill(0.5, K - 1)
    back  = fill(0.05, K - 1)
    return [flat, front, back]
end

"""
    spike_weights(m::FixedSpikes, p) -> Vector{Float64}

The K weights implied by the K-1 stick-breaking fractions. Sums to 1.
"""
function spike_weights(::FixedSpikes{K}, params::AbstractVector{Float64}) where {K}
    w = Vector{Float64}(undef, K)
    rest = 1.0
    @inbounds for i in 1:(K - 1)
        w[i] = rest * params[i]
        rest -= w[i]
    end
    w[K] = rest
    return w
end

# All the mass is discrete.
dfe_density!(out::_FVec, ::FixedSpikes, ::_FVec, ::_FVec) = fill!(out, 0.0)
discrete_components(m::FixedSpikes, p::AbstractVector{Float64}) = zip(m.positions, spike_weights(m, p))

"""
    WithLethal(base)

Adds a lethal class of weight `lethal_prop` to any DFE: mutations that never segregate
and never fix (GRAPES' `-p_lethal`, Galtier & Rousselle 2020). Appends the parameter
`lethal_prop` in [0, 0.99] and scales the base density, point masses and tail masses by
`1 - lethal_prop`, which is where GRAPES applies it. Under `SaturatedDivergence` it is
identified by the level of Pn against Ps only.
"""
struct WithLethal{M<:AbstractDFEModel} <: AbstractDFEModel
    base::M
end

dfe_name(m::WithLethal) = "WithLethal{" * dfe_name(m.base) * "}"
function dfe_parameters(m::WithLethal)
    base = dfe_parameters(m.base)
    # NamedTuple merging would silently drop a second `lethal_prop`.
    haskey(base, :lethal_prop) && throw(ArgumentError("WithLethal: the base model already has a lethal_prop parameter"))
    return (base..., lethal_prop = (0.0, 0.99))
end
dfe_log_scale(m::WithLethal) = (dfe_log_scale(m.base)..., false)
dfe_start_points(m::WithLethal) = [vcat(s, 0.1) for s in dfe_start_points(m.base)]

# A view, not a copy: this runs at every likelihood evaluation.
@inline _base_params(p::AbstractVector{Float64}) = @view p[1:(length(p) - 1)]

function dfe_density!(out::_FVec, m::WithLethal, p::_FVec, s_grid::_FVec)
    dfe_density!(out, m.base, _base_params(p), s_grid)
    w = 1.0 - p[end]
    @inbounds for i in eachindex(out)
        out[i] *= w
    end
    return out
end

discrete_components(m::WithLethal, p::AbstractVector{Float64}) =
    ((S, (1.0 - p[end]) * w) for (S, w) in discrete_components(m.base, _base_params(p)))

# -----------------------------------------------------------------------------
# Tail masses (deleterious side; unnormalised, they carry the mixture weight)
# -----------------------------------------------------------------------------

"""
    dfe_mass_below(m, p, b) -> Float64

Deleterious DFE mass with `|S| <= b` (point masses included). Closed form for the
gamma, exponential and lognormal families, `cdf(distribution(...), b)` for a model
with `distribution`, numerical quadrature of the density otherwise.
"""
function dfe_mass_below end

"""
    dfe_mass_above(m, p, a) -> Float64

Deleterious DFE mass with `|S| >= a`, by the same rules as `dfe_mass_below`.
"""
function dfe_mass_above end

# Scalar density at S = -x, for the numerical fallback. Two buffers, so a
# `dfe_density!` is never handed aliased input and output.
function _dfe_neg_scalar(m::AbstractDFEModel, params::AbstractVector{Float64})
    sbuf, obuf = zeros(1), zeros(1)
    return function (x::Float64)
        sbuf[1] = -x
        dfe_density!(obuf, m, params, sbuf)
        return obuf[1]
    end
end

function _dfe_pos_scalar(m::AbstractDFEModel, params::AbstractVector{Float64})
    sbuf, obuf = zeros(1), zeros(1)
    return function (x::Float64)
        sbuf[1] = x
        dfe_density!(obuf, m, params, sbuf)
        return obuf[1]
    end
end

function dfe_mass_below(m::AbstractDFEModel, params, b)
    _has_distribution(m) && return b <= 0 ? 0.0 : cdf(_dist(m, params), b)
    return _mass_below(_dfe_neg_scalar(m, params), b)
end
function dfe_mass_above(m::AbstractDFEModel, params, a)
    _has_distribution(m) && return ccdf(_dist(m, params), a)
    return _mass_above(_dfe_neg_scalar(m, params), a)
end

# Gamma family in (shape, MEAN), so scale = mean / shape.
_gamma_below(shape, mean, b) = gamma_inc(shape, b / (mean / shape))[1]
_gamma_above(shape, mean, a) = gamma_inc(shape, a / (mean / shape))[2]

dfe_mass_below(::GammaZero, p::_FVec, b::Float64) = _gamma_below(p[1], p[2], b)
dfe_mass_above(::GammaZero, p::_FVec, a::Float64) = _gamma_above(p[1], p[2], a)
dfe_mass_below(::GammaExpo, p::_FVec, b::Float64) = (1 - p[3]) * _gamma_below(p[1], p[2], b)
dfe_mass_above(::GammaExpo, p::_FVec, a::Float64) = (1 - p[3]) * _gamma_above(p[1], p[2], a)
dfe_mass_below(::GammaGamma, p::_FVec, b::Float64) = (1 - p[3]) * _gamma_below(p[1], p[2], b)
dfe_mass_above(::GammaGamma, p::_FVec, a::Float64) = (1 - p[3]) * _gamma_above(p[1], p[2], a)
dfe_mass_below(::ExpoZero, p::_FVec, b::Float64) = -expm1(-b / p[1])
dfe_mass_above(::ExpoZero, p::_FVec, a::Float64) = exp(-a / p[1])
dfe_mass_below(::ExpoExpo, p::_FVec, b::Float64) = (1 - p[2]) * -expm1(-b / p[1])
dfe_mass_above(::ExpoExpo, p::_FVec, a::Float64) = (1 - p[2]) * exp(-a / p[1])
dfe_mass_below(::ExpoGamma, p::_FVec, b::Float64) = (1 - p[2]) * -expm1(-b / p[1])
dfe_mass_above(::ExpoGamma, p::_FVec, a::Float64) = (1 - p[2]) * exp(-a / p[1])

function dfe_mass_below(::LogNormalZero, p::_FVec, b::Float64)
    b <= 0 && return 0.0
    mu = log(p[1]) - 0.5 * p[2]^2
    return 0.5 * erfc(-(log(b) - mu) / (p[2] * sqrt(2.0)))
end
dfe_mass_above(m::LogNormalZero, p::_FVec, a::Float64) = 1.0 - dfe_mass_below(m, p, a)

# The neutral class sits at |S| = 0, so it belongs to mass_below only.
dfe_mass_below(::NeuGamma, p::_FVec, b::Float64) = p[1] + (1 - p[1]) * _gamma_below(p[2], p[3], b)
dfe_mass_above(::NeuGamma, p::_FVec, a::Float64) = (1 - p[1]) * _gamma_above(p[2], p[3], a)

dfe_mass_below(m::FixedSpikes, p::_FVec, b::Float64) =
    sum(w for (S, w) in discrete_components(m, p) if S <= 0 && abs(S) <= b; init = 0.0)
dfe_mass_above(m::FixedSpikes, p::_FVec, a::Float64) =
    sum(w for (S, w) in discrete_components(m, p) if S <= 0 && abs(S) >= a; init = 0.0)

dfe_mass_below(m::WithLethal, p::_FVec, b::Float64) = (1.0 - p[end]) * dfe_mass_below(m.base, _base_params(p), b)
dfe_mass_above(m::WithLethal, p::_FVec, a::Float64) = (1.0 - p[end]) * dfe_mass_above(m.base, _base_params(p), a)

# -----------------------------------------------------------------------------
# CDF of the continuous part, for IntervalCDF (point masses excluded: they are
# carried exactly by `discrete_components`)
# -----------------------------------------------------------------------------

"""
    dfe_mass_ben_below(m, p, s) -> Float64

Beneficial mass with `0 < S <= s`. Zero by default (no beneficial class).
"""
dfe_mass_ben_below(::AbstractDFEModel, p, s) = 0.0
dfe_mass_ben_below(::GammaExpo, p::_FVec, s::Float64)  = p[3] * -expm1(-s / p[4])
dfe_mass_ben_below(::GammaGamma, p::_FVec, s::Float64) = p[3] * gamma_inc(p[4], s / (p[5] / p[4]))[1]
dfe_mass_ben_below(::ExpoExpo, p::_FVec, s::Float64)   = p[2] * -expm1(-s / p[3])
dfe_mass_ben_below(::ExpoGamma, p::_FVec, s::Float64)  = p[2] * gamma_inc(p[3], s / (p[4] / p[3]))[1]

_ben_below_numeric(m::AbstractDFEModel, p::AbstractVector{Float64}, s::Float64) =
    s <= 0 ? 0.0 : _mass_below(_dfe_pos_scalar(m, p), s)
dfe_mass_ben_below(m::DisplGamma, p::_FVec, s::Float64) = _ben_below_numeric(m, p, s)
dfe_mass_ben_below(m::ScaledBeta, p::_FVec, s::Float64) = _ben_below_numeric(m, p, s)
dfe_mass_ben_below(m::FGMBesselK, p::_FVec, s::Float64) = _ben_below_numeric(m, p, s)

# The numerical fallback's `_mass_above` is 0 at a = 0 by construction (its
# substitution s = a e^y), so a model on it integrates the half-line in two pieces.
_del_total_numeric(m, p) = dfe_mass_below(m, p, 1.0) + dfe_mass_above(m, p, 1.0)

const _CLOSED_FORM_DEL = Union{GammaZero,GammaExpo,GammaGamma,ExpoZero,ExpoExpo,ExpoGamma,
                               LogNormalZero,NeuGamma}

"""
    dfe_mass_del_total(m, p) -> Float64

Total continuous mass on the deleterious side.
"""
function dfe_mass_del_total(m::AbstractDFEModel, p)
    _has_distribution(m) && return ccdf(_dist(m, p), 0.0)
    return _del_total_numeric(m, p)
end
dfe_mass_del_total(m::_CLOSED_FORM_DEL, p::_FVec) = dfe_mass_above(m, p, 0.0)
dfe_mass_del_total(m::WithLethal, p::_FVec) = (1.0 - p[end]) * dfe_mass_del_total(m.base, _base_params(p))

"""
    dfe_cdf_continuous(m, p, S) -> Float64

`P(X <= S)` for the continuous part of the DFE on the signed S axis.
"""
function dfe_cdf_continuous(m::AbstractDFEModel, p, S)
    S < 0 && return dfe_mass_above(m, p, -S)
    return dfe_mass_del_total(m, p) + dfe_mass_ben_below(m, p, S)
end

# Closed forms where the tail masses are numerical.
function dfe_cdf_continuous(::ScaledBeta, p::_FVec, S::Float64)
    B = SCALED_BETA_BOUNDARY
    S <= -B && return 0.0
    S >= B && return p[3]
    return p[3] * beta_inc(p[1], p[2], S / (2B) + 0.5)[1]
end

# The gamma variable is g = -S - s0, so P(X <= S) = P(g >= -S - s0).
function dfe_cdf_continuous(::DisplGamma, p::_FVec, S::Float64)
    shape, mean, s0 = p[1], p[2], p[3]
    g = -S - s0
    g <= 0 && return 1.0
    return gamma_inc(shape, g / (mean / shape))[2]
end

dfe_cdf_continuous(m::WithLethal, p::_FVec, S::Float64) =
    (1.0 - p[end]) * dfe_cdf_continuous(m.base, _base_params(p), S)
dfe_cdf_continuous(::FixedSpikes, ::_FVec, ::Float64) = 0.0

"""
    dfe_cdf_continuous!(out, m, p, edges) -> out

`dfe_cdf_continuous` at every edge (ascending). The form the engine calls; it hoists
the total deleterious mass out of the loop.
"""
function dfe_cdf_continuous!(out, m::AbstractDFEModel, p, edges)
    if _has_distribution(m)
        d = _dist(m, p)
        tot = ccdf(d, 0.0)
        @inbounds for i in eachindex(edges)
            S = edges[i]
            out[i] = S < 0 ? ccdf(d, -S) : tot
        end
        return out
    end
    tot_del = dfe_mass_del_total(m, p)
    @inbounds for i in eachindex(edges)
        S = edges[i]
        out[i] = S < 0 ? dfe_mass_above(m, p, -S) : tot_del + dfe_mass_ben_below(m, p, S)
    end
    return out
end

_cdf_loop!(out, m, p, edges) =
    (@inbounds for i in eachindex(edges); out[i] = dfe_cdf_continuous(m, p, edges[i]); end; out)
dfe_cdf_continuous!(out::_FVec, m::ScaledBeta, p::_FVec, e::_FVec) = _cdf_loop!(out, m, p, e)
dfe_cdf_continuous!(out::_FVec, m::DisplGamma, p::_FVec, e::_FVec) = _cdf_loop!(out, m, p, e)
dfe_cdf_continuous!(out::_FVec, ::FixedSpikes, ::_FVec, ::_FVec) = fill!(out, 0.0)
dfe_cdf_continuous!(out::_FVec, m::WithLethal, p::_FVec, edges::_FVec) =
    (dfe_cdf_continuous!(out, m.base, _base_params(p), edges); out .*= (1.0 - p[end]); out)

# FGMBesselK: no analytic CDF and the adaptive tail quadrature is not accurate
# enough here, but IntervalCDF only uses differences, so integrate each interval
# directly with GL10 (the additive constant is irrelevant).
function dfe_cdf_continuous!(out::_FVec, m::FGMBesselK, p::_FVec, e::_FVec)
    buf_s, buf_o = zeros(10), zeros(10)
    out[1] = 0.0
    @inbounds for i in 1:(length(e) - 1)
        lo, hi = e[i], e[i+1]
        acc = 0.0
        if hi > lo
            c, h = (lo + hi) / 2, (hi - lo) / 2
            for q in 1:10
                buf_s[q] = c + h * _GL10_X[q]
            end
            dfe_density!(buf_o, m, p, buf_s)
            for q in 1:10
                acc += h * _GL10_W[q] * buf_o[q]
            end
        end
        out[i+1] = out[i] + max(acc, 0.0)
    end
    return out
end

# -----------------------------------------------------------------------------
# Rescaling to the ancestral Ne (divergence side)
# -----------------------------------------------------------------------------

# Parameters on the S scale (multiplied by r), and parameters with no S unit.
# `beta_mass` (GRAPES' med_prop) is the one exception: ScaledBeta's support is fixed at [-25, 25],
# so GRAPES divides its weight by r instead.
const _S_SCALE_NAMES = (:del_mean, :ben_mean, :del_shift, :fgm_scale)
const _UNITLESS_NAMES = (:del_shape, :ben_shape, :ben_prop, :neu_prop, :lethal_prop, :beta_a, :beta_b,
                         :del_sigma, :fgm_m, :fgm_nsz2, :fgm_sigma)

"""
    dfe_ne_multipliers(m::AbstractDFEModel, r) -> Vector{Float64}

Per-parameter factors that rescale the DFE from the current to the ancestral Ne
`r * Ne`: S-scale parameters times `r` (`del_mean`, `ben_mean`, `del_shift`,
`fgm_scale`), `beta_mass` divided by `r` (GRAPES' med_prop rule), unitless parameters unchanged. A
model with other parameter names must define this method to be used with an Ne
ratio other than 1.
"""
function dfe_ne_multipliers(m::AbstractDFEModel, r::Float64)
    nm = dfe_param_names(m)
    mult = ones(length(nm))
    r == 1.0 && return mult
    for (i, k) in enumerate(nm)
        if k in _S_SCALE_NAMES
            mult[i] = r
        elseif k === :beta_mass
            mult[i] = 1.0 / r
        elseif !(k in _UNITLESS_NAMES)
            throw(ArgumentError("$(dfe_name(m)): parameter $k has no known scaling with Ne; " *
                                "define MKtest.dfe_ne_multipliers for this model to use an Ne ratio != 1"))
        end
    end
    return mult
end

dfe_ne_multipliers(m::FixedSpikes, r::Float64) = r == 1.0 ? ones(_nparams(m)) :
    throw(ArgumentError("FixedSpikes: an Ne ratio != 1 moves the spike positions, which are not " *
                        "parameters; rescaling them is not implemented"))

# The built-in shapes have no `distribution`; saying so statically keeps the generic
# mass and CDF methods free of a runtime `applicable` lookup for them.
_has_distribution(::Union{GammaZero,GammaExpo,GammaGamma,ExpoZero,ExpoExpo,ExpoGamma,DisplGamma,
                          ScaledBeta,LogNormalZero,NeuGamma,FGMBesselK,FixedSpikes,WithLethal}) = false

# -----------------------------------------------------------------------------
# Parameter names: the scheme, GRAPES' names, and the names used before 2026-10-15
# -----------------------------------------------------------------------------

# Every name dfe_ml used before the `<component>_<quantity>` scheme (2026-10-15), mapped to
# the current one. GRAPES' own names are among them. Results written before that date
# (Task 2 and 3 outputs, the frozen test fixtures) use the old names.
const LEGACY_PARAM_NAMES = Dict{String,Symbol}(
    "negGshape" => :del_shape, "negGmean" => :del_mean,
    "pos_prop" => :ben_prop, "posGmean" => :ben_mean, "posGshape" => :ben_shape,
    "s0" => :del_shift, "negLNmean" => :del_mean, "negLNsigma" => :del_sigma,
    "p_neu" => :neu_prop, "Ba" => :beta_a, "Bb" => :beta_b, "med_prop" => :beta_mass,
    "BKm" => :fgm_m, "BKnsz2" => :fgm_nsz2, "BKsigma" => :fgm_sigma, "BKscale" => :fgm_scale,
    "p_lethal" => :lethal_prop, "syno_orientation_error" => :misorient,
    "a" => :nb_dispersion, "Ne_ratio" => :ne_ratio)

"""
    from_legacy_name(name) -> Symbol

The current dfe_ml parameter name for an old (or GRAPES) one: `"negGmean"` ->
`:del_mean`, `"w3"` -> `:spike_w3`. Current names map to themselves.
"""
function from_legacy_name(name::Union{AbstractString,Symbol})
    s = String(name)
    haskey(LEGACY_PARAM_NAMES, s) && return LEGACY_PARAM_NAMES[s]
    m = match(r"^w(\d+)$", s)
    m !== nothing && return Symbol("spike_w", m[1])
    return Symbol(s)
end

"""
    grapes_param_name(name) -> String

GRAPES' name of a dfe_ml parameter (for `MKtest.grapes()` output columns `"<Model>:<name>"`
and `-fixed_param` files): `:del_mean` -> `"negGmean"`, `:misorient` ->
`"syno_orientation_error"`. Parameters GRAPES does not have keep their dfe_ml name.
"""
grapes_param_name(name::Union{AbstractString,Symbol}) = get(_TO_GRAPES, Symbol(name), String(name))

const _TO_GRAPES = Dict{Symbol,String}(
    :del_shape => "negGshape", :del_mean => "negGmean", :ben_prop => "pos_prop",
    :ben_mean => "posGmean", :ben_shape => "posGshape", :del_shift => "s0",
    :beta_a => "Ba", :beta_b => "Bb", :beta_mass => "med_prop",
    :fgm_m => "BKm", :fgm_nsz2 => "BKnsz2", :fgm_sigma => "BKsigma", :fgm_scale => "BKscale",
    :lethal_prop => "p_lethal", :misorient => "syno_orientation_error")

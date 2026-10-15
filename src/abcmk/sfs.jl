################################
######    Polymorphism    ######
################################
# Expected number of polymorphism above frequency x from the standard diffusion theory
# f(x) = ∫(s) θs* (1/x(1-x)) * ( (ℯ^(Ns)*(1-ℯ^(-4Ns(1-x))) ) / (ℯ^(4Ns-1)))
# Convolote with the binomial to obtain the downsampled sfs
# E[P(x)] = ∑_x=x*→x=1 fB(x*)
############Neutral#############

"""
    sfs_neut(param, binom) -> Vector{Float64}

Expected synonymous (neutral) SFS of the sample, at derived allele counts `1:nn-1`,
under background selection (`Ne = B N`):

```math
\\mathbb{E}[P_s(i)] = B\\, \\theta_{coding}\\, p_{syn} \\sum_{j} \\mathrm{binom}[i, j] \\frac{1}{j}
```

`binom` is the projection matrix of [`binom_op`](@ref) for `param.B`. Per-frequency
counts, not cumulative.
"""
function sfs_neut(param::parameters, binom::SparseMatrixCSC{Float64,Int64})
    # binom * (1/j) is cached by binom_op; computed here if it is missing.
    binom_nn = size(binom, 1)  # = nn - 1 rows
    neut_key = (param.NN, binom_nn, param.B)
    if haskey(_NEUT_CACHE, neut_key)
        return param.B * param.θ_coding * param.syn_nonsyn_ratio .* _NEUT_CACHE[neut_key]
    end

    NN2 = ceil(Int, param.NN * param.B)
    x = 0:NN2
    solved_neutral_sfs = 1 ./ x
    replace!(solved_neutral_sfs, Inf => 0.0)
    v = binom * solved_neutral_sfs
    _NEUT_CACHE[neut_key] = v
    return param.B * param.θ_coding * param.syn_nonsyn_ratio .* v
end

############Positive############
# Variable gamma in function changed to s to avoid problem with exported SpecialFunctions.gamma
"""
    sfs_pos(param, s::Int64, p::Float64, binom) -> Vector{Float64}

Expected SFS of beneficial nonsynonymous mutations with coefficient `s` (`2N` times the
heterozygote coefficient), at derived allele counts `1:nn-1`, under background selection.

Background selection reduces `Ne` to `B N`: the spectrum is scaled by `B`, like
[`sfs_neut`](@ref) and [`sfs_neg`](@ref), the sojourn time is evaluated at `γ B`, and
[`Φ`](@ref) multiplies the whole spectrum. `p` is the probability that a new mutation
belongs to this class (`ppos_l` or `ppos_h`); `p = 0` returns zeros. `binom` is the
projection matrix of [`binom_op`](@ref) for `param.B`.
"""
function sfs_pos(
    param::parameters,
    s::Int64,
    p::Float64,
    binom::SparseMatrixCSC{Float64,Int64},
)
    if p == 0.0
        out = zeros(Float64, param.nn + 1)
        out = out[2:(end-1)]
    else
        red_plus = Φ(param, s)

        NN2 = ceil(Int, param.NN * param.B)
        xa1 = 0:NN2
        xa2 = @. xa1 / NN2

        g = s * param.B                       # sojourn at γ·B
        solved_positive_sfs::Vector{Float64} = @. (1.0 / (NN2)) * (p * 0.5 * positive_sfs(xa2, g))
        replace!(solved_positive_sfs, NaN => 0.0)

        out::Vector{Float64} =
            param.B * (param.θ_coding) * red_plus * (1 - param.syn_nonsyn_ratio) * (binom * solved_positive_sfs)
    end

    return out
end

function positive_sfs(x::Float64, gam::Real)

    if isinf(exp(gam * 2))
        # For large gam, approximate to avoid overflow
        numerator = 1 - exp(-2 * gam * (1 - x))
        denominator = x * (1 - x)
        return (numerator / denominator)
    else
        exp_2gam = exp(2 * gam)
        numerator = exp_2gam * (1 - exp(-2 * gam * (1 - x)))
        denominator = (exp_2gam - 1) * x * (1 - x)
        return (numerator / denominator)
    end
end


# Population-frequency vectors for rates(): sfs_pos and sfs_neg before the binomial
# projection and the scalar factor. rates() applies `TailWeights` to them. Keep the
# arithmetic order: test/fixtures pins the results bit for bit.

"""
    _src_pos!(scratch, param, g, p, ncols) -> view

Beneficial source vector of length `ncols = NN2 + 1`, written into `scratch`. `g` is the
coefficient of the sojourn time, `s·B`. The caller applies the projection and the factor
`B * θ_coding * Φ(param, s) * (1 - syn_nonsyn_ratio)`.
"""
function _src_pos!(
    scratch::AbstractVector{Float64},
    param::parameters,
    g::Real,
    p::Float64,
    ncols::Int,
)
    NN2 = ncols - 1
    inv_NN2 = 1.0 / NN2
    coef = p * 0.5 * inv_NN2
    src = view(scratch, 1:ncols)
    if p == 0.0
        fill!(src, 0.0)
        return src
    end
    # x = 0 and x = 1 are zero, as in sfs_pos.
    @inbounds src[1] = 0.0
    @inbounds src[ncols] = 0.0
    @inbounds for i in 2:NN2
        x = (i - 1) * inv_NN2
        src[i] = coef * positive_sfs(x, g)
    end
    return src
end

"""
    _src_neg!(scratch, param, p, ncols) -> view

Deleterious source vector of length `ncols = NN2 + 1`, as `_src_pos!`. The caller applies
the projection and the factor `B * θ_coding * (1 - syn_nonsyn_ratio)`. The loop matches
`sfs_neg`, including its early stop, so both give the same values.
"""
function _src_neg!(
    scratch::AbstractVector{Float64},
    param::parameters,
    p::Float64,
    ncols::Int,
)
    NN2 = ncols - 1
    beta = param.scale / (1.0 * param.B)
    zeta_const = zeta(param.shape, (2.0 + beta) / 2.0)
    coeff = (1.0 - p) * (2.0^-param.shape) * (beta^param.shape)
    inv_NN2 = 1.0 / NN2

    src = view(scratch, 1:ncols)
    fill!(src, 0.0)

    peak = 0.0
    tol = 1e-12
    @inbounds for i in 2:NN2
        x = (i - 1) * inv_NN2
        val = coeff * (-zeta(param.shape, x + beta / 2.0) + zeta_const) /
              ((-1.0 + x) * x)
        src[i] = val
        abs_val = abs(val)
        if abs_val > peak
            peak = abs_val
        elseif peak > 0.0 && abs_val < tol * peak
            break
        end
    end
    @inbounds @simd for i in eachindex(src)
        src[i] *= inv_NN2
    end
    return src
end

######Slightly deleterious######
"""
    sfs_neg(param, p::Float64, binom) -> Vector{Float64}

Expected SFS of deleterious nonsynonymous mutations, integrated over the gamma DFE
(`param.shape`, `param.scale / B`), at derived allele counts `1:nn-1`, under background
selection (`Ne = B N`).

`p` is the probability that a new nonsynonymous mutation is beneficial, so `1 - p` are
deleterious (callers pass `ppos_l + ppos_h`). The sum over population frequencies stops
once a term falls below `1e-12` of the largest one. `binom` is the projection matrix of
[`binom_op`](@ref) for `param.B`.
"""
function sfs_neg(param::parameters, p::Float64, binom::SparseMatrixCSC{Float64,Int64})
    beta = param.scale / (1.0 * param.B)
    NN2 = ceil(Int, param.NN * param.B)


    zeta_const = zeta(param.shape, (2.0 + beta) / 2.0)
    coeff = (1.0 - p) * (2.0^-param.shape) * (beta^param.shape)

    xa1 = 0:NN2
    xa2 = @. xa1 / NN2

    solve_z = zeros(Float64, length(xa2))

    peak = 0.0
    tol = 1e-12

    @inbounds for i in 2:(length(xa2)-1)
        x = xa2[i]
        val = coeff * (-zeta(param.shape, x + beta / 2.0) + zeta_const) /
              ((-1.0 + x) * x)
        solve_z[i] = val

        # Stop once a term falls below 1e-12 of the peak.
        abs_val = abs(val)
        if abs_val > peak
            peak = abs_val
        elseif peak > 0.0 && abs_val < tol * peak
            break
        end
    end
    # x = 0 and x = 1 stay zero.

    solved_negative = (1.0 / NN2) .* solve_z

    out = param.B * (param.θ_coding) * (1 - param.syn_nonsyn_ratio) * (binom * solved_negative)

    return out
end

################################
######     Fixations      ######
################################

# E[Dn]  = LT(E[Dn+] + E[Dn-] + E[Dns])

# Neutral fixation corrected by background selection
"""
    fix_neut(param) -> Float64

Expected rate of synonymous (neutral) fixations per new coding mutation: the
synonymous fraction times the neutral fixation probability at `Ne = B N`,

```math
\\mathbb{E}[D_{s}] = p_{syn} \\frac{1}{2 B N}
```

with `p_syn = param.f_syn`. `rates` and `analytical_alpha` multiply it by `B`.
"""
function fix_neut(param::parameters)
    # Synonymous probabilty * (fixation probabilty corrected by BGS value)
    out::Float64 = (param.f_syn) * (1.0 / (param.B * param.NN))
    return out
end

# Negative fixations corrected by background selection
"""
    fix_neg(param, ppos::Float64) -> Float64

Expected rate of deleterious nonsynonymous fixations per new coding mutation, integrated
over the gamma DFE (`param.shape`, `param.scale`) at `Ne = B N`:

```math
\\mathbb{E}[D_{n-}] = (1 - p_{syn})(1 - p_{+})\\, 2^{-a} B^{-a} b^{a}
\\left(\\zeta\\left(a, \\tfrac{1}{2}\\left(2 - \\tfrac{1}{N B} + \\tfrac{b}{B}\\right)\\right) - \\zeta\\left(a, 1 + \\tfrac{b}{2B}\\right)\\right)
```

with `a` the shape, `b` the scale, `ζ` the Hurwitz zeta function and `p₊ = ppos` the
probability that a nonsynonymous mutation is beneficial (callers pass
`0.5 p_weak + 0.5 p_strong`).
"""
function fix_neg(param::parameters, ppos::Float64)
    # Non-synonymous proportion * negative alleles probability * fixation probability from gamma distribution
    out::Float64 =
        (1 - param.f_syn) *
        (1 - ppos) *
        (2^(-param.shape)) *
        (param.B^(-param.shape)) *
        (param.scale^param.shape) *
        (
            -zeta(param.shape, 1.0 + param.scale / (2.0 * param.B)) + zeta(
                param.shape,
                0.5 * (2 - 1.0 / (param.N * param.B) + param.scale / param.B),
            )
        )
    return out
end

# Positive fixations
"""
    p_fix(param, gam::Int64) -> Float64

Fixation probability of a new beneficial mutation with coefficient `gam` (`2N` times the
heterozygote coefficient, `s = gam/NN`):

```math
p_{fix} = \\frac{1 - e^{-2s}}{1 - e^{-2\\gamma}}
```

For `s >= 0.1` it uses instead the survival probability of a branching process with
Poisson(`1 + s`) offspring (200 fixed-point iterations), following Uricchio et al. 2014.
Background selection is not included: see [`fix_pos_sim`](@ref).
"""
function p_fix(param::parameters, gam::Int64)

    # Fixation probability
    s::Float64 = gam / (param.NN)
    pfix::Float64 = (1.0 - ℯ^(-2.0 * s)) / (1.0 - ℯ^(-2.0 * gam))

    # Correcting p_fix for large s following Uricchio et al. 2014
    if s >= 0.1
        pfix = ℯ^(-(1.0 + s))
        lim = 0
        while (lim < 200)
            pfix = ℯ^((1.0 + s) * (pfix - 1.0))
            lim = lim + 1
        end
        pfix = 1 - pfix
    end

    return pfix
end

# Positive fixations after apply Φ. reduction of positive fixations due deleterious linkage given a value B of background selection
"""
    fix_pos_sim(param, gamma::Int64, ppos::Float64) -> Float64

Expected rate of beneficial fixations per new coding mutation for the class with
coefficient `gamma`, reduced by background selection and linkage:

```math
\\mathbb{E}[D_{n+}] = (1 - p_{syn})\\, p_{+}\\, \\Phi(\\gamma)\\, p_{fix}(\\gamma)
```

with `p₊ = ppos` (callers pass `0.5 p_weak` or `0.5 p_strong`), [`Φ`](@ref) the linkage
reduction and [`p_fix`](@ref) the fixation probability.
"""
function fix_pos_sim(param::parameters, gamma::Int64, ppos::Float64)
    red_plus = Φ(param, gamma)
    # Non-synonymous * positive alleles probability * B reduction * fixation probility
    out::Float64 = (1 - param.f_syn) * ppos * red_plus * p_fix(param, gamma)
    return out
end

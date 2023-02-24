################################
####### Define parameters ######
################################

"""
Mutable structure containing the variables required to solve the analytical approach. All the functions are solve using the internal values of the structure. You should declare a mutable structure to the perform the analytical estimations.

# Parameters
 - `gam_flanking::Int64`: Selection coefficient for deleterious alleles
 - `gL::Int64`: Selection coefficient for weakly benefical alleles
 - `gH::Int64`: Selection coefficient for strongly benefical alleles
 - `al_low::Float64`: Proportion of α due to weak selection
 - `al_tot::Float64`: α
 - `θ_flanking::Float64`: Mutation rate defining BGS strength
 - `θ_coding::Float64`: Mutation rate on coding region
 - `al::Float64`: DFE shape parameter 
 - `be::Float64`: DFE scale parameter
 - `B::Float64`: BGS strength
 - `B_bins::Array{Float64,1}`: BGS values to simulate
 - `ppos_l::Float64`: Fixation probabily of weakly beneficial alleles
 - `ppos_h::Float64`: Fixation probabily of strongly beneficial alleles
 - `N::Int64`: Population size
 - `n::Int64`: Sample size
 - `Lf::Int64`: Flanking region length
 - `ρ::Float64`: Recombination rate
 - `TE::Float64`
"""
@with_kw mutable struct parameters
    gL::Int64 = 10
    gH::Int64 = 500
    al_low::Float64 = 0.2
    al_tot::Float64 = 0.4
    θ_flanking::Float64 = 1e-3
    gam_flanking::Int64 = -457
    θ_coding::Float64 = 1e-3
    shape::Float64 = 0.184
    gam_dfe::Int64 = -457
    scale::Float64 = shape / abs(gam_dfe)
    B::Float64 = 0.999
    B_bins::Array{Float64,1} = push!(collect(0.1:0.025:0.975), 0.999)
    ppos_l::Float64 = 0
    ppos_h::Float64 = 0
    N::Int64 = 1000
    n::Int64 = 500
    Lf::Int64 = 2 * 10^5
    ρ::Float64 = 0.001
    TE::Float64 = 5.0
    diploid::Bool = false
    cutoff::Vector{Float64} = [0.0, 1.0]
    isolines::Bool = false

    NN::Int64 = 2 * N
    nn::Int64 = 2 * n
    dac::Array{Int64,1} = [2, 4, 5, 10, 20, 50, 200, 500, 700]
end

function assertion_params(param::parameters)
    @assert (param.cutoff[1] >= 0.0) & (param.cutoff[end] <= 1.0) "Frequency cutoff must be at the interval [0,1]"
    @assert all(param.dac .< param.nn) "Selected DAC is larger than sample size"

    freq_dac = hcat(round.(collect(1:(param.nn-1)) ./ param.nn, digits = 4), 1:(param.nn-1))
    freq_dac = freq_dac[
        (freq_dac[:, 1].>=param.cutoff[1]).&(freq_dac[:, 1].<=param.cutoff[end]),
        :,
    ]

    @assert all(in(freq_dac[:, 2]).(param.dac)) "Please select a DAC according to your frequency cutoff"
end

################################
###### Solving parameters ######
################################
"""
	Br(Lmax,theta)

Expected reduction in nucleotide diversity. Explored at [Charlesworth B., 1994](https://doi.org/10.1017/S0016672300032365):

```math
\\frac{\\pi}{\\pi_{0}} = e^{\\frac{4\\mu L}{2rL+t}}
```

# Arguments
 - `param::parameters`
 - `theta::Float64`
# Returns
 - `Float64`: expected reduction in diversity given a non-coding length, mutation rate and defined recombination.
"""
function Br(param::parameters, theta::Float64)
    ρ::Float64 = param.ρ
    t::Float64 = -1.0 * param.gam_flanking / (param.NN + 0.0)
    μ::Float64 = theta / (2.0 * param.NN)
    r::Float64 = ρ / (2.0 * param.NN)

    out::Float64 = ℯ^(-4 * μ * param.Lf / (2 * param.Lf * r + t))
    return out
end

# Set mutation rate given the expected reduction in nucleotide diversity (B value) in a locus.
"""
	set_θ!(param)

Find the optimum mutation given the expected reduction in nucleotide diversity (B value) in a locus.

# Returns
 - `adap.θ_flanking::Float64`: changes adap.θ_flanking value.
"""
function set_θ!(param::parameters)
    i(θ, p = param) = Br(p, θ) - p.B
    θ_flanking = find_zero(i, 0.0)
    param.θ_flanking = θ_flanking
end

function alpha_exp_sim_low(param::parameters, ppos_l::Float64, ppos_h::Float64)
    return fix_pos_sim(param, param.gL, 0.5 * ppos_l) / (
        fix_pos_sim(param, param.gL, 0.5 * ppos_l) +
        fix_pos_sim(param, param.gH, 0.5 * ppos_h) +
        fix_neg(param, 0.5 * ppos_l + 0.5 * ppos_h)
    )
end

function alpha_exp_sim_tot(param::parameters, ppos_l::Float64, ppos_h::Float64)
    return (
        fix_pos_sim(param, param.gL, 0.5 * ppos_l) +
        fix_pos_sim(param, param.gH, 0.5 * ppos_h)
    ) / (
        fix_pos_sim(param, param.gL, 0.5 * ppos_l) +
        fix_pos_sim(param, param.gH, 0.5 * ppos_h) +
        fix_neg(param, 0.5 * ppos_l + 0.5 * ppos_h)
    )
end

function solv_eqns(param, config)
    ppos_l, ppos_h = config
    return (
        alpha_exp_sim_tot(param, ppos_l, ppos_h) - param.al_tot,
        alpha_exp_sim_low(param, ppos_l, ppos_h) - param.al_low,
    )
end

"""
	set_ppos!(param)

Find the probabilty of positive selected alleles given the model. It solves a equation system taking into account fixations probabilities of weakly and strong beneficial alleles.

# Returns
 - `Tuple{Float64,Float64}`: weakly and strong beneficial alleles probabilites.
"""
function set_ppos!(param::parameters)
    function f!(F, x, param = param)
        F[1] = alpha_exp_sim_tot(param, x[1], x[2]) - param.al_tot
        F[2] = alpha_exp_sim_low(param, x[1], x[2]) - param.al_low
    end

    ppos_l, ppos_h = nlsolve(f!, [0.0; 0.0]).zero

    if ppos_l < 0.0
        ppos_l = 0.0
    end
    if ppos_h < 0.0
        ppos_h = 0.0
    end

    param.ppos_l, param.ppos_h = ppos_l, ppos_h
end

"""
	binom_op(param)

Binomial convolution to sample the allele frequencies probabilites depending on background selection values and sample size.

# Arguments
 - `param::parameters`
 - `binom::binomialDict`

# Returns
 - `Array{Float64,2}`: convoluted SFS for each B value defined in the model (param.B_bins). The estimations are saved at *convolutedBn.bn*.
"""
function binom_op(param::parameters)
    bn = Dict(
        param.B_bins[i] => spzeros(param.nn + 1, param.NN) for i = 1:length(param.B_bins)
    )

    for b in param.B_bins
        NN2 = convert(Int64, ceil(param.NN * b))
        samples = collect(1:(param.nn-1))
        p_size = collect(0:NN2)
        samples_freqs = permutedims(p_size / NN2)
        neutral_sfs = @. 1 / p_size
        replace!(neutral_sfs, Inf => 0.0)

        f(x, y = param.nn) = Binomial(y, x)
        z = f.(samples_freqs)

        out = pdf.(z, samples)
        out = round.(out, digits = 10)
        out_sparse = SparseArrays.dropzeros(SparseArrays.sparse(out))
        bn[b] = out_sparse
    end
    return (bn)
end

function binom_op(NN, nn, B)
    NN2 = convert(Int64, ceil(NN * B))
    samples = collect(1:(nn-1))
    p_size = collect(0:NN2)
    samples_freqs = permutedims(p_size / NN2)

    f(x, y = nn) = Binomial(y, x)
    z = f.(samples_freqs)

    out = pdf.(z, samples)
    out = round.(out, digits = 10)
    out_sparse = SparseArrays.dropzeros(SparseArrays.sparse(out))
    return (out_sparse)
end

"""

	Φ(param,s)


Reduction in fixation probabilty due to background selection and linkage. The formulas used have been subjected to several theoretical works ([Charlesworth B., 1994](https://doi.org/10.1017/S0016672300032365), [Hudson et al., 1995](https://www.genetics.org/content/141/4/1605), [Nordborg et al. 1995](https://doi.org/10.1017/S0016672300033619), [Barton NH., 1995](https://www.genetics.org/content/140/2/821)).

The fixation probabilty of selected alleles are reduce by a factor ``\\phi``:

```math
\\phi(t,s) = e^{[\\frac{-2\\mu}{t(1+\\frac{rL}{t}+\\frac{2s}{t})}]}
```

Multiplying across all deleterious linkes sites, we find:

```math
\\Phi = \\prod_{1}^{L} = \\phi(t,s) = e^{\\frac{-2t\\mu(\\psi[1,\\frac{r+2s+L}{r}] - \\psi[1,\\frac{r(L+1)+2s+t}{r}])}{r^2}}
```

```math
\\phi(t,s) = e^{\\frac{-2t\\mu(\\psi[1,\\frac{r+2s+L}{r}] - \\psi[1,\\frac{r(L+1)+2s+t}{r}])}{r^2}}
```

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model
 - `gamma::Int64`: selection coefficient.

# Returns
 - `Float64`: expected rate of positive fixations under background selection.

"""
function Φ(param::parameters, s::Int64)
    S::Float64 = abs(param.gam_flanking / (1.0 * param.NN))
    r::Float64 = param.ρ / (2.0 * param.NN)
    μ::Float64 = param.θ_flanking / (2.0 * param.NN)
    s::Float64 = s / (param.NN * 1.0)

    Ψ0::Float64 = polygamma(1, (s + S) / r)
    Ψ1::Float64 = polygamma(1, (r + param.Lf * r + s + S) / r)
    CC::Float64 = 1.0

    out::Float64 = (ℯ^(-2.0 * S * μ * (Ψ0 - Ψ1) / (r^2)))
    return out
end

"""
	analytical_alpha(param)

Analytical α(x) estimation. Solve α(x) generally. We used the expected fixation rates and SFSto approach the asympotic value accouting for background selection, weakly and strong positive selection. α(x) can be estimated taking into account the role of positive selected alleles or not. In this way we explore the role of linkage to deleterious alleles in the coding region.

```math
\\mathbb{E}[\\alpha_{x}] =  1 - \\left(\\frac{\\mathbb{E}[D_{s}]}{\\mathbb{E}[D_{N}]}\\frac{\\mathbb{E}[P_{N}]}{\\mathbb{E}[P_{S}]}\\right)
```

# Arguments
 - `param::parameters`
# Returns
 - `Tuple{Vector{Float64},Vector{Float64}}`: α(x) accounting for weak adaptation, α(x) non-accouting for weak adaptation.
"""
function analytical_alpha(param::parameters)
    binom = binom_op(param.NN, param.nn, param.B)
    ################################################################
    # Solve the model similarly to original python mktest  #	
    ################################################################
    B = param.B

    set_θ!(param)
    θ_flanking = param.θ_flanking
    # Solve the probabilities of fixations without BGS
    ## Set non-bgs
    param.B = 0.999
    ## Solve the mutation rate
    set_θ!(param)
    ## Solve the probabilities
    set_ppos!(param)
    # Return to the original values
    param.θ_flanking = θ_flanking
    param.B = B

    ################################################################
    # Accounting for positive alleles segregating due to linkage # 
    ################################################################

    # Fixation
    fN = param.B * fix_neut(param)
    fNeg = param.B * fix_neg(param, 0.5 * param.ppos_h + 0.5 * param.ppos_l)
    fPosL = fix_pos_sim(param, param.gL, 0.5 * param.ppos_l)
    fPosH = fix_pos_sim(param, param.gH, 0.5 * param.ppos_h)

    ds = fN
    dn = fNeg + fPosL + fPosH

    ## Polymorphism
    neut = sfs_neut(param, binom)

    selH::Array{Float64,1} = if isinf(exp(param.gH * 2))
        sfs_pos_float(param, param.gH, param.ppos_h, binom)
    else
        sfs_pos(param, param.gH, param.ppos_h, binom)
    end

    selL::Array{Float64,1} = sfs_pos(param, param.gL, param.ppos_l, binom)
    selN::Array{Float64,1} = sfs_neg(param, param.ppos_h + param.ppos_l, binom)

    function split_columns(matrix::Array{Float64,2})
        (view(matrix, :, i) for i = 1:size(matrix, 2))
    end
    tmp = cumulative_sfs(hcat(neut, selH, selL, selN), false)

    neut, selH, selL, selN = split_columns(tmp)
    sel = (selH + selL) + selN

    # ps = @. neut / (sel+neut)
    # pn = @. sel / (sel+neut)

    ## Outputs
    α = @. 1 - ((ds / dn) * (sel / neut))

    #############################################################
    #		Accounting for neutral and deleterious alleles		#	#############################################################
    ## Fixation
    fN_nopos = fN * (param.θ_coding / 2.0) * param.TE * param.NN
    fNeg_nopos = fNeg * (param.θ_coding / 2.0) * param.TE * param.NN
    fPosL_nopos = fPosL * (param.θ_coding / 2.0) * param.TE * param.NN
    fPosH_nopos = fPosH * (param.θ_coding / 2.0) * param.TE * param.NN

    ds_nopos = fN_nopos
    dn_nopos = fNeg_nopos + fPosL_nopos + fPosH_nopos

    ## Polymorphism
    sel_nopos = selN
    ps_nopos = @. neut / (sel_nopos + neut)
    pn_nopos = @. sel_nopos / (sel_nopos + neut)

    α_nopos = 1 .- (ds_nopos / dn_nopos) .* (sel_nopos ./ neut)

    ##########
    # Output #
    ##########
    return (α, α_nopos)
end

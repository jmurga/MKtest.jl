################################
######    Polymorphism    ######
################################
# Expected number of polymorphism above frequency x from the standard diffusion theory
# f(x) = ∫(s) θs* (1/x(1-x)) * ( (ℯ^(Ns)*(1-ℯ^(-4Ns(1-x))) ) / (ℯ^(4Ns-1)))
# Convolote with the binomial to obtain the downsampled sfs
# E[P(x)] = ∑_x=x*→x=1 fB(x*)
############Neutral#############

"""

	sfs_neut(param,binom)

Expected rate of neutral allele frequency reduce by background selection. The spectrum depends on the number of individual []

```math
\\mathbb{E}[Ps_{(x)}] = \\sum{x^{*}=x}{x^{*}=1}f_{B}(x)
```
# Input
 - `param::parameters`: mutable structure containing the model.
 - `SparseMatrixCSC{Float64, Int64}`: convoluted SFS given the B value.
# Return:
 - `Vector{Float64}`: expected rate of neutral alleles frequencies.
"""
function sfs_neut(param::parameters, binom::SparseMatrixCSC{Float64,Int64})
    NN2 = ceil(Int,param.NN * param.B)

    # Allocating variables
    x = 0:NN2
    solved_neutral_sfs = 1 ./ x
    replace!(solved_neutral_sfs, Inf => 0.0)

    out::Vector{Float64} = param.B * (param.θ_coding) * 0.25 * (binom * solved_neutral_sfs)

    return out
end

############Positive############
# Variable gamma in function changed to s to avoid problem with exported SpecialFunctions.gamma
"""

	sfs_pos(param,s,p,binom)

Expected rate of positive selected allele frequency reduce by background selection. The spectrum depends on the number of individuals.

# Arguments
 - `param::parameters`: mutable structure containing the model.
 - `s::Int64`: selection coefficient.
 - `p::Float64`: positive selected alleles probabilty.
 - `binom::SparseMatrixCSC{Float64, Int64}`: convoluted SFS given the B value.
# Return:
 - `Array{Float64}`: expected positive selected alleles frequencies.
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

        # Solving sfs
        NN2 = ceil(Int,param.NN * param.B)
        xa1 = 0:NN2
        xa2 = @. xa1 / NN2

        # Solving float precision performance using exponential rule. Only one BigFloat estimation.
        s_corrected = s * param.B

        s_exp1 = exp(s_corrected * 2)
        s_exp2 = exp(s_corrected * -2)

        function positive_sfs(
            i::Float64,
            g1::Float64 = s_exp1,
            g2::Float64 = s_exp2,
            p::Float64 = p,
        )
            Float64(p * 0.5 * (g1 * (1 - g2^(1.0 - i)) / ((g1 - 1.0) * i * (1.0 - i))))
        end

        # Original
        # p*0.5*(ℯ^(2*s_corrected)*(1-ℯ^(-2.0*s_corrected*(1.0-i)))/((ℯ^(2*s_corrected)-1.0)*i*(1.0-i)))

        # Allocating outputs
        solved_positive_sfs::Vector{Float64} = (1.0 / (NN2)) * (positive_sfs.(xa2))
        replace!(solved_positive_sfs, NaN => 0.0)

        out::Vector{Float64} =
            (param.θ_coding) * red_plus * 0.75 * (binom * solved_positive_sfs)
    end

    return out
end

function sfs_pos_float(
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

        # Solving sfs
        NN2 = ceil(Int,param.NN * param.B)
        xa1 = 0:NN2
        xa2 = @. xa1 / NN2

        s_corrected = s * param.B
        s_exp1::Quadmath.Float128 = exp(Quadmath.Float128(s_corrected * 2))
        s_exp2::Quadmath.Float128 = exp(Quadmath.Float128(s_corrected * -2))

        function positive_sfs(
            i::Float64,
            g1::Quadmath.Float128 = s_exp1,
            g2::Quadmath.Float128 = s_exp2,
            p::Float64 = p,
        )
            Float64(p * 0.5 * (g1 * (1 - g2^(1.0 - i)) / ((g1 - 1.0) * i * (1.0 - i))))
        end

        # Allocating outputs
        solved_positive_sfs::Vector{Float64} = (1.0 / (NN2)) * (positive_sfs.(xa2))
        replace!(solved_positive_sfs, NaN => 0.0)
        out::Vector{Float64} =
            (param.θ_coding) * red_plus * 0.75 * (binom * solved_positive_sfs)
        # out = out[2:end-1]

    end

    return out
end

######Slightly deleterious######
"""

	sfs_neg(param,p,binom)

Expected rate of positive selected allele frequency reduce by background selection. Spectrum drawn on a gamma DFE. It depends on the number of individuals.

# Arguments
 - `param::parameters`: mutable structure containing the model
 - `p::Float64`: positive selected alleles probabilty.
 - `binom::SparseMatrixCSC{Float64, Int64}`: convoluted SFS given the B value.
# Return:
 - `Array{Float64}`: expected negative selected alleles frequencies.
"""
function sfs_neg(param::parameters, p::Float64, binom::SparseMatrixCSC{Float64,Int64})
    beta = param.scale / (1.0 * param.B)
    NN2 = ceil(Int,param.NN * param.B)
    xa1 = 0:NN2
    xa2 = @. xa1 / NN2

    function z(x::Float64, p::Float64 = p)
        (1.0 - p) *
        (2.0^-param.shape) *
        (beta^param.shape) *
        (-zeta(param.shape, x + beta / 2.0) + zeta(param.shape, (2 + beta) / 2.0)) /
        ((-1.0 + x) * x)
    end

    solve_z = z.(xa2)

    if (solve_z[1] == Inf || isnan(solve_z[1]))
        solve_z[1] = 0.0
    end
    if (solve_z[lastindex(solve_z)] == Inf || isnan(solve_z[lastindex(solve_z)]))
        solve_z[lastindex(solve_z)] = 0.0
    end

    solved_negative = 1.0 / (NN2 + 0.0) .* solve_z

    out = param.B * (param.θ_coding) * 0.75 * (binom * solved_negative)

    return out
end

"""
	cumulative_sfs(sfs_tmp)

Changing SFS considering all values above a frequency *x*. The original aMK approach takes Pn(x) and Ps(x) as the number of polymorphic sites at frequency *x* rather than above *x*, but this approach scales poorly as sample size increases. We define the polymorphic spectrum as stated above since these quantities trivially have the same asymptote but are less affected by changing sample size.

# Arguments
 - `sfs_tmp::Union{Array,SubArray}`: SFS vector.
 - `freqs::Bool`: true/false wheter the input SFS vector containing a first columns with the frequencies.
# Output 
 - `Vector{Float64}`: SFS vector.
"""#=
function cumulative_sfs!(sfs_tmp::Matrix{Float64},freqs::Bool=true)

    idx = ifelse(freqs, 2, 1)

    for c ∈ eachcol(@view sfs_tmp[:,2:end])
        cumulative_vector!(c)
    end

    return sfs_tmp
end
=#
function cumulative_vector!(v::Union{SubArray,Vector{T}}) where {T<:Number}
    n = length(v)
    @inbounds for i = 1:n
        v[i] = sum(@view v[i:end])
    end
    return v
end

function cumulative_sfs(sfs_tmp::Union{Array, SubArray}, freqs::Bool = true)
    n, m = size(sfs_tmp)
    out = similar(sfs_tmp)

    idx = ifelse(freqs, 2, 1)

    out[1, idx:end] = sum(sfs_tmp[:, idx:end], dims=1)

    @inbounds @simd for i in 2:n
        app = view(out,i-1, idx:size(out,2)) .- view(sfs_tmp,i-1, idx:size(out,2))

        if sum(app) > 0.0
            out[i, idx:end] .= app
        else
            out[i, idx:end] .= 0.0
        end
    end

    if freqs
        out[:, 1] .= sfs_tmp[:, 1]
    end

    return out
end


"""
	reduce_sfs(sfs_tmp,bins)

Function to reduce the SFS into N bins.

# Arguments
 - `sfs_tmp::Vector{Float64}`: SFS vector.
 - `bins::Int64`: bin size to collapse the SFS.
# Output 
 - `Vector{Float64}`: Binned SFS vector.
"""
function reduce_sfs(sfs_tmp::Array, bins::Int64)
    freq = collect(0:(size(sfs_tmp, 1)-1)) / size(sfs_tmp, 1)
    h1 = StatsBase.fit(StatsBase.Histogram, freq, 0:(1/(bins-1)):1)
    xmap1 = StatsBase.binindex.(Ref(h1), freq)

    tmp = hcat(sfs_tmp, xmap1)
    out = Array{Float64}(undef, bins - 1, size(sfs_tmp, 2))
    vIter = convert(Array, unique(xmap1)')
    @simd for i in eachindex(vIter)
        @inbounds out[i, 2:end] = sum(tmp[tmp[:, end].==i, 2:(end-1)], dims = 1)
    end

    out[:, 1] = collect(1:(bins-1)) ./ bins

    return (out)
end

"""
	α_x(sfs,divergence)
Function to estimate alpha_x
"""
function α_x(sfs::Matrix{Float64},divergence::Matrix{Float64};cumulative::Bool=true)

	if cumulative
		sfs = cumulative_sfs(sfs)
	end

	alpha_x = @. 1 - (sfs[:,2]/sfs[:,3] * divergence[2]/divergence[1])

	return(alpha_x)
end
"""
	aMK(x,columns)

Function to estimate the asymptotic value of α(x).

# Arguments
 - `sfs::Matrix{Float64}`: SFS array.
 - `divergence::Matrix{Float64}`: SFS array.
# Returns
 - `Array{Float64,2}`: Array of array containing asymptotic values, lower confidence interval, higher confidence interval
"""
function aMK(
    param::parameters,
    α::Vector{Float64};
    dac_rm::Bool = false,
    na_rm::Bool = false)

    @unpack nn, dac = param

    f = collect(1:length(α)) / nn


    if dac_rm
        f = dac / nn
        α = α[dac]
    end

    if na_rm
        flt = .!isnan.(α) .&& .!isinf.(α) .&& α .!= 1
        α = α[flt]
        f = f[flt]
    end

    # Model
    model(x, p) = @. p[1] + p[2] * exp(-x * p[3])

    # Fit values

    fitted1 = curve_fit(
        model,
        f,
        α,
        [-1.0, -1.0, 1.0];
        lower = [-1.0, -1.0, 1.0],
        upper = [1.0, 1.0, 10.0],
    )
    fitted2 = curve_fit(model, f, α, fitted1.param)
    asymp = model(1, fitted2.param)

    ci = try
        [confidence_interval(fitted2)[1][1], confidence_interval(fitted2)[1][2]]
    catch
        [0.0, 0.0]
    end

    return (asymp, ci, fitted2.param)
end
"""
	imputedMK(sfs,divergence,m,cutoff)

Function to estimate the imputedMK

# Arguments
 - `sfs::Matrix{Float64}`: SFS array
 - `divergence::Array`: divegence count array
 - `m::Union{Nothing,Array}`: non-synonymous and synonymous sites
 - `cutoff{Float64}`: frequency cutoff to perform imputedMK
# Returns
 - `Dict: Dictionary containing results
"""
function imputedMK(sfs::Matrix{Float64},divergence::Matrix{Int64},cutoff::Float64=0.15;m::T=nothing) where {T<:Union{Nothing,Array}}

	output = OrderedDict{String,Float64}()

	pn = sum(sfs[:,2])
	ps = sum(sfs[:,3])
	dn = divergence[1]
	ds = divergence[2]
	
	deleterious = 0
	### Estimating slightly deleterious with pn/ps ratio
	flt_low = (sfs[:, 1] .<= cutoff)
	pn_low   = sum(sfs[flt_low,2])
	ps_low   = sum(sfs[flt_low,3])

	flt_inter = (sfs[:, 1] .>= cutoff) .& (sfs[:, 1] .<= 1)
	pn_inter = sum(sfs[flt_inter,2])
	ps_inter = sum(sfs[flt_inter,3])

	ratio_ps       = ps_low / ps_inter
	deleterious   = pn_low - (pn_inter * ratio_ps)

	if (deleterious > pn) || (deleterious < 0)
		deleterious = 0
		pn_neutral = round(pn - deleterious,digits=3)
	else
		pn_neutral = round(pn - deleterious,digits=3)
	end

	output["alpha"] = round(1 - ((pn_neutral/ps) * (ds / dn)),digits=3)

	#  method = :minlike same results R, python two.sides
	output["pvalue"] = pvalue(FisherExactTest(Int(ps),Int(ceil(pn_neutral)),Int(ds),Int(dn)))


	if (!isnothing(m))

		mn = m[1]; ms = m[2]
		# ## Estimation of b: weakly deleterious
		output["b"] = (deleterious / ps) * (ms / mn)

		## Estimation of f: neutral sites
		output["f"] = (ms * pn_neutral) / (mn * ps)

		## Estimation of d, strongly deleterious sites
		output["d"] = 1 - (output["f"] + output["b"])

		ka      = dn / mn
		ks       = ds / ms
		output["omega"]    = ka/ ks

		# Omega A and Omega D
		output["omegaA"] = output["omega"] * output["alpha"]
		output["omegaD"] = output["omega"] - output["omegaA"]	    
	end

	return output
end

"""
	fwwMK(sfs,divergence,m,cutoff)

Function to estimate the Fay, Wycoff and Wu's MK extension.

# Arguments
 - `sfs::Matrix{Float64}`: SFS array
 - `divergence::Array`: divegence count array
 - `m::Union{Nothing,Array}`: non-synonymous and synonymous sites
 - `cutoff{Float64}`: frequency cutoff to perform fwwMK
# Returns
 - `Dict: Dictionary containing results
"""
function fwwMK(;sfs::Array,divergence::Array,m::T,cutoff::Float64=0.15) where {T<:Union{Nothing,Array}}

	output = OrderedDict{String,Float64}()

	ps = sum(sfs[:,2])
	pn = sum(sfs[:,3])
	dn = divergence[1]
	ds = divergence[2]
	
	
	deleterious = 0
	### Estimating slightly deleterious with pn/ps ratio
	flt_inter = (sfs[:, 1] .>= cutoff) .& (sfs[:, 1] .<= 1)
	pn_inter = sum(sfs[flt_inter,2])
	ps_inter = sum(sfs[flt_inter,3])


	# output['alpha'] = 1 - (((pn - deleterious) / ps) * (ds / dn))
	output["alpha"] = round(1 - ((pn_inter/ps_inter) * (ds / dn)),digits=3)
	#  method = :minlike same results R, python two.sides
	output["pvalue"] = pvalue(FisherExactTest(Int(ps_inter),Int(ceil(pn_inter)),Int(ds),Int(dn)))


	if (!isnothing(m))

		mn = m[1]; ms = m[2]
		ka      = dn / mn
		ks       = ds / ms
		output["omega"]    = ka/ ks


		# Omega A and Omega D
		output["omegaA"] = output["omega"] * output["alpha"]
		output["omegaD"] = output["omega"] - output["omegaA"]	    
	end

	return output
end

"""
	standardMK(x,columns)

Function to estimate the original α value

# Arguments
 - `sfs::Matrix{Float64}`: SFS array
 - `divergence::Array`: divegence count array
 - `m::Union{Nothing,Array}`: non-synonymous and synonymous sites# Returns

# Output
 - `Dict: Dictionary containing results
"""
function standardMK(;sfs::Array,divergence::Array,m::T=nothing) where {T<:Union{Nothing,Array}}

	output = OrderedDict{String,Float64}()

	pn = sum(sfs[:,2])
	ps = sum(sfs[:,3])
	dn = divergence[1]
	ds = divergence[2]
	
	
	output["alpha"] = round(1 - ((pn/ps) * (ds / dn)),digits=3)
	#  method = :mnnlike same results R, python two.sides
	output["pvalue"] = pvalue(FisherExactTest(Int(ps),Int(ceil(pn)),Int(ds),Int(dn)))

	if (!isnothing(m))

		mn = m[1]; ms = m[2]

		ka      = dn / mn
		ks       = ds / ms
		output["omega"]    = ka/ ks


		# Omega A and Omega D
		output["omegaA"] = output["omega"] * output["alpha"]
		output["omegaD"] = output["omega"] - output["omegaA"]	    
	end

	return output
end

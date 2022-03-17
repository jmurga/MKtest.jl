################################
######     Fixations      ######
################################

# E[Dn]  = LT(E[Dn+] + E[Dn-] + E[Dns])

# Neutral fixation corrected by background selection
"""

 - fix_neut(param)

Expected neutral fixations rate reduce by B value.

```math
\\mathbb{E}[D_{s}] = (1 - p_{-} - p_{+}) B \\frac{1}{2N}
```
# Input
  - `param::parameters`: mutable structure containing the model
# Returns
 - `Float64`: expected rate of neutral fixations.

"""
function fix_neut(param::parameters)
	# Synonymous probabilty * (fixation probabilty corrected by BGS value)
	out::Float64 = 0.25*(1.0/(param.B*param.NN))
	return out
end

"""

	fix_neg_b(ppos)

Expected fixation rate from negative DFE.

```math
\\mathbb{E}[D_{n-}] =  p_{-}\\left(2^-\\alpha\\beta^\\alpha\\left(-\\zeta[\\alpha,\\frac{2+\\beta}{2}] + \\zeta[\\alpha,1/2(2-\\frac{1}{N+\\beta})]\\right)\\right)
```

# Arguments
 - `param::parameters`: mutable structure containing the model
 - `ppos::Float64`: positive selected alleles probabilty.

# Returns
 - `Float64`: expected rate of fixations from negative DFE.

"""
function fix_neg_b(param::parameters,ppos::Float64)
	# Non-synonymous proportion * negative alleles probability * fixation probability from gamma distribution
	out::Float64 = 0.75*(1-ppos)*(2^(-param.al))*(param.B^(-param.al))*(param.be^param.al)*(-zeta(param.al,1.0+param.be/(2.0*param.B))+zeta(param.al,0.5*(2-1.0/(param.N*param.B)+param.be/param.B)))
	return out
end

# Positive fixations
"""
	p_fix(param,s)

Expected positive fixation rate.

```math
\\mathbb{E}[D_{n+}] =  p_{+} \\cdot B \\cdot (1 - e^{(-2s)})
```

# Arguments
 - `param::parameters`: mutable structure containing the model
 - `s::Float64`: Selection coefficient

# Returns
 - `Float64`: probability of fixation.

"""
function p_fix(param::parameters,s::Int64)

	# Fixation probability
	s_corrected::Float64    = s/(param.NN)
	pfix::Float64 = (1.0-ℯ^(-2.0*s))/(1.0-ℯ^(-2.0*s))

	# Correcting p_fix for large s following Uricchio et al. 2014
	if s_corrected >= 0.1
		pfix = ℯ^(-(1.0+s_corrected))
		lim = 0
		while(lim < 200)
			pfix = ℯ^((1.0+s_corrected)*(pfix-1.0))
			lim  = lim + 1
		end
		pfix = 1 -pfix
	end

	return pfix
end

"""

	fix_pos_sim(param,s,ppos)

Expected positive fixations rate reduced due to the impact of background selection and linkage. The probabilty of fixation of positively selected alleles is reduced by a factor Φ across all deleterious linked sites [`Analytical.Φ`](@ref).

```math
\\mathbb{E}[D_{n+}] =  \\Phi \\cdot \\mathbb{E}[D_{n+}]
```

# Arguments
 - `param::parameters`: mutable structure containing the model
 - `s::Float64`: selection coefficient
 - `ppos::Float64`: positive selected alleles probabilty.

# Returns
 - `Float64`: expected rate of positive fixations under background selection

"""
function fix_pos_sim(param::parameters,s::Int64,ppos::Float64)


	red_plus = Φ(param,s)

	# Non-synonymous * positive alleles probability * B reduction * fixation probility
	out::Float64 = 0.75*ppos*red_plus*p_fix(param,s)
	return out
end

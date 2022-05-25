"""
	rates(param::parameters,iterations::Int64,divergence::Array,sfs::Array)

Function to solve randomly *N* scenarios. The function will create *N* models, defined by ```Analytical.parameters()```, to estimate analytically fixation and polymorphic rates for each model. The rates will be used to compute summary statistics required at ABC. The function output a HDF5 file containing the solved models, the selected DAC and the analytical rates. 

If ρ and/or theta are set to ```nothing```, the function will input random values given the range 0.0005:0.0005:0.01. Otherwise you can fix the values.

If gL is set to ```nothing```, the function will not account the role of the weakly selected alleles in the estimation.

# Arguments
 - `param::parameters`: mutable structure containing the model
 - `binom::binomialDict` : structure containing the binomial convolution
 - `gH::Array{Int64,1}` : Range of strong selection coefficients
 - `gL::Union{Array{Int64,1},Nothing}`: Range of weak selection coefficients
 - `gam_neg::Array{Int64,1}` : Range of deleterious selection coefficients
 - `theta::Union{Float64,Nothing}` : Population-scaled mutation rate on coding region
 - `ρ::Union{Float64,Nothing}` : Population-scaled recombination rate
 - `shape::Float64=0.184` : DFE shape parameter
 - `iterations::Int64` : Number of solutions
 - `output::String` : File to output HDF5 file
# Returns
 - `Array`: summary statistics
 - `Output`: HDF5 file containing models solved and rates.
"""
function rates(;param::parameters,
				gH::S,
				gL::S,
				gam_neg::S,
				alpha::Vector{Float64}=[0.1,0.9],
				theta::Union{Float64,Nothing}=0.001,
				rho::Union{Float64,Nothing}=0.001,
				iterations::Int64,
				output::String,
				threads::Bool=false) where S <: Union{Array{Int64,1},UnitRange{Int64},Nothing}
	
	@assert (alpha[1] >= 0) & (alpha[end] <= 1) "α values must be at the interval [0,1]"

	assertion_params(param)
	
	# temp = mktempdir(homedir());
	# temp = temp .* "/tmp_" .* string.(1:iterations) .* ".txt"
	# Iterations = models to solve

	# Factor to modify input Γ(shape) parameter. Flexible Γ distribution over negative alleles
	fac     = rand(-2:0.05:2,iterations)
	afac    = @. param.al*(2^fac)
	
	# Deleting shape > 1. Negative alpha_x values
	idx = findall(afac .> 1)
	if !isempty(idx)
		afac[idx] = rand(afac[afac .< 1],size(idx,1))
	end

	# Random α values
	nTot    = rand(alpha[1]:0.01:alpha[end],iterations)
	
	# Defining αW. It is possible to solve non-accounting for weak fixations
	if isnothing(gL)
		# Setting αW to 0 for all estimations
		nLow    = fill(0.0,iterations)
		# Random strong selection coefficients
		ngl     = rand(repeat([1],iterations),iterations);
	else
		# Setting αW as proportion of α
		lfac    = rand(0.0:0.05:0.9,iterations)
		nLow    = @. nTot * lfac
		# Random weak selection coefficients
		ngl     = rand(repeat(gL,iterations),iterations);
	end

	# Random strong selection coefficients
	ngh     = rand(repeat(gH,iterations),iterations);
	# Random negative selection coefficients
	ngam_neg = rand(repeat(gam_neg,iterations),iterations);

	if !isnothing(theta)
		θ = fill(theta,iterations)
	else
		θ = rand(0.0005:0.0005:0.01,iterations)
	end

	# Random rho on coding regions
	if !isnothing(rho)
		ρ = fill(rho,iterations)
	else
		ρ = rand(0.0005:0.0005:0.05,iterations)
	end
	
	# Creating N models to iter in threads. Set N models (paramerters) and sampling probabilites (binomialDict)
	binom  = binom_op!(param)
	# Estimations to distributed workers
	if threads
		@time out = ThreadsX.mapi( (α_strong,α_weak,γ_strong,γ_weak,γ_neg,shape,θ_n,ρ_n) -> iter_rates(param, binom, α_strong,α_weak,γ_strong,γ_weak,γ_neg,shape,θ_n,ρ_n),nTot, nLow, ngh, ngl, ngam_neg, afac, θ, ρ);
	else
		@time m,r_ps,r_pn,r_f = unzip(ParallelUtilities.pmapbatch( (α_strong,α_weak,γ_strong,γ_weak,γ_neg,shape,θ_n,ρ_n) -> iter_rates(param, binom, α_strong,α_weak,γ_strong,γ_weak,γ_neg,shape,θ_n,ρ_n),nTot, nLow, ngh, ngl, ngam_neg, afac, θ, ρ));
	end

	# Reducing output array
	# df = vcat(out...)
	
	# # Saving models and rates
	# models = DataFrame(df[:,1:8],[:B,:al_low,:al_tot,:gam_neg,:gL,:gH,:al,:ρ])
	# neut   = df[:,9:(8+size(param.dac,1))]
	# sel    = df[:,(9+size(param.dac,1)):(8+size(param.dac,1)*2)]
	# dsdn   = Array(df[:,(end-3):end])

	df = vcat(m...)
	
	# Saving models and rates
	models = DataFrame(df,[:B,:al_low,:al_tot,:gam_neg,:gL,:gH,:al,:ρ]);
	neut   = vcat(r_ps...);
	sel    = vcat(r_ps...);
	dsdn   = vcat(r_f...);

	# Saving multiple summary statistics
	n = OrderedDict{Int,Array}()
	s = OrderedDict{Int,Array}()
	for i in eachindex(param.dac)
		n[param.dac[i]] = neut[:,i]
		s[param.dac[i]] = sel[:,i]
	end

	# Writting HDF5 file
	string_cutoff = "cutoff=[" * string(param.cutoff[1]) * "," * string(param.cutoff[end]) * "]"
	JLD2.jldopen(output, "a+") do file
		file[string(param.N)* "/" * string(param.n) * "/" * string_cutoff *  "/models"] = models;
		file[string(param.N)* "/" * string(param.n) * "/" * string_cutoff *  "/neut"]   = n;
		file[string(param.N)* "/" * string(param.n) * "/" * string_cutoff *  "/sel"]    = s;
		file[string(param.N)* "/" * string(param.n) * "/" * string_cutoff *  "/dsdn"]   = dsdn;
		file[string(param.N)* "/" * string(param.n) * "/" * string_cutoff *  "/dac"]    = param.dac;
	end;
end

"""
	iter_rates(param::parameters,afac::Float64,bfac::Float64,al_tot::Float64,al_low::Float64,divergence::Array,sfs::Array)

Estimating rates given a model for all B range.

# Arguments
 - `param::parameters`
 - `binom::Dict`
 - `al_tot::Float64`
 - `al_low::Float64`
 - `gH::Int64`
 - `gL::Int64`
 - `gam_neg::Int64`
 - `afac::Float64`
 - `ρ::Float64`
 - `θ::Float64`
# Output
 - `Array{Float64,2}`
"""
function iter_rates(param::parameters,binom::Dict{Float64, SparseMatrixCSC{Float64, Int64}},al_tot::Float64,al_low::Float64,gH::Int64,gL::Int64,gam_neg::Int64,afac::Float64,θ::Float64,ρ::Float64)

	# Creating model to solve
	# Γ distribution
	param.al    = afac; param.be = abs(afac/gam_neg); param.gam_neg = gam_neg
	# α, αW
	param.al_low = al_low; param.al_tot = al_tot;
	# Positive selection coefficients
	param.gH    = gH;param.gL = gL
	# Mutation rate and recomb
	param.θ_coding = θ; param.ρ = ρ
	#=param.θ_coding = θ; param.θᵣ .= θᵣ; param.ρ = ρ=#
	# Solving θ on non-coding region and probabilites to get α value without BGS
	param.B = 0.999
	set_θ!(param)
	set_ppos!(param)

	# Allocate array to solve the model for all B values
	m = zeros(size(param.B_bins,1),8)
	r_ps = zeros(size(param.B_bins,1),size(param.dac,1))
	r_pn = zeros(size(param.B_bins,1),size(param.dac,1))
	r_f = zeros(size(param.B_bins,1),4)
	# r = zeros(size(param.B_bins,1),(size(param.dac,1) * 2) + 12)
	for j in eachindex(param.B_bins)
		# Set B value
		param.B = param.B_bins[j]
		# Solve θ non-coding for the B value.
		set_θ!(param)
		# Solve model for the B value
		x,y,z,w = try
			getting_rates(param,binom[param.B])
		catch e
			zeros(size(param.dac,1) *2+ 12)'
		end
		@inbounds m[j,:] = x
		@inbounds r_ps[j,:] = y
		@inbounds r_pn[j,:] = z
		@inbounds r_f[j,:] = w
	end

	return (m,r_ps,r_pn,r_f)
end

"""
	getting_rates(gammaL,gammaH,ppos_l,ppos_h,observedData,nopos)

Estimating analytical rates of fixation and polymorphism to approach α value accouting for background selection, weakly and strong positive selection. Output values will be used to sample from a Poisson distribution the total counts of polymorphism and divergence using observed data. 

# Arguments
 - `param::parameters`
 - `binom::SparseMatrixCSC{Float64,Int64}`
# Returns
 - `Array{Float64,2}` containing solved model, fixation and polymorphic rates
"""
function getting_rates(param::parameters,binom::SparseMatrixCSC{Float64,Int64})

	################################################
	# Subset rates accounting for positive alleles #
	################################################

	# Fixation
	fN       = param.B*fix_neut(param)
	fNeg     = param.B*fix_neg(param,0.5*param.ppos_h+0.5*param.ppos_l)
	fPosL    = fix_pos_sim(param,param.gL,0.5*param.ppos_l)
	fPosH    = fix_pos_sim(param,param.gH,0.5*param.ppos_h)

	ds       = fN
	dn       = fNeg + fPosL + fPosH

	# Polymorphism
	neut::Array{Float64,1} = sfs_neut(param,binom)
	selH::Array{Float64,1} = if isinf(exp(param.gH * 2))
		sfs_pos_float(param,param.gH,param.ppos_h,binom)
	else
		sfs_pos(param,param.gH,param.ppos_h,binom)
	end
	selL::Array{Float64,1} = sfs_pos(param,param.gL,param.ppos_l,binom)
	selN::Array{Float64,1} = sfs_neg(param,param.ppos_h+param.ppos_l,binom)
	# Cumulative rates
	if(sum(param.cutoff)!=1)
		freq = round.(collect(1:(param.nn-1))./param.nn,digits=4)

		tmp = hcat(freq,neut,selH,selL,selN)
		tmp = @view tmp[(tmp[:,1] .>= param.cutoff[1]) .& (tmp[:,1] .<= param.cutoff[end]),2:end]

		tmp = cumulative_sfs(tmp,false)
	else
		tmp = cumulative_sfs(hcat(neut,selH,selL,selN),false)
	end

	split_columns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));
	neut, selH, selL, selN = split_columns(tmp)
	sel = (selH+selL)+selN

	##########
	# Output #
	##########
	analytical_m::Matrix{Float64} = vcat(param.B,param.al_low,param.al_tot,param.gam_neg,param.gL,param.gH,param.al,param.θ_coding)'
	analytical_ps::Matrix{Float64} = neut[param.dac]'
	analytical_pn::Matrix{Float64} = sel[param.dac]'
	analytical_f::Matrix{Float64} = hcat(ds,dn,fPosL,fPosH)
	# analytical_values::Array{Float64,2} = vcat(param.B,param.al_low,param.al_tot,param.gam_neg,param.gL,param.gH,param.al,param.θ_coding,neut[param.dac],sel[param.dac],ds,dn,fPosL,fPosH)'

	return (analytical_m,analytical_ps,analytical_pn,analytical_f)
end

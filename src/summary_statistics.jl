################################
###    Summary statistics    ###
################################
"""
	poisson_fixation(observed_values,λds, λdn)

Divergence sampling from Poisson distribution. The expected neutral and selected fixations are subset through their relative expected rates ([`fix_neut`](@ref), [`fix_neg_b`](@ref), [`fix_pos_sim`](@ref)). Empirical values are used are used to simulate the locus *L* along a branch of time *T* from which the expected *Ds* and *Dn* raw count estimated given the mutation rate (``\\mu``). Random number generation is used to subset samples arbitrarily given the success rate ``\\lambda`` in the distribution.

```math
\\mathbb{E}[D_N] = X \\in Poisson\\left(\\lambda = D \\times \\left[\\frac{\\mathbb{E}[D_+] + \\mathbb{E}[D_-]}{\\mathbb{E}[D_+] + \\mathbb{E}[D_-] + \\mathbb{E}[D_0]}\\right]\\right)
```
```math
\\mathbb{E}[D_S] = X \\in Poisson\\left(\\lambda = D \\times \\left[\\frac{\\mathbb{E}[D_0]}{\\mathbb{E}[D_+] + \\mathbb{E}[D_-] + \\mathbb{E}[D_0]}\\right]\\right)
```
# Arguments
 - `observed_values::Array`: Array containing the total observed divergence.
 - ` λds::Float64`: expected neutral fixations rate.
 - ` λdn::Float64`: expected selected fixations rate.
# Returns
 - `Array{Int64,1}` containing the expected count of neutral and selected fixations.

"""
function poisson_fixation(;observed_values::Array, λds::Array, λdn::Array,λweak::Array,λstrong::Array)

	ds = @. λds / (λds + λdn) * observed_values
	dn = @. λdn / (λds + λdn) * observed_values
	dweak = @. λweak / (λds + λdn) * observed_values
	dstrong = @. λstrong / (λds + λdn) * observed_values

	sampledDs     = pois_rand.(ds)
	sampledDn     = pois_rand.(dn)
	sampledWeak   = pois_rand.(dweak)
	sampledStrong = pois_rand.(dstrong)

	alphas = @. [sampledWeak/sampledDn sampledStrong/sampledDn (sampledWeak+sampledStrong)/sampledDn]

	out = alphas,sampledDn, sampledDs
	return out
end

"""
	poisson_polymorphism(observed_values,λps,λpn)

Polymorphism sampling from Poisson distributions. The total expected neutral and selected polimorphism are subset through the relative expected rates at the frequency spectrum ([`fix_neut`](@ref), [`sfs_neut`](@ref),). Empirical sfs are used to simulate the locus *L* along a branch of time *T* from which the expected *Ps* and *Pn* raw count are estimated given the mutation rate (``\\mu``). Random number generation is used to subset samples arbitrarily from the whole sfs given each frequency success rate ``\\lambda`` in the distribution.

The success rate managing the Poisson distribution by the observed count each frequency.  We considered both sampling variance and process variance is affecting the number of variable alleles we sample from SFS. This variance arises from the random mutation-fixation process along the branch. To incorporate this variance we do one sample per frequency-bin and use the total sampled variation and the SFS for the summary statistics.

```math
\\mathbb{E}[P_N] = \\sum_{x=0}^{x=1} X \\in Poisson\\left(\\lambda = SFS_{(x)} \\times \\left[\\frac{\\mathbb{E}[P_{+(x)}] + \\mathbb{E}[P_{-(x)}]}{\\mathbb{E}[P_{+(x)}] + \\mathbb{E}[P_{-(x)}] + \\mathbb{E}[P_{0(x)}]}\\right]\\right)
```

```math
\\mathbb{E}[P_S] = \\sum_{x=0}^{x=1} X \\in Poisson\\left(\\lambda = SFS_{(x)} \\times \\left[\\frac{\\mathbb{E}[P_{0(x)}]}{\\mathbb{E}[P_{+(x)}] + \\mathbb{E}[P_{-(x)}] + \\mathbb{E}[P_{0(x)}]}\\right]\\right)
```

# Arguments
 - `observed_values::Array{Int64,1}`: Array containing the total observed divergence.
 - ` λps::Array{Float64,1} `: expected neutral site frequency spectrum rate.
 - ` λpn::Array{Float64,1} `: expected selected site frequency spectrum rate.
# Returns
 - `Array{Int64,2}` containing the expected total count of neutral and selected polymorphism.

"""
function poisson_polymorphism(;observed_values::Array, λps::Array, λpn::Array)

	# Neutral λ;
	λ1 = @. λps / (λps + λpn) * observed_values
	# Selected λ;
	λ2 = @. λpn / (λps + λpn) * observed_values

	# Replace negative and 0 rates values. Strange behaviour depending on θ and Γ values
	# Relative rates output NaN values due to 0 divisons.

	replace!(λ1,NaN=>1)	
	replace!(λ2,NaN=>1)

	sampledPs = pois_rand.(λ1)
	sampledPn = pois_rand.(λ2)

	return (sampledPn, sampledPs)
end

"""
	sampled_alpha(observed_values,λds, λdn)

Ouput the expected values from the Poisson sampling process. Please check [`poisson_fixation`](@ref) and [`poisson_polymorphism`](@ref) to understand the samplingn process. α(x) is estimated through the expected values of Dn, Ds, Pn and Ps.

# Arguments
 - `param::parameters`: Array containing the total observed divergence.
 - `d::Array`: observed divergence.
 - `afs::Array`: observed polymorphism.
 - ` λdiv::Array{Float64,2}`: expected fixations rate.
 - ` λdiv::Array{Float64,2}`: expected site frequency spectrum rates.
# Returns
α_summaries,exp_dn,exp_ds,exp_pn,exp_ps,ssAlpha
 - `Array{Int64,2}` containing α(x) values.
 - `Array{Int64,1}` expected non-synonymous divergence.
 - `Array{Int64,1}` expected synonymous divergence.
 - `Array{Int64,1}` expected non-synonymous polymorphism.
 - `Array{Int64,1}` expected synonymous polymorphism.
 - `Array{Int64,1}` expected synonymous polymorphism.
 - `Array{Int64,1}` expected synonymous polymorphism.
 - `Array{Int64,2}` containing α(x) binned values.

	sampled_from_rates(gammaL,gammaH,ppos_l,ppos_h,observedData,nopos)

"""
function sampled_from_rates(models::Matrix{Float64},fs::Vector{Float64},d::Vector{Int64},neut::Matrix{Float64},sel::Matrix{Float64},dsdn::Matrix{Float64})

	ds             = dsdn[:,1]
	dn             = dsdn[:,2]
	dweak          = dsdn[:,3]
	dstrong        = dsdn[:,4]
	gn             = abs.(models[:,4])
	sh             = round.(models[:,end-1],digits=5)

	## Outputs
	alphas, exp_dn, exp_ds = poisson_fixation(observed_values=d,λds=ds,λdn=dn,λweak=dweak,λstrong=dstrong);

	exp_pn, exp_ps         = poisson_polymorphism(observed_values=fs,λps=permutedims(neut),λpn=permutedims(sel));

	## Alpha from expected values. Used as summary statistics
	α_summaries = @. round(1 - ((exp_ds/exp_dn) * (exp_pn/exp_ps)'),digits=5);


	expected_values = hcat(round.(alphas,digits=5),gn,sh,α_summaries);

	return expected_values
end

"""
	summary_statistics(param::parameters,rates::JLD2.JLDFile,analysis_folder::String,summstat_size::Int64,replicas::Int64,bootstrap::Bool)

Estimate summary statistics using observed data and analytical rates. *analysis_folder* will check for the SFS and divergence file and will be used to output summary statistics

# Arguments
 - `param::parameters` : Mutable structure containing the models
 - `h5_file::String` : HDF5 containing solved models, fixation and polymorphic rates
 - `analysis_folder::String` : Folder containing the SFS and divergence files. It will be used to output the observed data and summary estatistics.
 - `summstat_size::Int64` : Number of summary statistics
 - `replicas::Int64` : Number of bootstrap replicas
 - `bootstrap::Bool` : Boolean to perform or not bootstrapping
# Output
 - Summary statistics to ABC inference
"""
function summary_statistics(;param::parameters,h5_file::String,analysis_folder::String,summstat_size::Int64,replicas::Int64=1,bootstrap::Bool=false)

	#Opening files
	s_file   = filter(x -> occursin("sfs",x), readdir(analysis_folder,join=true));
	d_file   = filter(x -> occursin("div",x), readdir(analysis_folder,join=true));

	sfs,divergence,α = open_sfs_div(s_file,d_file,param.dac,replicas,bootstrap);

	#Open rates
	h = jldopen(h5_file);
	tmp    = h[string(param.N) * "/" * string(param.n)]

	#Subset index
	idx    = sample(1:size(tmp["models"],1),summstat_size,replace=false)
	models = Array(view(tmp["models"],idx,:));
	dsdn   = Array(view(tmp["dsdn"],idx,:));

	# Filtering polymorphic rate by dac
	n    = hcat(map(x -> view(tmp["neut"][x],:),param.dac)...);
	s    = hcat(map(x -> view(tmp["sel"][x],:),param.dac)...);
	neut = Array(view(n,idx,:));
	sel  = Array(view(s,idx,:));

	#Making summaries
	expected_values = sampled_from_rates(models,sfs[1],divergence[1],neut,sel,dsdn);

	w(x,name) = CSV.write(name,DataFrame(x,:auto),delim='\t',header=false);

	# Controling outlier cases
	flt_inf(e) = replace!(e, -Inf=>NaN)
	flt_nan(e) = e[vec(.!any(isnan.(e),dims=2)),:]
	
	expected_values = flt_inf(expected_values)
	expected_values = flt_nan(expected_values)
	expected_values = expected_values[(expected_values[:,3] .> 0 ) .& (expected_values[:,3] .<1 ),:]
	
	# Writting ABCreg input
	w(vcat(α...), analysis_folder * "/alphas.txt");
	w(expected_values, analysis_folder * "/summstat.txt");

	# pmapbatch(w, α, analysis_folder * "/alphas_" .* string.(1:size(sfs,1)) .* ".txt");
	# pmapbatch(w, expected_values,  analysis_folder * "/summstat_" .* string.(1:size(sfs,1)) .* ".txt");

	return(expected_values)
end

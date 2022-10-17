########################
# Parse data functions # 
########################

"""
	parse_sfs(parameters;sample_size::Int64,data::String)

Function to parse polymorphism and divergence by subset of genes. The input data is based on supplementary material described at [Uricchio et al. 2019](https://doi.org/10.1038/s41559-019-0890-6). Please be sure the file is tabulated.

| GeneId | Pn | DAF seppareted by commas | Ps | DAF separated by commas | Dn | Ds |
|--------|----|--------------------------|----|-------------------------|----|----|
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |

# Arguments
 - `sample_size::Int64`: data sample size.
 - `data`: File containing polymorphic and divergence data
 - `gene_list`: File containing gene IDs to subset. You can perform multiple subset input a file where eachrow contain the gene IDs to subset.
 - `sfs_columns::Array{Int64,1}`: non-synonymous and synonymous daf columns. Please introduce first the non-synonymous number.
 - `div_columns::Array{Int64,1}`: non-synonymous and synonymous divergence columns. Please introduce first the non-synonymous number.
 - `bins::Array{Int64,1}`: bin to collapse the SFS from `sample_size` to new `bin` size

# Returns
 - `Array{Float64,1}`: α values
 - `Array{Float64,2}`: Site Frequency Spectrum
 - `Array{Float64,1}`: Synonymous and non-synonymous divergence counts
"""
function parse_sfs(
	param::parameters;
	data::S,
	gene_list::Union{Nothing,S} = nothing,
	sfs_columns::Vector{Int64} = [3, 5],
	div_columns::Vector{Int64} = [6, 7],
	m_columns::Vector{Int64} = [8, 9],
	bins::Union{Nothing,Int64} = nothing,
	isolines::Bool = false) where {S<:AbstractString}

	@unpack n, cutoff = param

	if isolines
		s_size = n
	else
		s_size = (n * 2)
	end

	df = CSV.read(data, header = false, delim = '\t', DataFrame)

	if (!isnothing(gene_list))
		ids = @view df[:, 1]

		gene_matrix = Array(CSV.read(gene_list, header = false, DataFrame))
		m, n = size(gene_matrix)

		if n == 1
			gene_matrix = permutedims(gene_matrix)
		end

		out = SubDataFrame[]
		for c ∈ eachrow(gene_matrix)
			tmp = @view df[filter(!isnothing, indexin(c, ids)), :]
			push!(out, tmp)
		end

		α, sfs, divergence, m = unzip(map(i -> get_pol_div(i,s_size,cutoff,sfs_columns,div_columns,m_columns,bins),out))
	else
		α, sfs, divergence, m = get_pol_div(df, s_size, cutoff, sfs_columns, div_columns, m_columns, bins)
		α                     = [α]
		sfs                   = [sfs]
		divergence            = [divergence]
		m                     = [m]
	end

	return α, sfs, divergence, m
end

function get_pol_div(
	df_subset::Union{DataFrame,SubDataFrame},
	s_size::Int64,
	cutoff::Vector{Float64},
	sfs_columns::Vector{Int64},
	div_columns::Vector{Int64},
	m_columns::Vector{Int64},
	bins::Union{Nothing,Int64})

	g(x)   = parse.(Float64,filter(y -> y!="", x))

	freq   = OrderedDict(round.(collect(1:(s_size-1)) / s_size, digits = 4) .=> 0)
	tmp    = split.(df_subset[:, sfs_columns], ",")

	pn     = vcat(g.(tmp[:, 1])...)
	ps     = vcat(g.(tmp[:, 2])...)

	# Round SFS frequencies to the lowest floating value of the dataset independtly of the sample size. Needed to countmap and merge.
    dgts = length(string(minimum(vcat(pn,ps))))-2

    pn   = sort(OrderedDict(countmap(round.(pn, digits = dgts))))
    ps   = sort(OrderedDict(countmap(round.(ps, digits = dgts))))

	# Dn, Ds, Pn, Ps, sfs
	Dn     = sum(df_subset[:, div_columns[1]])
	Ds     = sum(df_subset[:, div_columns[2]])

	sfs_pn = reduce(vcat, values(merge(+, freq, pn)))
	sfs_ps = reduce(vcat, values(merge(+, freq, ps)))

	if (!isnothing(bins))
		sfs_pn = reduce_sfs(hcat(collect(1:(s_size-1)), sfs_pn), bins)[:, 2]
		sfs_ps = reduce_sfs(hcat(collect(1:(s_size-1)), sfs_ps), bins)[:, 2]
		sfs    = hcat(collect(1:(bins-1)) ./ bins, sfs_pn, sfs_ps, 1:(bins-1))
	else
		sfs    = hcat(freq.keys, sfs_pn, sfs_ps, (1:s_size-1))
	end

	# Filtering SFS and changing frequency to DAC
	sfs   = sfs[sfs[:, 1].>=cutoff[1] .&& sfs[:, 1].<=cutoff[2], [4, 2, 3]]

	scumu = cumulative_sfs(sfs)

	α     = round.(1 .- (Ds / Dn .* scumu[:, 2] ./ scumu[:, 3]), digits = 5)

	m::Matrix{Int64} = try
		[sum(df_subset[:, m_columns[1]]) sum(df_subset[:, m_columns[2]])]
	catch
		m = [0 0]
	end

	return (α, sfs, [Dn Ds], m)
end

function data_to_poisson(
	sfs::Vector{Matrix{Float64}},
	divergence::Vector{Matrix{Int64}},
	dac::Vector{Int64})


	f(x::Matrix{Float64},d::Vector{Int64}=dac) = map(z -> sum(x[x[:, 1].==z, 2:3]), d)
	al(a, b)  = hcat(a[:, 1], @. round(1 - (b[2] / b[1] * a[:, 2] / a[:, 3]), digits = 5))

	scumu     = cumulative_sfs.(sfs)
	s_poisson = f.(scumu)
	d_poisson = [[sum(divergence[i][1:2])] for i ∈ eachindex(divergence)]
	α_x       = al.(scumu, divergence)
	α_observed        = map(x -> permutedims(x[in(dac).(x[:, 1]), 2]), α_x)

	return(α_observed, s_poisson, d_poisson)
end

#################
# ABC functions #
#################

q_lower(x::Vector{Float64}) = quantile(x, [0.5])
q_upper(x::Vector{Float64}) = quantile(x, [0.95])
write_files(x::Matrix{Float64}, name::String) = CSV.write(name, Tables.table(x), delim = '\t', header = false);

"""
	ABCreg(analysis_folder, S, P, tol, abcreg)

Performing ABC inference using ABCreg. Please, be sure your analysis_folder contain the files alphas.txt and summaries.txt produced by Analytical.summaryStatsFromRates()
# Arguments
 - `analysis_folder::String` : Folder containing the observed data and summary estatistics. It will be used to output the posterior distributions
 - `P::Int64` : Number of parameters to perform the inference.
 - `S::Int64` : Number of summary stastitics to perform the inference.
 - `tol::Float64` : Tolerance value. It define the number of accepted values at ABC inference
 - `abcreg::String` : Path to ABCreg binary
# Output
Files containing posterior distributions from ABCreg
"""
function ABCreg(;
	analysis_folder::String,
	S::Int64,
	P::Int64 = 5,
	tol::Float64,
	abcreg::String)

	# List alphas and summstat files
	a_file   = filter(x -> occursin("alphas", x), readdir(analysis_folder, join = true))
	sum_file = filter(x -> occursin("summstat", x), readdir(analysis_folder, join = true))

	# Creating output names
	out = analysis_folder .* "/out_" .* string.(1:size(a_file, 1))

	@info "Running ABCreg"
	# Using mapi to limit tasks. ThreadsX.map Not working properly with bash function. Try change to pipeline or something else.
	r(a::String, s::String, o::String, abcreg::String=abcreg,P::Int64=P,S::Int64=S,tol::Float64=tol) = run(`$abcreg -d $a -p $s -P $P -S $S -t $tol -b $o`);

	# Using mapi instead map. Bash interaction not working as expected
	ThreadsX.mapi(
		(x, y, z) -> r(x, y, z),
		a_file,
		sum_file,
		out,
		ntasks = Threads.nthreads()
	)

	@info "Opening and filtering posteriors distributions"
	out        = filter(x -> occursin("post", x), readdir(analysis_folder, join = true))
	out        = filter(x -> !occursin(".1.", x), out)
	
	# Control outlier inference. 2Nes non negative values
	open(x)    = Array(CSV.read(x, DataFrame))
	flt(x)     = x[(x[:, 4].>0).&(x[:, 1].>0).&(x[:, 2].>0).&(x[:, 3].>0), :]
	posteriors = flt.(open.(out))

	return posteriors
end

"""
	get_mode(x)
Estimating mode from posterior following R package abc

# Arguments
 - `posterior::Matrix{Float64}` : posterior distribution.
# Output
 - `{Float64}` : mode estimation by column.
Files containing posterior distributions from ABCreg	
"""
function get_mode(posterior::Matrix{Float64})
	
	# Allocating outputs
	out = zeros((1, size(posterior, 2)))
	
	for j::Int64=1:size(posterior, 2)
		y      = kde(posterior[:, j])
		m      = collect(y.x)[y.density.==maximum(y.density)][1]
		out[j] = m
	end

	return (out)
end

function summary_abc(
	posteriors::Vector{Matrix{Float64}};
	stat::String="Mean",
	plot_path::Union{String,Nothing}=nothing)

	@info "Computing statistics over posteriors distributions"

	p_min    = vcat(map(x -> minimum(x, dims = 1), posteriors)...)
	p_max    = vcat(map(x -> maximum(x, dims = 1), posteriors)...)
	p_mean   = vcat(map(x -> mean(x, dims = 1), posteriors)...)
	p_median = vcat(map(x -> median(x, dims = 1), posteriors)...)
	p_mode   = vcat(map(get_mode, posteriors)...)
	p_lower  = vcat(map(x -> mapslices(q_lower, x, dims = 1), posteriors)...)
	p_upper  = vcat(map(x -> mapslices(q_upper, x, dims = 1), posteriors)...)


	if length(posteriors) == 1

		stats = DataFrame(
			:Stats => repeat(
				["Min.", "0.5% Perc.", "Median", "Mean", "Mode", "99.5% Perc.", "Max"],
				inner = length(posteriors)
			)
		)

		tmp = DataFrame(
			vcat(p_min, p_lower, p_median, p_mean, p_mode, p_upper, p_max),
			[:α_weak, :α_strong, :α, :γ, :β]
		)

		df = hcat(stats, tmp)

		@info df

		if !isnothing(plot_path)
			plot_posterior!(posteriors[1], plot_path)
		end

		return df, df[df.Stats.==uppercasefirst(stat), 2:end]
	else

		df = DataFrame[]
		for i::Int64 = 1:length(posteriors)
			stats = DataFrame(
				:Stats => [
					"Min.",
					"0.5% Perc.",
					"Median",
					"Mean",
					"Mode",
					"95% Perc.",
					"Max",
				]
			)

			tmp = DataFrame(
				hcat(
					p_min[i, :],
					p_lower[i, :],
					p_median[i, :],
					p_mean[i, :],
					p_mode[i, :],
					p_upper[i, :],
					p_max[i, :],
				)',
				[:α_weak, :α_strong, :α, :γ, :β]
			)
			push!(df, hcat(stats, tmp))
		end

		df_vcat = vcat(df...)

		return df, df_vcat[df_vcat.Stats .== uppercasefirst(stat), 2:end]
	end
end




##################
# Dofe functions #
##################

"""
	Bootstrap data following polyDFE manual
"""
function bootstrap_data!(
	sfs::Matrix{Float64},
	divergence::Matrix{Float64},
	bootstrap::Int64)

	pr(x::Matrix{Float64}) = hcat(x[:, 1], pois_rand.(x[:, 2:end]))

	sfs         = repeat(sfs, bootstrap)
	divergence  = repeat(divergence, bootstrap)
	sfs[2:end] .= pr.(sfs[2:end])

	return sfs, divergence
end

function grapes(
	sfs::Vector{Matrix{Float64}},
	divergence::Vector{Matrix{Int64}},
	m::Vector{Matrix{Int64}},
	model::String,
	folder::String,
	bins::Int64)

	if (isnothing(CondaPkg.which("grapes")))
		CondaPkg.add("grapes-static", channel = "genomedk")
	end

	grapes = CondaPkg.which("grapes")

	sfs    = reduce_sfs.(sfs, bins)

	pn     = map(x -> permutedims(x[:, 2]), sfs)
	ps     = map(x -> permutedims(x[:, 3]), sfs)

	dn     = map(x -> x[1], divergence)
	ds     = map(x -> x[2], divergence)

	mn     = map(x -> x[1], m)
	ms     = map(x -> x[2], m)

	idx    = string.(collect(1:length(sfs)))

	# Temporal function to broadcast pol and div data
	f(pn, ps, dn, ds, mn, ms, w) =
		DataFrame(hcat("dofe_" * string(w), bins, mn, pn, ms, ps, mn, dn, ms, ds...), :auto)

	dofe          = @. f(pn, ps, dn, ds, mn, ms, idx)
	h             = fill(DataFrame(["" ""; "#unfolded" ""], :auto), length(sfs))
	output_dofe   = @. folder * "/dofe_" * idx * ".txt"
	output_grapes = @. folder * "/dofe_" * idx * "." * model


	@. write_files(h, output_dofe)
	@. write_files(dofe, output_dofe, fill(true, length(sfs)))

	r(d, o, m = model, gr = grapes) = run(`$gr -in $d -out $o -model $m`)

	@suppress_out begin
		progress_pmap(r, output_dofe, output_grapes)
	end

	df = @suppress begin
		CSV.read.(output_grapes, DataFrame, footerskip = 1, skipto = 3)
	end

	return (vcat(df...))
end

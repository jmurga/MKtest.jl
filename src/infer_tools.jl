"""
	parse_sfs(;sample_size::Int64,data::String)

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
function parse_sfs(;param::parameters,data::S,gene_list::Union{Nothing,S}=nothing,sfs_columns::Vector{Int64}=[3,5],div_columns::Vector{Int64}=[6,7],m_columns::Vector{Int64}=[8,9],bins::Union{Nothing,Int64}=nothing,isolines::Bool=false) where S <: AbstractString

	@unpack n, cutoff = param

	if isolines
		s_size = n
	else
		s_size = (n*2)
	end
	
	df = CSV.read(data,header=false,delim='\t',DataFrame)	

	if(!isnothing(gene_list))
		ids = @view df[:,1];

		gene_matrix = Array(CSV.read(gene_list,header=false,DataFrame))
		m,n = size(gene_matrix)

		if n == 1
			gene_matrix = permutedims(gene_matrix)
		end
		
		out = SubDataFrame[]
		for c in eachrow(gene_matrix)
			tmp = @view df[filter(!isnothing,indexin(c,ids)),:];
			push!(out,tmp);
		end

		α, sfs, divergence, m = unzip(map(i-> get_pol_div(i,s_size,cutoff,sfs_columns,div_columns,m_columns,bins),out));

	else
		α, sfs, divergence,m = get_pol_div(df,s_size,cutoff,sfs_columns,div_columns,m_columns,bins);
		α = [α]; sfs = [sfs]; divergence = [divergence]; m = [m]
	end

	return α, sfs, divergence, m
end

function get_pol_div(df_subset::Union{DataFrame,SubDataFrame},s_size::Int64,cutoff::Vector{Float64},sfs_columns::Vector{Int64},div_columns::Vector{Int64},m_columns::Vector{Int64},bins::Union{Nothing,Int64})
	
	g(x) = parse.(Float64,x[2:end-1])
	
	freq = OrderedDict(round.(collect(1:(s_size-1))/s_size,digits=4) .=> 0)
	# freq = OrderedDict(collect(1:(s_size-1)) .=> 0)

	tmp  = split.(df_subset[:,sfs_columns], ",")

	pn   = sort!(OrderedDict(countmap(round.(reduce(vcat,tmp[:,1] .|> g),digits=4))))
	ps   = sort!(OrderedDict(countmap(round.(reduce(vcat,tmp[:,2] .|> g),digits=4))))

	# Dn, Ds, Pn, Ps, sfs
    Dn     = sum(df_subset[:,div_columns[1]])
    Ds     = sum(df_subset[:,div_columns[2]])
    Pn     = sum(values(pn))
    Ps     = sum(values(ps))
    sfs_pn = reduce(vcat,values(merge(+,freq,pn)))
    sfs_ps = reduce(vcat,values(merge(+,freq,ps)))

	if(!isnothing(bins))
		sfs_pn = reduce_sfs(hcat(collect(1:(s_size-1)),sfs_pn),bins)[:,2]
		sfs_ps = reduce_sfs(hcat(collect(1:(s_size-1)),sfs_ps),bins)[:,2]

		sfs   = hcat(collect(1:(bins-1))./bins,sfs_pn,sfs_ps,1:(bins-1))
	else
		sfs   = hcat(freq.keys,sfs_pn,sfs_ps,(1:s_size-1))
	end
	
	# Filtering SFS and changing frequency to DAC
	sfs = sfs[sfs[:,1] .>= cutoff[1].&& sfs[:,1] .<= cutoff[2],[4,2,3]]
 
	scumu = cumulative_sfs(sfs)
	α     = round.(1 .- (Ds/Dn .*  scumu[:,2] ./scumu[:,3]),digits=5)

	m = try
		[sum(df_subset[:,m_columns[1]]) sum(df_subset[:,m_columns[2]])]
	catch
		m = [0 0]
	end

	return (α,sfs,[Dn Ds],m)
end

"""
	ABCreg(analysis_folder, replicas, S, tol, abcreg)

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
function ABCreg(;analysis_folder::String,S::Int64,P::Int64=5,tol::Float64,abcreg::String)
	
	# List alphas and summstat files
	a_file     = filter(x -> occursin("alphas",x), readdir(analysis_folder,join=true));
	sum_file   = filter(x -> occursin("summstat",x), readdir(analysis_folder,join=true));

	# Creating output names
	# out = [analysis_folder .* "/out"]
	out = analysis_folder .* "/out_" .* string.(1:size(a_file,1))

	r(a,s,o,abcreg=abcreg,P=P,S=S,tol=tol) = run(`$abcreg -d $a -p $s -P $P -S $S -t $tol -b $o`)
	
	progress_pmap(r,a_file,sum_file,out);

end

"""
	Bootstrap data following polyDFE manual
"""
function data_to_poisson(sfs::Vector{Matrix{Float64}},divergence::Vector{Matrix{Int64}},dac::Vector{Int64},bootstrap::Union{Bool,Int64})
	
	if bootstrap != false
		sfs         = repeat(sfs,bootstrap)
		divergence  = repeat(divergence,bootstrap)
		pr(x)       = hcat(x[:,1],pois_rand.(x[:,2:end]))
		sfs[2:end] .= pr.(sfs[2:end])
	end

    scumu      = cumulative_sfs.(sfs)

    f(x,d=dac) = map(z->sum(x[x[:,1].==z,2:3]),d)
    al(a,b)    = hcat(a[:,1],@. round(1 - (b[2]/b[1] * a[:,2]/a[:,3]),digits=5))

    s          = f.(scumu)
    d          = [[sum(divergence[i][1:2])] for i in eachindex(divergence)]

    α_x        = al.(scumu,divergence)
    α          = map(x -> permutedims(x[in(dac).(x[:,1]),2]),α_x)

	return(α,s,d)
end

function open_sfs_div(x::Array{String,1},y::Array{String,1},dac::Vector{Int64},bootstrap::Union{Bool,Int64})

	sfs = Array.(CSV.read.(x,DataFrame,header=false))
	divergence = Array.(CSV.read.(y,DataFrame,header=false))

	if bootstrap != false
		sfs = repeat(sfs,bootstrap)
		divergence = repeat(divergence,bootstrap)
		pr(x) = hcat(x[:,1],pois_rand.(x[:,2:end]))
		sfs[2:end] .= pr.(sfs[2:end])
	end

	scumu = cumulative_sfs.(sfs)
	f(x,d=dac) = sum(x[:,2:3],dims=2)[d]
	s = f.(scumu)

	d = [[sum(divergence[i][1:2])] for i in eachindex(divergence)]
	al(a,b,c=dac) = @. round(1 - (b[2]/b[1] * a[:,2]/a[:,3])[c],digits=5)
	α = permutedims.(al.(scumu,divergence))
	return(s,d,α)
end

"""
	Function to download and source plotMap function. We do not include at part of the module to avoid external dependecies. Once the function is execute properly you will have a function called *plotMap which used R to 
		estimate and plot the Maximum A Posterior following ABCreg example. It uses locfit and ggplot2 libraries.
#Arguments
- `script_path::String`: Path and name to save MAP scritp.
"""
function source_plot_map_r(script_path::String)

	try
		download("https://raw.githubusercontent.com/jmurga/MKtest.jl/main/scripts/plot_map.jl",script_path)
		include(script_path)
	catch
		println("\nPlease be sure you have R and the pacakges data.table, ggplot2 and locfit installed in your system. You can install them using Conda inside julia with the following commands:\n\n\tusing Pkg\n\tPkg.add(\"conda\")\n\tENV[\"R_HOME\"]=\"*\"\n\tPkg.add(\"Conda\")\n\tusing Conda\n\tConda.add(\"r-base\",channel=\"conda-forge\")\n\tConda.add([\"r-locfit\",\"r-ggplot2\",\"r-data.table\",\"r-r.utils\"],channel=\"conda-forge\")\n\tPkg.add(\"RCall\")\n")
	end
end


function grapes(sfs::Vector{Matrix{Float64}},divergence::Vector{Matrix{Int64}},m::Vector{Matrix{Int64}},model::String,output::String,grapes::String) 
	
    sfs = reduce_sfs.(sfs,20);
	
    pn  = map(x-> permutedims(x[:,2]),sfs);
    ps  = map(x-> permutedims(x[:,3]),sfs);
	
    dn  = map(x->x[1],divergence);
    ds  = map(x->x[2],divergence);

    mn  = map(x->x[1],m);
    ms  = map(x->x[2],m);

    idx = string.(collect(1:length(sfs)));

	f(pn,ps,dn,ds,mn,ms,w) = DataFrame(hcat("dofe_"*string(w),20,mn,pn,ms,ps,mn,dn,ms,ds...),:auto)

    dofe          = f.(pn,ps,dn,ds,mn,ms,idx);
    h            = fill(DataFrame(["" ""; "#unfolded" ""],:auto),length(sfs))
    output_dofe  = output .* "_" .* idx .* ".txt"
    output_grapes = output .* "_" .* idx .* "." .* model

	w(x,name,a=false) = CSV.write(name,x,delim='\t',header=false,append=a);

	w.(h,output_dofe)
	w.(dofe,output_dofe,fill(true,length(sfs)))

	r(d,o,model=model,grapes=grapes) = run(`$grapes -in $d -out $o -model $model`)
	
	@suppress_out begin
		progress_pmap(r,output_dofe,output_grapes);   
	end

	df = @suppress begin
		CSV.read.(output_grapes,DataFrame,footerskip=1,skipto=3);
	end

	return(vcat(df...))
end

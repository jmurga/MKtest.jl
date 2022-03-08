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
 - `data`: String or Array of strings containing files names with full path.
 - `gene_list`: String or Array of strings containing files names to subset data.
 - `sfs_columns::Array{Int64,1}`: non-synonymous and synonymous daf columns. Please introduce first the non-synonymous number.
 - `div_columns::Array{Int64,1}`: non-synonymous and synonymous divergence columns. Please introduce first the non-synonymous number.
 - `bins::Array{Int64,1}`: bin to collapse the SFS from `sample_size` to new `bin` size

# Returns
 - `Array{Float64,1}`: α values
 - `Array{Float64,2}`: Site Frequency Spectrum
 - `Array{Float64,1}`: Synonymous and non-synonymous divergence counts
"""
function parse_sfs(;sample_size::Int64,data::Union{String,DataFrame},gene_list::Union{Nothing,Vector{String},Matrix{String}}=nothing,sfs_columns::Array{Int64,1}=[3,5],div_columns::Array{Int64,1}=[6,7],bins::Union{Nothing,Int64}=nothing,isolines::Bool=false)

	g(x) = parse.(Float64,x[2:end-1])
	
	if isolines
		s = sample_size
	else
		s = (sample_size*2)
	end
	
	freq = OrderedDict(round.(collect(1:(s-1))/s,digits=4) .=> 0)

	if typeof(data) == String
		df   = CSV.read(data,header=false,delim='\t',DataFrame)
	else
		df = data
	end

	if(!isnothing(gene_list))
		df =  vcat([ df[df[:,1] .==i,:]  for i in gene_list]...);
	end

	#=if(!isnothing(B))
		df = df[df[:,end] .== B,:]
		println(nrow(df))
		tmp  = split.(df[:,sfs_columns], ",")
	else
	end=#
	
	tmp  = split.(df[:,sfs_columns], ",")

	pn   = sort!(OrderedDict(round.(reduce(vcat,tmp[:,1] .|> g),digits=4) |> countmap))
	ps   = sort!(OrderedDict(round.(reduce(vcat,tmp[:,2] .|> g),digits=4) |> countmap))

	# Dn, Ds, Pn, Ps, sfs
	Dn           = sum(df[:,div_columns[1]])
	Ds           = sum(df[:,div_columns[2]])
	Pn           = sum(values(pn))
	Ps           = sum(values(ps))
	sfs_pn        = reduce(vcat,values(merge(+,freq,pn)))
	sfs_ps        = reduce(vcat,values(merge(+,freq,ps)))

	if(!isnothing(bins))
        sfs_pn = reduce_sfs(hcat(collect(1:(s-1)),sfs_pn),bins)[:,2]
        sfs_ps = reduce_sfs(hcat(collect(1:(s-1)),sfs_ps),bins)[:,2]

        sfs   = reduce_sfs(hcat(freq.keys,merge(+,freq,pn).vals,merge(+,freq,ps).vals),bins)
	else
        sfs   = hcat(freq.keys,merge(+,freq,pn).vals,merge(+,freq,ps).vals)
	end

	scumu = cumulative_sfs(sfs)
    α     = round.(1 .- (Ds/Dn .*  scumu[:,2] ./scumu[:,3]),digits=5)

	return (α,sfs,[Dn,Ds])
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
    out = [analysis_folder .* "/out"]
    # out = analysis_folder .* "/out_" .* string.(1:size(a_file,1))

	r(a,s,o,abcreg=abcreg,P=P,S=S,tol=tol) = run(`$abcreg -d $a -p $s -P $P -S $S -t $tol -b $o`)

	r.(a_file,sum_file,out);
end

"""
	Bootstrap data following polyDFE manual
"""
function bootstrap_data(s_file::Array{Float64,2},d_file::Array{Float64,2},replicas::Int64,outputFolder::String)
	
	# Open Data
	sfs        = Array(CSV.read(s_file,DataFrame))
	divergence = fill(Array(CSV.read(d_file,DataFrame)),replicas)
	scumu      = fill(cumulative_sfs(sfs[:,2:end]),replicas)

	# Bootstraping
	b(x)       = pois_rand.(x)
	bootstrap  = b.(scumu)

	# Output
	outSfs = @. output * "sfs" * string(1:replicas) * ".tsv"
	outDiv = @. output * "div" * string(1:replicas) * ".tsv"
	for i  = 1:replicas
		tmp = hcat(sfs[:,1],bootstrap[i])
		CSV.write(out,DataFrame(tmp),header=false)
		CSV.write(out,DataFrame(divergence[i]),header=false)
	end
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

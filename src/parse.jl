########################
# Parse data functions # 
########################

"""
	parse_sfs(param;data,gene_list,sfs_columns,div_columns,l_columns,bins,isolines)

Function to parse polymorphism and divergence by subset of genes. The input data is based on supplementary material described at [Uricchio et al. 2019](https://doi.org/10.1038/s41559-019-0890-6). Please be sure the file is tabulated.

| GeneId | Pn | DAF seppareted by commas | Ps | DAF separated by commas | Dn | Ds |
|--------|----|--------------------------|----|-------------------------|----|----|
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |
| XXXX   | 0  | ,0,0,0,0,0,0,0,0         | 0  | ,0,0,0,0,0,0,0,0        | 0  | 0  |

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model.
 - `data::S`: path to polymorphic and divergence raw data.
 - `gene_list:S`: path to gene list. 
 - `sfs_columns::Vector{Int64}`: column numbers corresponding to Pn and Ps frequencies. Default values [3,5]
 - `div_columns::Vector{Int64}`: column numbers corresponding to Dn and Ds frequencies. Default values [6,7]
 - `l_columns::Vector{Int64}`: column numbers corresponding to total number of non-synonymous and synonymous sites. Not require at data. Default values [7,8]
 - `bins::Union{Nothing,Int64}`: size to bin the SFS independtly of the sample size.
# Returns
 - `Vector{Vector{Float64}}`: α_x values. Estimated using the cumulative SFS
 - `Vector{Matrix{Float64}}`: Site Frequency Spectrum (non-cumulative)
 - `Vector{Matrix{Int64}}`: Synonymous and non-synonymous divergence counts
"""
function parse_sfs(
    param::parameters;
    data::S,
    gene_list::Union{Nothing,Vector{S},S} = nothing,
    sfs_columns::Vector{Int64} = [3, 5],
    div_columns::Vector{Int64} = [6, 7],
    l_columns::Vector{Int64} = [8, 9],
    bins::Union{Nothing,Int64} = nothing,
) where {S<:AbstractString}
    @unpack n, cutoff, isolines = param

    if isolines
        s_size = n
    else
        s_size = (n * 2)
    end

    df = CSV.read(data, header = false, delim = '\t', DataFrame)

    if (!isnothing(gene_list))
        ids = @view df[:, 1]

        @assert gene_list isa AbstractString || gene_list isa Vector{String} "Please input a valid gene_list file or vector containing gene ids"

        if gene_list isa AbstractString
            gene_matrix = Array(CSV.read(gene_list, header = false, DataFrame))
        else
            gene_matrix = permutedims(gene_list)
        end

        m, n = size(gene_matrix)

        if n == 1
            out = [filter_df(df,ids,gene_matrix)]
        else
            out = Vector{SubDataFrame}(undef,m)
            Threads.@threads for i=1:m
                out[i] = filter_df(df,ids,view(gene_matrix,i,:))
            end
        end

        α, sfs, divergence = unzip(
            map(
                i -> get_pol_div(
                    i,
                    s_size,
                    cutoff,
                    sfs_columns,
                    div_columns,
                    l_columns,
                    bins,
                ),
                out,
            ),
        )
    else
        α, sfs, divergence =
            get_pol_div(df, s_size, cutoff, sfs_columns, div_columns, l_columns, bins)
        α = [α]
        sfs = [sfs]
        divergence = [divergence]
    end

    return α, sfs, divergence
end

# Not duplicated entries in ids or df, so get first index using indexin is valid
filter_df(df::DataFrame,ids::SubArray,row::Union{SubArray,Array}) = view(df,filter(!isnothing, indexin(row, ids)), :)


function get_pol_div(
    df_subset::Union{DataFrame,SubDataFrame},
    s_size::Int64,
    cutoff::Vector{Float64},
    sfs_columns::Vector{Int64},
    div_columns::Vector{Int64},
    l_columns::Vector{Int64},
    bins::Union{Nothing,Int64},
)::Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}}

    pn, ps  = eachcol(view(df_subset,:, sfs_columns))
    pn = split(join(pn),",");
    ps = split(join(ps),",");
    pn = parse.(Float64,view(pn,pn.!=""));
    ps = parse.(Float64,view(ps,ps.!=""));
    pn = view(pn,pn.!=0);
    ps = view(ps,ps.!=0);

    # Round SFS frequencies to the lowest floating value of the dataset independtly of the sample size. Needed to countmap and merge.

    freq = OrderedDict(round.(collect(1:(s_size-1)) / s_size, digits = 4) .=> 0)

    c_pn = sort(countmap(round.(pn, digits = 4)))
    c_ps = sort(countmap(round.(ps, digits = 4)))

    # Dn, Ds, Pn, Ps, sfs
    Dn = sum(view(df_subset,:, div_columns[1]))
    Ds = sum(view(df_subset,:, div_columns[2]))
    sfs_pn = merge(+, freq, c_pn).vals
    sfs_ps = merge(+, freq, c_ps).vals

    sfs = if (!isnothing(bins))
        sfs_pn_binned::Vector{Float64} = reduce_sfs(hcat(collect(1:(s_size-1)), sfs_pn), bins)[:, 2]
        sfs_ps_binned::Vector{Float64} = reduce_sfs(hcat(collect(1:(s_size-1)), sfs_ps), bins)[:, 2]
        hcat(collect(1:(bins-1)) ./ bins, sfs_pn, sfs_ps, 1:(bins-1))
    else
        hcat(freq.keys, sfs_pn, sfs_ps, (1:(s_size-1)))
    end;

    # Filtering SFS and changing frequency to DAC
    sfs_flt = sfs[view(sfs,:, 1).>=cutoff[1] .&& view(sfs,:, 1) .<=cutoff[2], [4, 2, 3]];

    scumu = cumulative_sfs(sfs_flt);

    α = round.(1 .- (Ds / Dn .* scumu[:, 2] ./ scumu[:, 3]), digits = 5);

    l = try
        [round(sum(df_subset[:, l_columns[1]]),digits=2) round(sum(df_subset[:, l_columns[2]]),digits=2)]
    catch
        [0 0]
    end;

    return (α, sfs_flt, hcat([Dn Ds], l))
end


function data_to_poisson(
    sfs::Vector{Matrix{Float64}},
    divergence::Vector{Matrix{Float64}},
    dac::Vector{Int64}
)
    scumu = cumulative_sfs.(sfs)
    s_poisson = fill(zeros(length(dac)),length(sfs))

    for (i, sf) in enumerate(scumu)
        s_poisson[i] = map(z -> sum(sf[sf[:, 1] .== z, 2:3]), dac)
    end

    d_poisson = map(x -> [sum(x[1:2]), x[3], x[4]], divergence);

    α_observed = fill(zeros(1,length(dac)),length(sfs))
    for (i, sf) in enumerate(scumu)
        tmp = sf[in(dac).(view(sf,:,1)),:];
        α_observed[i] = permutedims(@. round(1 - (divergence[i][2] / divergence[i][1] * tmp[:, 2] / tmp[:, 3]), digits=5))
    end

    return (α_observed, s_poisson, d_poisson)
end

"""
	bootstrap_data!(sfs,divergence,bootstrap)

Bootstrap data following polyDFE manual.

# Arguments
 - `sfs::Matrix{Float64}`: SFS data.
 - `divergence::Matrix{Float64}`: divergence data.
 - `bootstrap::Int64`: number of bootstrapping replicas
# Output
 - `sfs::Vector{Matrix{Float64}}`: boostrapped sfs.
 - `divergence::Vector{Matrix{Float64}}`: divergence data.
"""
function bootstrap_data!(
    sfs::Matrix{Float64},
    divergence::Matrix{Float64},
    bootstrap::Int64,
)
    pr(x::Matrix{Float64}) = hcat(x[:, 1], pois_rand.(x[:, 2:end]))

    sfs = repeat(sfs, bootstrap)
    divergence = repeat(divergence, bootstrap)
    sfs[2:end] .= pr.(sfs[2:end])

    return sfs, divergence
end


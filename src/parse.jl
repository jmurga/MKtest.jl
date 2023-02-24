########################
# Parse data functions # 
########################

"""
	parse_sfs(param;data,gene_list,sfs_columns,div_columns,m_columns,bins,isolines)

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
 - `m_columns::Vector{Int64}`: column numbers corresponding to total number of non-synonymous and synonymous sites. Not require at data. Default values [7,8]
 - `bins::Union{Nothing,Int64}`: size to bin the SFS independtly of the sample size.
# Returns
 - `Vector{Vector{Float64}}`: α_x values. Estimated using the cumulative SFS
 - `Vector{Matrix{Float64}}`: Site Frequency Spectrum (non-cumulative)
 - `Vector{Matrix{Int64}}`: Synonymous and non-synonymous divergence counts
"""
function parse_sfs(param::parameters;
                   data::S,
                   gene_list::Union{Nothing, Vector{S}, S} = nothing,
                   sfs_columns::Vector{Int64} = [3, 5],
                   div_columns::Vector{Int64} = [6, 7],
                   m_columns::Vector{Int64} = [8, 9],
                   bins::Union{Nothing, Int64} = nothing,
                   isolines::Bool = false) where {S <: AbstractString}
    @unpack n, cutoff = param

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
            gene_matrix = permutedims(gene_matrix)
        end

        out = SubDataFrame[]
        for c in eachrow(gene_matrix)
            tmp = @view df[filter(!isnothing, indexin(c, ids)), :]
            push!(out, tmp)
        end

        α, sfs, divergence, m = unzip(map(i -> get_pol_div(i, s_size, cutoff, sfs_columns,
                                                           div_columns, m_columns, bins),
                                          out))
    else
        α, sfs, divergence, m = get_pol_div(df, s_size, cutoff, sfs_columns, div_columns,
                                            m_columns, bins)
        α = [α]
        sfs = [sfs]
        divergence = [divergence]
        m = [m]
    end

    return α, sfs, divergence, m
end

function get_pol_div(df_subset::Union{DataFrame, SubDataFrame},
                     s_size::Int64,
                     cutoff::Vector{Float64},
                     sfs_columns::Vector{Int64},
                     div_columns::Vector{Int64},
                     m_columns::Vector{Int64},
                     bins::Union{Nothing, Int64})
    g(x) = parse.(Float64, filter(y -> y != "", x))

    tmp = split.(df_subset[:, sfs_columns], ",")

    pn = vcat(g.(tmp[:, 1])...)
    ps = vcat(g.(tmp[:, 2])...)

    pn = pn[pn .!= 0]
    ps = ps[ps .!= 0]
    # Round SFS frequencies to the lowest floating value of the dataset independtly of the sample size. Needed to countmap and merge.
    freq = OrderedDict(round.(collect(1:(s_size - 1)) / s_size, digits = 4) .=> 0)

    pn = sort(OrderedDict(countmap(round.(pn, digits = 4))))
    ps = sort(OrderedDict(countmap(round.(ps, digits = 4))))

    # Dn, Ds, Pn, Ps, sfs
    Dn = sum(df_subset[:, div_columns[1]])
    Ds = sum(df_subset[:, div_columns[2]])

    sfs_pn = reduce(vcat, values(merge(+, freq, pn)))
    sfs_ps = reduce(vcat, values(merge(+, freq, ps)))

    if (!isnothing(bins))
        sfs_pn = reduce_sfs(hcat(collect(1:(s_size - 1)), sfs_pn), bins)[:, 2]
        sfs_ps = reduce_sfs(hcat(collect(1:(s_size - 1)), sfs_ps), bins)[:, 2]
        sfs = hcat(collect(1:(bins - 1)) ./ bins, sfs_pn, sfs_ps, 1:(bins - 1))
    else
        sfs = hcat(freq.keys, sfs_pn, sfs_ps, (1:(s_size - 1)))
    end

    # Filtering SFS and changing frequency to DAC
    sfs_flt = sfs[sfs[:, 1] .>= cutoff[1] .&& sfs[:, 1] .<= cutoff[2], [4, 2, 3]]

    scumu = cumulative_sfs(sfs_flt)

    α = round.(1 .- (Ds / Dn .* scumu[:, 2] ./ scumu[:, 3]), digits = 5)

    m::Matrix{Int64} = try
        [sum(df_subset[:, m_columns[1]]) sum(df_subset[:, m_columns[2]])]
    catch
        m = [0 0]
    end

    return (α, sfs_flt, [Dn Ds], m)
end

function data_to_poisson(sfs::Vector{Matrix{Float64}},
                         divergence::Vector{Matrix{Int64}},
                         dac::Vector{Int64})
    f(x::Matrix{Float64}, d::Vector{Int64} = dac) = map(z -> sum(x[x[:, 1] .== z, 2:3]), d)
    al(a, b) = hcat(a[:, 1], @. round(1 - (b[2] / b[1] * a[:, 2] / a[:, 3]), digits = 5))

    scumu = cumulative_sfs.(sfs)
    s_poisson = f.(scumu)
    d_poisson = [[sum(divergence[i][1:2])] for i in eachindex(divergence)]
    α_x = al.(scumu, divergence)
    α_observed = map(x -> permutedims(x[in(dac).(x[:, 1]), 2]), α_x)

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
function bootstrap_data!(sfs::Matrix{Float64},
                         divergence::Matrix{Float64},
                         bootstrap::Int64)
    pr(x::Matrix{Float64}) = hcat(x[:, 1], pois_rand.(x[:, 2:end]))

    sfs = repeat(sfs, bootstrap)
    divergence = repeat(divergence, bootstrap)
    sfs[2:end] .= pr.(sfs[2:end])

    return sfs, divergence
end

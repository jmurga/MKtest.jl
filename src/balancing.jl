"""
    Parse private and shared non-synonymous and synonymous SNPs
"""
function parse_private_share(param::parameters;
    data::S,
    gene_list::Union{Nothing,Vector{S},S} = nothing,
    private_columns::Vector{Int64} = [10,11],
    shared_columns::Vector{Int64} = [12,13],
    ) where {S<:AbstractString}
    @unpack n, cutoff, isolines = param

    if isolines
        s_size = n
    else
        s_size = (n * 2)
    end

    df = CSV.read(data, header = false, delim = '\t', DataFrame)

    # In case m is not defined need to redifine columns positions 
    if size(df,2) <= 11
        private_columns = @. private_columns - 2
        shared_columns  = @. shared_columns -  2
    end

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

        private = map(i -> get_r_s(i,cutoff,private_columns),out)
        shared = map(i -> get_r_s(i,cutoff,shared_columns),out)
        r_s = Vector{Vector{Int64}}(undef, size(gene_matrix, 1))
        for i=1:size(gene_matrix,1)
            r_s[i] = vcat(private[i],shared[i])
        end
    else
        rₙ, rₛ  =   get_r_s(df, cutoff, private_columns)
        sₙ, sₛ  =   get_r_s(df, cutoff, shared_columns)
        r_s = [rₙ,rₛ,sₙ,sₛ]
    end
    return (r_s)
end

function get_r_s(
    df_subset::Union{DataFrame,SubDataFrame},
    cutoff::Vector{Float64},
    columns::Vector{Int64},
)
    g(x) = parse.(Float64, filter(y -> y != "", x))

    tmp = split.(df_subset[:, columns], ",")

    pn = vcat(g.(tmp[:, 1])...)
    ps = vcat(g.(tmp[:, 2])...)

    # pn = pn[pn.!=0]
    # ps = ps[ps.!=0]
    
    
    pn = @. ifelse(pn > 0.5, 1 - pn, pn)
    ps = @. ifelse(ps > 0.5, 1 - ps, ps)

    pn = sum((pn.> cutoff[1]) .& (pn.<=cutoff[2]))
    ps = sum((ps.> cutoff[1]) .& (ps.<=cutoff[2]))

    return [pn,ps]
end

"""
    α_b estimation from Vivak et al. (2022)
"""
function α_b(r_s::Vector{Int64})
    
    rₙ,rₛ,sₙ,sₛ = r_s

    @info "Estimating αᵦ"
    Z = z(r_s)
    @info "Estimating confidence intervals"
    ci = bootstrapper(r_s)

    αᵦ = ifelse(Z > 1, 1 - (1 / Z), 0)
    α_low = ifelse(ci[2] > 1, 1 - (1 / ci[2]), 0)
    α_high = ifelse(ci[end] > 1, 1 - (1 / ci[end]), 0)

    pnᵦ = @. ifelse(αᵦ != 0, αᵦ * sₙ, 0)

    out = DataFrame(:Z=>Z,:αᵦ=>αᵦ,:α_low=>α_low,:α_high=>α_high,:pnᵦ=>pnᵦ,:rₙ=>rₙ,:rₛ=>rₛ,:sₙ=>sₙ,:sₛ=>sₛ,:variance=>ci[1])
    return(out)
end

"""
    α_b estimation from Vivak et al. (2022)
"""
function α_b(r_s::Vector)

    @info "Estimating αᵦ"
    out = @suppress begin
        ThreadsX.map(x -> MKtest.α_b(x),r_s[1:10])
    end
    out = vcat(out...)
    return out

end


"""
    Z estimation from Vivak et al. (2022)
"""
function z(r_s::Vector{Int64})
    rₙ,rₛ,sₙ,sₛ = r_s
    z = if (sₛ>0) & (rₙ>0) & (rₛ>0)
        (sₙ/sₛ)/(rₙ/rₛ)
    else
        0
    end

    return(z)
end

function bootstrapper(r_s::Vector{Int64})

    rₙ,rₛ,sₙ,sₛ = r_s
    total_snps=sum(r_s)

    b_z_list = zeros(1000)

    for i=1:1000

        btstrp = sample(["sN", "sS", "rN", "rS"], Weights([sₙ/total_snps, sₛ/total_snps, rₙ/total_snps, rₛ/total_snps]), total_snps; replace=true)

        b_counts  = countmap(btstrp)

        b_z_list[i] = if(length(b_counts) == 4)
        (b_counts["sN"]/b_counts["sS"]) / (b_counts["rN"]/b_counts["rS"])
        else
            0
        end

    end

    sort!(b_z_list)
    lp = b_z_list[250]
    up = b_z_list[750]
    v  = var(b_z_list)

    # out = DataFrame(:variance=>v,:CI_low=>lp,:CI_high=>up)
    out = vcat(v,lp,up)
    return(out)
end

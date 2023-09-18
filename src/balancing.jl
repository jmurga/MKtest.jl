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

        rₙ, rₛ = unzip(map(i -> get_r_s(i,cutoff,private_columns),out))
        sₙ, sₛ = unzip(map(i -> get_r_s(i,cutoff,shared_columns),out))
    else
        rₙ, rₛ  = get_r_s(df, cutoff, private_columns)
        sₙ, sₛ  = get_r_s(df, cutoff, shared_columns)
    end
    return ([rₙ,rₛ,sₙ,sₛ])
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

    pn = pn[pn.!=0]
    ps = ps[ps.!=0]
    
    
    pn = @. ifelse(pn > 0.5, 1 - pn, pn)
    ps = @. ifelse(ps > 0.5, 1 - ps, ps)

    pn = sum((pn.>=cutoff[1]) .& (pn.<=cutoff[2]))
    ps = sum((ps.>=cutoff[1]) .& (ps.<=cutoff[2]))

    return(pn,ps)
end

"""
    α_b estimation from Vivak et al. (2022)
"""
function α_b(r_s::Vector{Int64})
    
    rₙ,rₛ,sₙ,sₛ = r_s

    @info "Estimating αᵦ"
    Z = z.(rₙ,rₛ,sₙ,sₛ)
    @info "Estimating confidence intervals"
    ci = bootstrapper(rₙ,rₛ,sₙ,sₛ)

    αᵦ = @. ifelse(Z > 1, 1 - (1 / Z), 0)
    α_low = @. ifelse(ci.CI_low > 1, 1 - (1 / ci.CI_low), 0)
    α_high = @. ifelse(ci.CI_high > 1, 1 - (1 / ci.CI_high), 0)

    pnᵦ = @. ifelse(αᵦ != 0, αᵦ * sₙ, 0)

    out = DataFrame(:Z=>Z,:αᵦ=>αᵦ,:α_low=>α_low,:α_high=>α_high,:pnᵦ=>pnᵦ,:rₙ=>rₙ,:rₛ=>rₛ,:sₙ=>sₙ,:sₛ=>sₛ,:variance=>ci.variance)
    return(out)
end

"""
    α_b estimation from Vivak et al. (2022)
"""
function α_b(r_s::Vector)
    
    rₙ,rₛ,sₙ,sₛ = r_s

    Z = z.(rₙ,rₛ,sₙ,sₛ)

    ci = ThreadsX.map((x,y,z,v) -> bootstrapper(x,y,z,v) ,rₙ,rₛ,sₙ,sₛ)
    ci = vcat(ci...)

    αᵦ = @. ifelse(Z > 1, 1 - (1 / Z), 0)
    α_low = @. ifelse(ci.CI_low > 1, 1 - (1 / ci.CI_low), 0)
    α_high = @. ifelse(ci.CI_high > 1, 1 - (1 / ci.CI_high), 0)

    pnᵦ = @. ifelse(αᵦ != 0, αᵦ * sₙ, 0)

    out = DataFrame(:Z=>Z,:Z_low=>ci.CI_low,:Z_high=>ci.CI_high,:αᵦ=>αᵦ,:α_low=>α_low,:α_high=>α_high,:pnᵦ=>pnᵦ,:rₙ=>rₙ,:rₛ=>rₛ,:sₙ=>sₙ,:sₛ=>sₛ,:variance=>ci.variance)
    return(out)
end

"""
    Z estimation from Vivak et al. (2022)
"""
function z(rₙ::Int64,rₛ::Int64,sₙ::Int64,sₛ::Int64)

    z = if (sₛ>0) & (rₙ>0) & (rₛ>0)
        (sₙ/sₛ)/(rₙ/rₛ)
    else
        0
    end

    return(z)
end

"""
    Bootstrap z
"""
function bootstrapper(rₙ::Int64,rₛ::Int64,sₙ::Int64,sₛ::Int64)

    total_snps=rₙ+rₛ+sₙ+sₛ

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
    lp = (b_z_list[250])
    up = (b_z_list[750])
    v  = (var(b_z_list))

    out = DataFrame(:variance=>v,:CI_low=>lp,:CI_high=>up)
    return(out)
end
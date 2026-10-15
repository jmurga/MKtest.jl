########################
# Parse data functions #
########################

"""
    parse_sfs(data; sample_size, cutoff = [0.0, 1.0], gene_list = nothing, mode = :case,
              sfs_columns = [3, 5], div_columns = [6, 7], l_columns = [8, 9], ploidy = 2)

Read a gene table into a raw SFS, the divergence counts and the alpha(x) table.

`data` is a tab-separated table with one row per gene and no header (Uricchio et al.
2019, https://doi.org/10.1038/s41559-019-0890-6), or the same table as a `DataFrame`:

| column | content |
|---|---|
| 1 | gene ID |
| 2 | Pn |
| 3 | derived allele frequencies of the nonsynonymous polymorphisms, `,0.0008,0.0015,` |
| 4 | Ps |
| 5 | derived allele frequencies of the synonymous polymorphisms |
| 6, 7 | Dn, Ds |
| 8, 9 | Ln, Ls, the nonsynonymous and synonymous site counts (optional) |

Frequencies are rounded to 4 digits (3 when `sample_size * ploidy < 100`) and matched
to the grid `i / (sample_size * ploidy)`; zeros are ignored. A frequency off the grid
throws a `DimensionMismatch`, except with a multi-set `gene_list`, where it is skipped.
Without the site-count columns, Ln and Ls are 0.

# Keywords
 - `sample_size::Int64`: number of sampled individuals; the full SFS has
   `sample_size * ploidy - 1` rows.
 - `ploidy::Int`: chromosomes per individual, `1` for inbred lines or haploids.
 - `cutoff::Vector{Float64}`: frequency range `[lo, hi]` kept, both ends included. Rows
   are removed without renumbering, so column 1 stays the true DAC. A trimmed SFS cannot
   be given to [`project`](@ref) or `dfe_ml`.
 - `gene_list`: a `Vector` of gene IDs, or the path of a file (plain or `.gz`) with one
   ID per line, or several space-separated IDs per line (one set per line, detected
   automatically). IDs absent from the table are skipped silently.
 - `mode::Symbol`: `:case` (one SFS for all selected genes) or `:gene` (one SFS per
   gene plus their total as the last element). Ignored for a multi-set file; any other
   value acts as `:case`.
 - `sfs_columns`, `div_columns`, `l_columns`: column numbers of the Pn and Ps
   frequencies, of Dn and Ds, and of Ln and Ls.

# Returns
Three values, in this order:
 - `sfs::Matrix{Float64}`: the raw (non-cumulative) SFS `[DAC Pn Ps]`.
 - `divergence::Vector{Float64}`: `[Dn, Ds, Ln, Ls]`.
 - `alpha_x::Matrix{Float64}`: `[DAC frequency alpha(x)]`, alpha(x) from the tail sums of
   the (cut) SFS, `1 - (Ds / Dn) Pn(>= x) / Ps(>= x)`, rounded to 5 digits.

With a multi-set `gene_list` or `mode = :gene`, each of the three is a `Vector` with one
element per set or gene (plus the total for `:gene`).
"""
function parse_sfs(
    data::Union{AbstractString,DataFrame};
    sample_size::Int64,
    cutoff::Vector{Float64}=[0.0, 1.0],
    gene_list::Union{Nothing,Vector{<:AbstractString},AbstractString}=nothing,
    mode::Symbol=:case,
    sfs_columns::Vector{Int64}=[3, 5],
    div_columns::Vector{Int64}=[6, 7],
    l_columns::Vector{Int64}=[8, 9],
    ploidy::Int=2,
)

    s_size = sample_size * ploidy

    if data isa AbstractString
        df = CSV.read(data, header=false, delim='\t', DataFrame)
    elseif data isa DataFrame
        df = data
    end

    if !isnothing(gene_list)
        ids = @view df[:, 1]

        @assert gene_list isa AbstractString || gene_list isa Vector{<:AbstractString} "Please input a valid gene_list file or vector containing gene ids"

        if gene_list isa AbstractString
            gene_matrix = Matrix(CSV.read(gene_list, header=false, DataFrame))
        else
            gene_matrix = reshape(gene_list, :, 1)
        end

        m, n = size(gene_matrix)

        if n > 1
            # -- :control (auto-detected): each row is a bootstrap sample -------
            pp = _preparse_df(df, s_size, cutoff, sfs_columns, div_columns, l_columns)
            α, sfs, divergence = unzip(
                ThreadsX.mapi(
                    i -> _get_pol_div_preparsed(view(gene_matrix, i, :), pp), 1:m,
                ),
            )

        elseif mode == :gene
            # -- :gene: per-gene SFS + aggregated total --------------------------
            all_genes = vec(gene_matrix)
            df_filtered = filter_df(df, ids, all_genes)

            rows = collect(eachrow(df_filtered))
            α, sfs, divergence = unzip(
                ThreadsX.mapi(
                    row -> get_pol_div(row, s_size, cutoff, sfs_columns, div_columns, l_columns), rows,
                ),
            )

            α_total, sfs_total, div_total = get_pol_div(
                df_filtered, s_size, cutoff, sfs_columns, div_columns, l_columns,
            )
            push!(α, α_total)
            push!(sfs, sfs_total)
            push!(divergence, div_total)

        else
            # -- :case (default): single merged SFS for all listed genes ---------
            all_genes = vec(gene_matrix)
            df_filtered = filter_df(df, ids, all_genes)
            α, sfs, divergence = get_pol_div(
                df_filtered, s_size, cutoff, sfs_columns, div_columns, l_columns,
            )
        end

    elseif mode == :gene
        # -- :gene, no gene_list: per-row SFS over the whole df + aggregated total -
        rows = collect(eachrow(df))
        α, sfs, divergence = unzip(
            ThreadsX.mapi(
                row -> get_pol_div(row, s_size, cutoff, sfs_columns, div_columns, l_columns), rows,
            ),
        )
        α_total, sfs_total, div_total = get_pol_div(
            df, s_size, cutoff, sfs_columns, div_columns, l_columns,
        )
        push!(α, α_total)
        push!(sfs, sfs_total)
        push!(divergence, div_total)

    else
        α, sfs, divergence =
            get_pol_div(df, s_size, cutoff, sfs_columns, div_columns, l_columns)
    end

    return sfs, divergence, α
end

# Not duplicated entries in ids or df, so get first index using indexin is valid
filter_df(df::DataFrame, ids::SubArray, row::Union{SubArray,Array,DataFrameRow}) = view(df, filter(!isnothing, indexin(row, ids)), :)

# Control sets (many genes per row): the DAF strings and the divergence and site columns
# are parsed once for the whole table, and each row only adds up counts. Much faster than
# calling get_pol_div per row, with the same result.

struct _PreparsedDF
    gene_index::Dict{String,Int}       # gene_id  → df row index
    pn_data::Vector{Vector{Float64}}   # non-zero rounded pn DAF per gene
    ps_data::Vector{Vector{Float64}}   # non-zero rounded ps DAF per gene
    dn::Vector{Float64}
    ds::Vector{Float64}
    ln::Vector{Float64}
    ls::Vector{Float64}
    freq_keys::Vector{Float64}         # round(i/s_size, 4) for i in 1:s_size-1
    freq_idx::Dict{Float64,Int}        # freq_value → histogram bin index
    dac_range::Vector{Int}             # DAC indices that survive the cutoff filter
    s_size::Int                        # haploid sample size; the frequency denominator
end

function _preparse_df(
    df::DataFrame,
    s_size::Int,
    cutoff::Vector{Float64},
    sfs_columns::Vector{Int64},
    div_columns::Vector{Int64},
    l_columns::Vector{Int64},
)
    ng = nrow(df)
    gene_ids = df[:, 1]
    gene_index = Dict{String,Int}(string(gene_ids[i]) => i for i in 1:ng)

    pn_data = Vector{Vector{Float64}}(undef, ng)
    ps_data = Vector{Vector{Float64}}(undef, ng)
    dn_arr = Vector{Float64}(undef, ng)
    ds_arr = Vector{Float64}(undef, ng)
    ln_arr = Vector{Float64}(undef, ng)
    ls_arr = Vector{Float64}(undef, ng)

    # Frequency grid shared across all genes, rounded as in get_pol_div
    dgts = s_size >= 100 ? 4 : 3
    freq_keys = round.(collect(1:(s_size-1)) ./ s_size, digits=dgts)
    freq_idx = Dict{Float64,Int}(freq_keys[i] => i for i in eachindex(freq_keys))

    Threads.@threads for i in 1:ng
        row = df[i, :]
        pn_raw = split(row[sfs_columns[1]], ",")
        ps_raw = split(row[sfs_columns[2]], ",")
        pn = parse.(Float64, filter(!isempty, pn_raw))
        ps = parse.(Float64, filter(!isempty, ps_raw))
        pn_data[i] = round.(filter(!iszero, pn), digits=dgts)
        ps_data[i] = round.(filter(!iszero, ps), digits=dgts)
        dn_arr[i] = row[div_columns[1]]
        ds_arr[i] = row[div_columns[2]]
        try
            ln_arr[i] = row[l_columns[1]]
            ls_arr[i] = row[l_columns[2]]
        catch
            ln_arr[i] = 0.0
            ls_arr[i] = 0.0
        end
    end

    dac_range = [i for i in eachindex(freq_keys)
                 if freq_keys[i] >= cutoff[1] && freq_keys[i] <= cutoff[2]]

    return _PreparsedDF(gene_index, pn_data, ps_data, dn_arr, ds_arr, ln_arr, ls_arr,
        freq_keys, freq_idx, dac_range, s_size)
end

function _get_pol_div_preparsed(
    gene_row::AbstractVector,
    pp::_PreparsedDF,
)::Tuple{Matrix{Float64},Matrix{Float64},Vector{Float64}}
    nb = length(pp.freq_keys)
    pn_counts = zeros(Int, nb)
    ps_counts = zeros(Int, nb)
    Dn = 0.0
    Ds = 0.0
    Ln = 0.0
    Ls = 0.0

    for g in gene_row
        idx = get(pp.gene_index, string(g), 0)
        idx == 0 && continue
        for v in pp.pn_data[idx]
            fi = get(pp.freq_idx, v, 0)
            fi > 0 && (pn_counts[fi] += 1)
        end
        for v in pp.ps_data[idx]
            fi = get(pp.freq_idx, v, 0)
            fi > 0 && (ps_counts[fi] += 1)
        end
        Dn += pp.dn[idx]
        Ds += pp.ds[idx]
        Ln += pp.ln[idx]
        Ls += pp.ls[idx]
    end

    dr = pp.dac_range
    sfs_flt = Float64.(hcat(dr, pn_counts[dr], ps_counts[dr]))

    scumu = cumulative_sfs(sfs_flt)
    α = round.(1 .- (Ds / Dn .* scumu[:, 2] ./ scumu[:, 3]), digits=5)
    return (_alpha_table(sfs_flt, α, pp.s_size), sfs_flt,
        [Dn, Ds, round(Ln, digits=2), round(Ls, digits=2)])
end


# alpha(x) with its coordinates: [dac, frequency, alpha_x]. Column 1 of the SFS is the
# derived allele count (a cutoff removes rows without renumbering), so frequency = dac / s_size.
function _alpha_table(sfs_flt::AbstractMatrix, α::AbstractVector, s_size::Integer)
    dac = Float64.(view(sfs_flt, :, 1))
    return hcat(dac, dac ./ s_size, Float64.(α))
end

# A column as a vector, for one gene (DataFrameRow) or several (DataFrame)
_col_as_vec(row::DataFrameRow, col::Int) = [row[col]]
_col_as_vec(df::Union{DataFrame,SubDataFrame}, col::Int) = df[:, col]

function get_pol_div(
    df_subset::Union{DataFrame,SubDataFrame,DataFrameRow},
    s_size::Int64,
    cutoff::Vector{Float64},
    sfs_columns::Vector{Int64},
    div_columns::Vector{Int64},
    l_columns::Vector{Int64},
)::Tuple{Matrix{Float64},Matrix{Float64},Vector{Float64}}

    pn_col = _col_as_vec(df_subset, sfs_columns[1])
    ps_col = _col_as_vec(df_subset, sfs_columns[2])
    # Comma-separated DAF strings
    pn_raw = collect(Iterators.flatten(split(s, ",") for s in skipmissing(pn_col)))
    ps_raw = collect(Iterators.flatten(split(s, ",") for s in skipmissing(ps_col)))
    pn = parse.(Float64, view(pn_raw, pn_raw .!= ""))
    ps = parse.(Float64, view(ps_raw, ps_raw .!= ""))
    pn = view(pn, pn .!= 0)
    ps = view(ps, ps .!= 0)

    # Round frequencies to match the sample size grid before counting
    dgts = s_size >= 100 ? 4 : 3
    freq = OrderedDict(round.(collect(1:(s_size-1)) / s_size, digits=dgts) .=> 0)

    c_pn = countmap(round.(pn, digits=dgts))
    c_ps = countmap(round.(ps, digits=dgts))

    # Dn, Ds, Pn, Ps, sfs
    Dn = sum(_col_as_vec(df_subset, div_columns[1]))
    Ds = sum(_col_as_vec(df_subset, div_columns[2]))
    sfs_pn = mergewith!(+, copy(freq), c_pn).vals
    sfs_ps = mergewith!(+, copy(freq), c_ps).vals

    sfs = hcat(freq.keys, sfs_pn, sfs_ps, (1:(s_size-1)))

    # Filtering SFS and changing frequency to DAC
    sfs_flt = sfs[view(sfs, :, 1).>=cutoff[1].&&view(sfs, :, 1).<=cutoff[2], [4, 2, 3]]

    scumu = cumulative_sfs(sfs_flt)

    α = round.(1 .- (Ds / Dn .* scumu[:, 2] ./ scumu[:, 3]), digits=5)
    α = _alpha_table(sfs_flt, α, s_size)

    l = try
        [round(sum(_col_as_vec(df_subset, l_columns[1])), digits=2), round(sum(_col_as_vec(df_subset, l_columns[2])), digits=2)]
    catch
        [0, 0]
    end

    return (α, sfs_flt, vcat([Dn, Ds], l))
end

function data_to_poisson(
    sfs::Vector{Matrix{Float64}},
    divergence::Vector,
    dac::Vector{Int64}
)
    scumu = cumulative_sfs.(sfs)
    s_poisson = fill(zeros(length(dac)), length(sfs))

    for (i, sf) in enumerate(scumu)
        s_poisson[i] = map(z -> sum(sf[sf[:, 1].==z, 2:3]), dac)
    end

    d_poisson = try
        map(x -> [sum(x[1:2]), x[3], x[4]], divergence)
    catch
        map(x -> [sum(x[1:2])], divergence)
    end

    α_observed = fill(zeros(1, length(dac)), length(sfs))
    for (i, sf) in enumerate(scumu)
        tmp = sf[in(dac).(view(sf, :, 1)), :]
        α_observed[i] = permutedims(@. round(1 - (divergence[i][2] / divergence[i][1] * tmp[:, 2] / tmp[:, 3]), digits=5))
    end

    return (α_observed, s_poisson, d_poisson)
end


"""
    bootstrap_data(sfs, divergence, bootstrap)

Poisson resamples of one SFS, as in the polyDFE manual: every count of columns 2 to end
is replaced by a Poisson draw with that mean; column 1 (the DAC) is kept.

The divergence is not resampled: each replicate gets the input `divergence` (the same
object, not a copy).

# Arguments
 - `sfs::Matrix{Float64}`: raw SFS `[DAC Pn Ps]`.
 - `divergence::AbstractVecOrMat`: `[Dn, Ds, Ln, Ls]`.
 - `bootstrap::Int64`: number of replicates.

# Returns
 - `(sfs_b, divergence_b)`: a `Vector{Matrix{Float64}}` of `bootstrap` resampled SFS and
   a `Vector` holding `divergence` `bootstrap` times.
"""
function bootstrap_data(
    sfs::Matrix{Float64},
    divergence::AbstractVecOrMat,
    bootstrap::Int64,
)
    pr(x::Matrix{Float64}) = hcat(x[:, 1], pois_rand.(x[:, 2:end]))

    sfs_b = fill(sfs, bootstrap)
    divergence = fill(divergence, bootstrap)
    sfs_b .= pr.(sfs_b)

    return sfs_b, divergence
end

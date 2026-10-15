"""
    bootstrap_parameters(; data, annotations, distance, max_reps, tolerance, iterations,
                         confounding_factors, output, filter_hla_hist)

Settings of `bootstrap`, the matched control-set builder of Enard & Petrov (2020) and
Di et al. (2021). A keyword (`Parameters.@with_kw`) mutable struct; every field has a
default.

# Fields
 - `data::String`: file with the case gene IDs, one per line, no header. Default
   `"ensembl_list.txt"`.
 - `annotations::String`: tab-separated BED file without header, at least five columns
   (chromosome, start, end, gene ID, one more column). The gene IDs in column 4 define the
   gene universe. Default `"ensembl_gene_coords_v109.bed"`.
 - `distance::Int64`: minimum distance in bp from a candidate control gene to the closest
   case gene, as reported by `bedtools closest -d` (see `get_distance`). Default `1`.
 - `max_reps::Int64`: a control gene may make up at most `max_reps / n_case` of the pooled
   control genes of one iteration, i.e. on average at most `max_reps` copies per control
   set. Default `3`.
 - `tolerance::Float64`: relative tolerance of the matching: the running average of every
   confounding factor must stay within `[1 - tolerance, 1 + tolerance]` times the case average.
   Default `0.05`.
 - `iterations::Int64`: number of independent runs; each gives 100 control sets. Default `10`
   (1000 sets).
 - `confounding_factors::String`: file without header, gene ID then one numeric column per
   factor. Default `"confounding_factors.txt"`.
 - `output::String`: output prefix, without extension; `bootstrap` writes
   `<output>_case.txt` and `<output>_control.txt`. Default `"test"`.
 - `filter_hla_hist::Bool`: remove the human HLA and histone genes (Ensembl gene IDs, lists
   `hla_genes` and `hist_genes` in `src/bootstrap.jl`) from case and controls. Default `true`.
"""
@with_kw mutable struct bootstrap_parameters
    data::String = "ensembl_list.txt"
    annotations::String = "ensembl_gene_coords_v109.bed"
    distance::Int64 = 1
    max_reps::Int64 = 3
    tolerance::Float64 = 0.05
    iterations::Int64 = 10
    confounding_factors::String = "confounding_factors.txt"
    output::String = "test"
    filter_hla_hist::Bool = true
end


function filter_case_control(x::Vector{String}, y::Vector{String}, distance::Vector{String})

    # setdiff!(x, hla_hist)
    # x_case::Vector{String} = intersect(x, y[:, 1])
    y_control::Vector{String} = setdiff(y, x)

    y_control = intersect(x,intersect(y, distance))

    return (string.(y_control))
end

"""
    get_valid(x, y, distance)

Unused internal helper: computes `x ∩ y ∩ distance` and returns `nothing`.
"""
function get_valid(x::Vector{<:AbstractString}, y::Vector{<:AbstractString}, distance::Vector{<:AbstractString})

    y_control = intersect(x,intersect(y, distance))

    return
end

"""
    get_factors(case_set, control_set, factors_raw)

Internal: the case-set average of each confounding factor, the control factor matrix, a
gene -> factor vector `Dict` of the controls and a gene -> 0 use counter.
"""
function get_factors(
    case_set::Vector{String},
    control_set::Vector{String},
    factors_raw::Matrix,
)
    control_values = factors_raw[(@view factors_raw[:, 1]).∈[control_set], 1:end]
    # control_dict = Dict{String,Vector{Float64}}()

    # for row in eachrow(control_values)
    #     control_dict[row[1]] = row[2:end]
    # end

    control_number::Dict{String,Int64} = Dict(control_set .=> 0)
    control_dict::Dict{String,Vector{Float64}} =  Dict(control_values[:,1] .=> eachrow(control_values[:,2:end]));

    case_avg::Vector{Float64} =
        vec(mean(factors_raw[(@view factors_raw[:, 1]).∈[case_set], 2:end], dims = 1))

    return (
        case_avg,
        convert(Matrix{Float64}, @view control_values[:, 2:end]),
        control_dict,
        control_number,
    )
end

"""
    check_limits(control_avg, tol)

Internal: per-factor bounds `[(1 - tol) * avg, (1 + tol) * avg]` (a two-column matrix) and a
vector of `"start"` flags. Called with the case averages.
"""
function check_limits(control_avg::Vector{Float64}, tol::Float64)
    l_u::Matrix{Float64} = hcat(control_avg, control_avg)
    lims::Vector{String} = fill("start", length(control_avg))

    for i::Int64 in eachindex(control_avg)
        if lims[i] == "higher"
            l_u[i, 1] = (1 - tol) * control_avg[i]
            l_u[i, 2] = (1 + 1 * tol) * control_avg[i]
        else
            l_u[i, 1] = (1 - 1 * tol) * control_avg[i]
            l_u[i, 2] = (1 + tol) * control_avg[i]
        end
    end

    return (lims, l_u)
end

function get_samples(case_set::Vector{String}, control_boot::Vector{String})
    sc_t = length(control_boot)
    case_number = length(case_set)

    out = String[]
    for p::Int64 = 0:99
        sup_ind = sc_t - 1 - p * case_number
        inf_ind = sc_t - 1 - p * case_number - case_number + 1
        iter = control_boot[inf_ind:sup_ind]

        z = p + 1
        s = "sample_" * string(z) * " " * join(iter, " ")
        push!(out, s)
    end

    return (out)
end

"""
    get_distance(data::String, annotation::String) -> DataFrame

Distance in bp from every gene of `annotation` to the closest gene listed in `data`, with
`bedtools closest -d` (installed into the CondaPkg environment from bioconda on first use).

`data` is the case gene list (one ID per line), `annotation` a tab-separated BED file with
at least five columns (chromosome, start, end, gene ID, ...). Each gene is reduced to a 1 bp
interval at its midpoint `int((start + end) / 2)`, so the distance is between gene
midpoints, not gene bodies; case genes get 0. Case genes are selected with `fgrep -f`
(substring match of the IDs). The search runs per chromosome; a chromosome without case
genes is skipped, so its genes get no row. Needs `fgrep`, `grep -P`, `awk`, `sort` and `cut`
on the PATH.

# Returns
 - `DataFrame` without header names: `Column1` the gene ID, `Column2` the distance in bp
   (tied closest genes give one row each).
"""
function get_distance(data::String, annotation::String)
    if (isnothing(CondaPkg.which("bedtools")))
        CondaPkg.add("bedtools", channel = "bioconda")
    end

    bedtools = CondaPkg.which("bedtools")

    tmp = tempname()
    out = IOBuffer()

    run(
        pipeline(
            `fgrep -f $data $annotation`,
            `awk  -F'\t' 'BEGIN {OFS = FS} {print $1,int(($2+$3)/2),int(($2+$3)/2+1),$4,$5}'`,
            `sort -k1,2V`,
            tmp,
        ),
    )

    run(pipeline(`cut -f1 $annotation`, stdout = out))

    nchr = unique(split(String(take!(out)), "\n")[1:(end-1)])
    distance = tempname()
    gene_1 = tempname()
    gene_2 = tempname()

    for i in nchr
        try
            run(pipeline(`grep -P "^$i\t" $tmp`, stdout = gene_1))
            run(
                pipeline(
                    `grep -P "^$i\t" $annotation`,
                    `awk -F'\t' 'BEGIN {OFS = FS} {print $1,int(($2+$3)/2),int(($2+$3)/2+1),$4,$5}'`,
                    `sort -k2,3n`,
                    gene_2,
                ),
            )
            run(
                pipeline(
                    pipeline(
                        `$bedtools closest -d -a $gene_2 -b $gene_1`,
                        `cut -f4,11`,
                        `sort`,
                    ),
                    stdout = distance,
                    append = true,
                ),
            )
        catch e
            continue
        end
    end

    df = CSV.read(distance, DataFrame, delim = '\t', header = false)

    rm.([tmp, gene_1, gene_2])

    return df
end

"""
    get_bootstrap(case_set, control_set; factors_raw, max_reps, tolerance) -> Vector

Draw one pool of matched control genes for `case_set` from `control_set` (Enard & Petrov
2020).

Pairs of distinct candidate genes are drawn at random and accepted when the running
average of every confounding factor stays within `tolerance` (relative) of the case
average and neither gene exceeds a share of `max_reps / length(case_set)` of the pool. The
running average starts from a seed of 100 pseudo-genes at the case average; each accepted
pair removes one of them. Drawing stops at `100 + 110 * length(case_set)` genes, or after 5
million draws in a row without an accepted pair.

`bootstrap` calls the equivalent `get_bootstrap_noalloc`, not this function.

# Arguments
 - `case_set::Vector{String}`, `control_set::Vector{String}`: gene IDs.
 - `factors_raw::Matrix`: gene ID in column 1, one confounding factor per further column.
 - `max_reps::Int64`, `tolerance::Float64`: as in [`bootstrap_parameters`](@ref).

# Returns
 - The accepted control genes in order of acceptance (with repeats), or an empty vector if
   fewer than 1000 were accepted. `get_samples` cuts the `100 * length(case_set)` genes
   before the last one into 100 control sets; the earlier ones are burn-in.
"""
function get_bootstrap(
    case_set::Vector{String},
    control_set::Vector{String};
    factors_raw::Matrix,
    max_reps::Int64,
    tolerance::Float64,
)

    # Get cont
    case_number    = length(case_set)
    control_number = length(control_set)

    # Get averages
    case_avg, factors, factors_dict, used_number =
        get_factors(case_set, control_set, factors_raw)
    control_avg                                  = vec(mean(factors, dims = 1))

    # high_low
    high_low, inf_sup    = check_limits(case_avg, tolerance)
    flt_avg              = case_avg .> control_avg
    high_low[flt_avg]   .= "higher"
    high_low[.!flt_avg] .= "lower"

    # Increase control dataset to iter to create the control sets per se. Creating a big list speeds up the process because it avoids repeating the slow start of adding control genes to the "fake" seed.
    # gene_1, gene_2 = control_iteration(x, y)

    ######################
    # Iter combinations
    ######################
    previous_gn  = 0
    same_counter = 0

    #Control sets are created by packs of 100 to speed up control creation while keeping memory usage in check.
    init_fake = fake_seed_num = 100

    # to keep the control sets close the averages for confounding factors in the group of interest, we later require that there some level of alternance with the newly added control genes swithcing between increasing the overall average and decreasing it.
    lim             = size(factors, 2)

    direction       = fill("start", lim)
    current_factors = case_avg

    added_nonvips   = 0
    good_control    = String[]

    sc_gn           = length(good_control)

    while sc_gn < (init_fake + (10 + 100) * case_number)
        # the two control genes to choose to add to the current control set.
        gene_1 = rand(control_set)
        gene_2 = rand(control_set)

        while gene_1 == gene_2
            gene_2 = rand(control_set)
        end

        factor_value_1 = factors_dict[gene_1]
        factor_value_2 = factors_dict[gene_2]

        sc_fc = length(factor_value_1)
        sc_t = length(good_control)

        matching::Int64 = 0

        test_value_1 =
            (
                fake_seed_num * case_avg +
                current_factors * sc_t +
                factor_value_1 +
                factor_value_2
            ) / (sc_t + 2 + fake_seed_num)

        test_value_2 =
            (current_factors * sc_t + factor_value_1 + factor_value_2) / (sc_t + 2)

        matching =
            sum((test_value_1 .>= inf_sup[:, 1]) .&& (test_value_1 .<= inf_sup[:, 2]))

        if matching >= lim
            test_dir = similar(direction)
            flt_dir = test_value_1 .>= current_factors
            test_dir[flt_dir] .= "more"
            test_dir[.!flt_dir] .= "less"
            other_dir = sum(test_dir .!= direction)
            test_ln = length(good_control) + 1

            if (
                (other_dir >= 0) &&
                (used_number[gene_1] / test_ln <= (max_reps / case_number)) &&
                (used_number[gene_2] / test_ln <= (max_reps / case_number))
            )
                if (fake_seed_num > 0)
                    fake_seed_num -= 1
                end

                push!(good_control, gene_1)
                push!(good_control, gene_2)

                used_number[gene_1] += 1
                used_number[gene_2] += 1

                # Update current_factors
                current_factors = copy(test_value_1)

                # Update directions
                copy!(direction, test_dir)
            end
        end

        sc_gn = length(good_control)

        if (sc_gn == previous_gn)
            same_counter += 1
        elseif sc_gn > previous_gn
            same_counter = 0
        end

        if (same_counter >= 5000000)
            break
        end

        previous_gn = sc_gn
    end

    if (length(good_control) < 1000)
        return []
    else
        return (good_control)
    end
end


function get_bootstrap_noalloc(
    case_set::Vector{String},
    control_set::Vector{String};
    factors_raw::Matrix,
    max_reps::Int64,
    tolerance::Float64,
)


    case_number    = length(case_set)
    control_number = length(control_set)

    case_avg, factors, factors_dict, used_number =
        get_factors(case_set, control_set, factors_raw)
    control_avg                 = vec(mean(factors, dims = 1))
    high_low, inf_sup           = check_limits(case_avg, tolerance)
    flt_avg                     = case_avg .> control_avg
    high_low[flt_avg]  .= "higher"
    high_low[.!flt_avg] .= "lower"

    lim      = size(factors, 2)
    max_frac = max_reps / case_number

    # CHANGED: reserve capacity upfront for pushes
    N_fake   = 100
    N_target = N_fake + (10 + 100)*case_number
    good_control = Vector{String}()
    # CHANGED: sizehint to avoid reallocs
    sizehint!(good_control, N_target)

    #  CHANGED: pre-allocate/reuse scratch buffers outside the loop
    current_factors = copy(case_avg)             # holds the evolving control-avg
    direction       = fill("start", lim)         # initial directions
    test_value_1    = similar(case_avg)          # buffer for new avg w/ fake seed
    test_value_2    = similar(case_avg)          # buffer for new avg w/o fake seed
    test_dir        = similar(direction)         # buffer for per-factor direction

    prev_len  = 0
    same_ctr  = 0

    #  main loop
    while length(good_control) < N_target
        # pick two distinct controls
        g1 = rand(control_set)
        g2 = rand(control_set)
        while g2 === g1
            g2 = rand(control_set)
        end

        v1 = factors_dict[g1]
        v2 = factors_dict[g2]
        sc_t = length(good_control)

        # CHANGED: in-place vector math via @. to avoid temporaries
        @. test_value_1 = (N_fake*case_avg + sc_t*current_factors + v1 + v2) / (N_fake + sc_t + 2)
        @. test_value_2 = (sc_t*current_factors + v1 + v2) / (sc_t + 2)

        # CHANGED: manual loop for matching instead of array ops
        matches = 0
        for i in 1:lim
            matches += (test_value_1[i] ≥ inf_sup[i,1] && test_value_1[i] ≤ inf_sup[i,2])
        end

        if matches ≥ lim
            #  CHANGED: in-place fill of test_dir, no Boolean temporary
            for i in 1:lim
                test_dir[i] = test_value_1[i] ≥ current_factors[i] ? "more" : "less"
            end

            #  CHANGED: lazy zip+count instead of non-existent count(!=,…)
            other_dir = count(t -> t[1] != t[2], zip(test_dir, direction))

            next_len  = sc_t + 1

            if other_dir ≥ 0 &&
               used_number[g1]/next_len ≤ max_frac &&
               used_number[g2]/next_len ≤ max_frac

                if N_fake > 0
                    N_fake -= 1
                end

                # CHANGED: single push! of both genes, no reallocs due to hint
                push!(good_control, g1, g2)
                used_number[g1] += 1
                used_number[g2] += 1

                # CHANGED: update in-place, reuse buffers
                copy!(current_factors, test_value_1)
                copy!(direction, test_dir)
            end
        end

        # break-out to avoid infinite loop
        cur_len = length(good_control)
        same_ctr = (cur_len == prev_len) ? same_ctr+1 : 0
        prev_len = cur_len
        if same_ctr ≥ 5_000_000
            break
        end
    end

    return length(good_control) < 1000 ? String[] : good_control
end

"""
    bootstrap(param::bootstrap_parameters)

Build matched control sets for a case gene list and write them to disk (Enard & Petrov
2020; Di et al. 2021). See [`bootstrap_parameters`](@ref) for the settings.

Steps: read the case IDs (`data`), the confounding factors and the `annotations` gene IDs;
compute the distance of every gene to the closest case gene with `get_distance`; keep as
candidate controls the genes present in `annotations` and in the factor file, at least
`distance` bp from any case gene and not in the case set; keep only case genes present in
both files; drop HLA and histone genes if `filter_hla_hist`. Then run `iterations`
independent matched draws in parallel (`ThreadsX.map`, one per iteration), each cut into
100 control sets of `length(case)` genes, for `iterations * 100` sets in total.

Warns when there are fewer than 1.5 candidate controls per case gene.

# Files written
 - `<output>_control.txt`: one control set per line, `sample_i GENE GENE ...`
   (space-separated; `i` restarts at 1 in each iteration); the format `parse_sfs` reads as
   `gene_list`. An existing file is replaced.
 - `<output>_case.txt`: the case genes kept after the filters, one per line.

If any iteration fails to find enough matching controls (see `get_bootstrap`), it logs an
error and writes nothing; it does not throw.

# Returns
 - `nothing`.
"""
function bootstrap(param::bootstrap_parameters)
    @unpack data, annotations, distance, max_reps, tolerance, iterations, confounding_factors, output, filter_hla_hist = param

    @info "Opening data"

    case        = vec(CSV.read(data, header = false, delim = '\t', Tables.matrix, stringtype = String));

    factors_raw = CSV.read(confounding_factors, header = false, Tables.matrix);
    factors_id  = String.(@view factors_raw[:, 1]);
    factors     = Float64.(@view factors_raw[:, 2:end]);
    ann         = CSV.read(annotations,header=false,DataFrame)[:,4];

    @info "Estimating distance between genes"
    distance_raw = get_distance(data, annotations);

    #Filter case, control and distance
    distance_min = distance_raw[distance_raw[:, 2] .>= distance, 1];

    # Get all valid factors and gene set, case excluded unless dist == 0
    control = intersect(ann, factors_id, distance_min);
    case    = intersect(case,intersect(ann, factors_id));
    setdiff!(control,case);

    if filter_hla_hist
        setdiff!(case,hla_genes,hist_genes);
        setdiff!(control,hla_genes,hist_genes);
    end

    case_number    = length(case);
    control_number = length(control);

    @info "There are " *
          string(case_number) *
          " genes of interest and " *
          string(control_number) *
          " potential control genes at distance of at least " *
          string(distance) *
          " bases."

    if (control_number <= (1.5 * case_number))
        @warn "Warning! The number of control genes is less than 1.5 times the number of genes of interest. FDR may be high as a result."
    end

    # using ThreadsX
    n_case    = fill(case, iterations)
    n_control = fill(control, iterations)

    @info "Running bootstrap"
    good_control = ThreadsX.map(
        (x, y) -> get_bootstrap_noalloc(
            x,
            y,
            factors_raw = factors_raw,
            max_reps = max_reps,
            tolerance = tolerance,
        ),
        n_case,
        n_control,
    )


    if any(isempty.(good_control))
        @error "The bootstrap struggled too much to find matching controls. Stopping bootstrap test here. Please restart the entire pipeline for safety at input $data changing the parameters"
    else
        out_samples = ThreadsX.map(x -> get_samples(case, x), good_control);

        try
            rm(output * "_control.txt")
        catch e
            touch(output * "_control.txt")
        end

        io = Base.open(output * "_control.txt", "a+")
        for i in vcat(out_samples...)
            println(io, i)
        end

        close(io)

        io = Base.open(output * "_case.txt", "w+")
        for i in case
            println(io, i)
        end
        close(io)
    end
end

"""
    empirical_pvalues(observed_values::Vector{Float64}) -> Vector{Float64}

Rank-based empirical p-value of each value against the whole vector, upper tail:
`ordinalrank(-x) / n`. The largest value gets `1 / n`, the smallest `1`; ties are broken by
position. For the p-value of a case against its control sets use
[`empirical_test`](@ref).
"""
function empirical_pvalues(observed_values::Vector{Float64})
  return ordinalrank(observed_values * -1)/length(observed_values)
end

"""
    case_control_table(case, control; analysis)

Stack the estimates of one case set and its control sets into a long table, the input of
[`empirical_test`](@ref) and [`plot_violins_empirical`](@ref). Build one table per
analysis and `vcat` them.

# Arguments
 - `case`, `control`: [`abcmk`](@ref) outputs (their `:inference` entry is used) or
   DataFrames with one row per set (e.g. `fww_mk` on a vector of SFS). `case` must hold
   exactly one row: wrap a single result with `DataFrame([result])`.
 - `analysis`: name of the comparison, e.g. the gene category (stored as a String).

# Returns
 - `DataFrame`: `analysis`, `type` (`"case"`/`"control"`), `set` (1 for the case, 1..n for
   the controls), then the columns both inputs share (α, α_strong, α_weak, ... for ABC-MK).

Throws `ArgumentError` if `case` has more than one row, if the inputs share no column, or
if a Dict input has no `:inference` entry.
"""
function case_control_table(case, control; analysis)
    df_case = _inference_df(case)
    df_control = _inference_df(control)
    nrow(df_case) == 1 || throw(ArgumentError(
        "the case must be one set, got $(nrow(df_case)) rows; run each case separately"))
    cols = intersect(propertynames(df_case), propertynames(df_control))
    isempty(cols) && throw(ArgumentError("case and control share no parameter column"))
    n = nrow(df_control)
    meta = DataFrame(analysis = fill(string(analysis), n + 1),
                     type = vcat("case", fill("control", n)),
                     set = vcat(1, 1:n))
    return hcat(meta, vcat(df_case[:, cols], df_control[:, cols]))
end

_inference_df(x::AbstractDict) = haskey(x, :inference) ? DataFrame(x[:inference]) :
    throw(ArgumentError("expected an `abcmk` output with an :inference entry"))
_inference_df(x::AbstractDataFrame) = DataFrame(x)

# The three α or the three ω rows of Figure 4 style comparisons.
_empirical_params(omega::Bool) = omega ? [:ωₐ, :ωₐ_strong, :ωₐ_weak] : [:α, :α_strong, :α_weak]

"""
    empirical_test(tbl; omega, params, alternative)

Empirical p-value of each case against its own controls, per analysis and parameter:
`p = (1 + #{control >= case}) / (1 + n)` for `alternative = :greater` (Davison & Hinkley
1997), so it is never 0 and its minimum is `1 / (1 + n)`. `:less` counts `control <= case`;
`:two_sided` is `min(1, 2 * min(p_greater, p_less))`. Non-finite control values are dropped
and `n` counts the finite ones.

# Arguments
 - `tbl::DataFrame`: output of [`case_control_table`](@ref), one or more analyses.
 - `omega::Bool`: test ωₐ, ωₐ_strong, ωₐ_weak instead of α, α_strong, α_weak. Default `false`.
 - `params`: parameter columns to test; overrides `omega`.
 - `alternative::Symbol`: `:greater` (case above controls, default), `:less` or `:two_sided`.

# Returns
 - `DataFrame`: `analysis, param, case, control_median, n_control, p, stars`, with stars
   `*` p < 0.05, `**` p < 0.01, `***` p < 0.001.

Throws `ArgumentError` for an unknown `alternative`, a missing column, an analysis without
exactly one case row or without controls, or a parameter with no finite control value.
"""
function empirical_test(
    tbl::AbstractDataFrame;
    omega::Bool = false,
    params::AbstractVector{Symbol} = _empirical_params(omega),
    alternative::Symbol = :greater,
)
    alternative in (:greater, :less, :two_sided) || throw(ArgumentError(
        "alternative must be :greater, :less or :two_sided, got :$alternative"))
    for c in vcat([:analysis, :type], params)
        hasproperty(tbl, c) || throw(ArgumentError(
            "column $c is missing" * (c in _empirical_params(true) ?
            " (abcmk drops the ω columns when any ω is NaN)" : "")))
    end

    out = DataFrame(analysis = String[], param = String[], case = Float64[],
                    control_median = Float64[], n_control = Int[], p = Float64[],
                    stars = String[])
    for a in unique(tbl.analysis)
        sub = tbl[tbl.analysis .== a, :]
        is_case = sub.type .== "case"
        count(is_case) == 1 || throw(ArgumentError(
            "analysis \"$a\" has $(count(is_case)) case rows; it needs exactly one"))
        any(sub.type .== "control") || throw(ArgumentError("analysis \"$a\" has no controls"))
        for c in params
            x = Float64(sub[is_case, c][1])
            ctrl = filter(isfinite, Float64.(sub[sub.type .== "control", c]))
            n = length(ctrl)
            n > 0 || throw(ArgumentError("analysis \"$a\": no finite control value for $c"))
            p_hi = (1 + count(>=(x), ctrl)) / (1 + n)
            p_lo = (1 + count(<=(x), ctrl)) / (1 + n)
            pv = alternative == :greater ? p_hi :
                 alternative == :less ? p_lo : min(1.0, 2 * min(p_hi, p_lo))
            push!(out, (string(a), string(c), x, median(ctrl), n, pv, _stars(pv)))
        end
    end
    return out
end

_stars(p::Real) = p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : ""

const hla_genes = ["ENSG00000233095","ENSG00000223980","ENSG00000230254","ENSG00000204642","ENSG00000204632","ENSG00000206503","ENSG00000204592","ENSG00000114455","ENSG00000235220","ENSG00000235346","ENSG00000235657","ENSG00000229252","ENSG00000228299","ENSG00000228987","ENSG00000236884","ENSG00000236418","ENSG00000233209","ENSG00000223793","ENSG00000228813","ENSG00000243496","ENSG00000226264","ENSG00000241394","ENSG00000235744","ENSG00000229685","ENSG00000237710","ENSG00000206435","ENSG00000228964","ENSG00000234794","ENSG00000227357","ENSG00000228080","ENSG00000232062","ENSG00000231939","ENSG00000231526","ENSG00000226165","ENSG00000239457","ENSG00000241296","ENSG00000243215","ENSG00000232957","ENSG00000236177","ENSG00000226826","ENSG00000237508","ENSG00000226260","ENSG00000227826","ENSG00000229074","ENSG00000225890","ENSG00000235680","ENSG00000225824","ENSG00000231834","ENSG00000231823","ENSG00000229493","ENSG00000239329","ENSG00000243189","ENSG00000231558","ENSG00000228163","ENSG00000236632","ENSG00000229295","ENSG00000237022","ENSG00000224608","ENSG00000229698","ENSG00000132297","ENSG00000206509","ENSG00000206506","ENSG00000206505","ENSG00000206493","ENSG00000206452","ENSG00000206450","ENSG00000204525","ENSG00000234745","ENSG00000232962","ENSG00000224103","ENSG00000230708","ENSG00000206308","ENSG00000196101","ENSG00000206306","ENSG00000206305","ENSG00000206302","ENSG00000225103","ENSG00000196610","ENSG00000241674","ENSG00000242685","ENSG00000206292","ENSG00000206291","ENSG00000215048","ENSG00000233841","ENSG00000223532","ENSG00000204287","ENSG00000198502","ENSG00000196126","ENSG00000196735","ENSG00000179344","ENSG00000237541","ENSG00000232629","ENSG00000241106","ENSG00000137403","ENSG00000230413","ENSG00000224320","ENSG00000197568","ENSG00000242574","ENSG00000204257","ENSG00000204252","ENSG00000231389","ENSG00000223865","ENSG00000233904","ENSG00000234487","ENSG00000237216","ENSG00000229215","ENSG00000227715","ENSG00000230726","ENSG00000231021","ENSG00000228284","ENSG00000231286","ENSG00000225201","ENSG00000206301","ENSG00000228254","ENSG00000241910","ENSG00000242386","ENSG00000239463","ENSG00000230141","ENSG00000168384","ENSG00000230763","ENSG00000230463","ENSG00000233192","ENSG00000230675","ENSG00000227993","ENSG00000243612","ENSG00000231679","ENSG00000206240","ENSG00000206237","ENSG00000204276","ENSG00000224305","ENSG00000241386","ENSG00000225691","ENSG00000242092","ENSG00000242361","ENSG00000235844","ENSG00000236693","ENSG00000232126","ENSG00000234154","ENSG00000243719"];

const hist_genes = ["ENSG00000164508","ENSG00000146047","ENSG00000124610","ENSG00000198366","ENSG00000196176","ENSG00000124529","ENSG00000124693","ENSG00000137259","ENSG00000196226","ENSG00000196532","ENSG00000187837","ENSG00000197061","ENSG00000187475","ENSG00000180596","ENSG00000180573","ENSG00000168298","ENSG00000158373","ENSG00000197697","ENSG00000188987","ENSG00000197409","ENSG00000196866","ENSG00000197846","ENSG00000198518","ENSG00000187990","ENSG00000168274","ENSG00000196966","ENSG00000124575","ENSG00000198327","ENSG00000124578","ENSG00000256316","ENSG00000197459","ENSG00000256018","ENSG00000168242","ENSG00000158406","ENSG00000124635","ENSG00000196787","ENSG00000197903","ENSG00000198339","ENSG00000184825","ENSG00000185130","ENSG00000196747","ENSG00000203813","ENSG00000182611","ENSG00000196374","ENSG00000197238","ENSG00000197914","ENSG00000184348","ENSG00000233822","ENSG00000198374","ENSG00000184357","ENSG00000182572","ENSG00000198558","ENSG00000197153","ENSG00000233224","ENSG00000196331","ENSG00000168148","ENSG00000181218","ENSG00000196890","ENSG00000203818","ENSG00000203814","ENSG00000183598","ENSG00000183941","ENSG00000203811","ENSG00000183558","ENSG00000203812","ENSG00000203852","ENSG00000182217","ENSG00000184678","ENSG00000184260","ENSG00000184270","ENSG00000197837","ENSG00000265232","ENSG00000263376","ENSG00000265133","ENSG00000265198","ENSG00000266225","ENSG00000266725","ENSG00000264719","ENSG00000263521","ENSG00000263965"];

# const hla_genes = ["ENSG00000153029","ENSG00000179583","ENSG00000221887","ENSG00000101294","ENSG00000158497","ENSG00000204642","ENSG00000204632","ENSG00000206503","ENSG00000204592","ENSG00000204525","ENSG00000234745","ENSG00000204287","ENSG00000198502","ENSG00000196126","ENSG00000196735","ENSG00000179344","ENSG00000237541","ENSG00000232629","ENSG00000241106","ENSG00000242574","ENSG00000204257","ENSG00000204252","ENSG00000231389","ENSG00000223865","ENSG00000277263","ENSG00000276051","ENSG00000229215","ENSG00000230463","ENSG00000233192","ENSG00000230675","ENSG00000243612","ENSG00000242092","ENSG00000242361","ENSG00000235844","ENSG00000236693","ENSG00000137403","ENSG00000230413","ENSG00000224320","ENSG00000233904","ENSG00000233841","ENSG00000223532","ENSG00000227993","ENSG00000231679","ENSG00000206240","ENSG00000257473","ENSG00000206237","ENSG00000224305","ENSG00000241386","ENSG00000234154","ENSG00000243719","ENSG00000232962","ENSG00000224103","ENSG00000230708","ENSG00000235220","ENSG00000235346","ENSG00000235657","ENSG00000229252","ENSG00000228299","ENSG00000228987","ENSG00000236884","ENSG00000236418","ENSG00000233209","ENSG00000223793","ENSG00000228813","ENSG00000243496","ENSG00000226264","ENSG00000241394","ENSG00000235744","ENSG00000229685","ENSG00000237710","ENSG00000237508","ENSG00000235680","ENSG00000231834","ENSG00000236632","ENSG00000237022","ENSG00000224608","ENSG00000226260","ENSG00000227826","ENSG00000229074","ENSG00000225890","ENSG00000225824","ENSG00000231823","ENSG00000229493","ENSG00000239329","ENSG00000243189","ENSG00000231558","ENSG00000228163","ENSG00000229295","ENSG00000234487","ENSG00000237216","ENSG00000227715","ENSG00000225201","ENSG00000225691","ENSG00000232126","ENSG00000230726","ENSG00000231021","ENSG00000228284","ENSG00000231286","ENSG00000206301","ENSG00000228254","ENSG00000241910","ENSG00000242386","ENSG00000239463","ENSG00000230141","ENSG00000168384","ENSG00000230763","ENSG00000206509","ENSG00000206506","ENSG00000206505","ENSG00000206493","ENSG00000206452","ENSG00000206450","ENSG00000206308","ENSG00000196101","ENSG00000206306","ENSG00000206305","ENSG00000206302","ENSG00000225103","ENSG00000196610","ENSG00000241674","ENSG00000242685","ENSG00000206292","ENSG00000206291","ENSG00000215048","ENSG00000229698","ENSG00000233095","ENSG00000223980","ENSG00000230254","ENSG00000206435","ENSG00000228964","ENSG00000234794","ENSG00000227357","ENSG00000228080","ENSG00000232062","ENSG00000231939","ENSG00000231526","ENSG00000226165","ENSG00000239457","ENSG00000241296","ENSG00000243215","ENSG00000232957","ENSG00000236177","ENSG00000226826"]

# const hist_genes = ["ENSG00000116478","ENSG00000273213","ENSG00000203814","ENSG00000183598","ENSG00000270882","ENSG00000203811","ENSG00000288825","ENSG00000288859","ENSG00000203852","ENSG00000270276","ENSG00000184678","ENSG00000184260","ENSG00000184270","ENSG00000143379","ENSG00000116539","ENSG00000117222","ENSG00000163041","ENSG00000168148","ENSG00000181218","ENSG00000196890","ENSG00000152455","ENSG00000078403","ENSG00000099284","ENSG00000149308","ENSG00000188486","ENSG00000172273","ENSG00000197837","ENSG00000246705","ENSG00000188375","ENSG00000061273","ENSG00000187166","ENSG00000139718","ENSG00000136169","ENSG00000100601","ENSG00000099381","ENSG00000177602","ENSG00000109111","ENSG00000290320","ENSG00000108840","ENSG00000132475","ENSG00000154655","ENSG00000104885","ENSG00000105011","ENSG00000162961","ENSG00000128708","ENSG00000068024","ENSG00000185513","ENSG00000274559","ENSG00000234289","ENSG00000099954","ENSG00000100084","ENSG00000189060","ENSG00000100395","ENSG00000100429","ENSG00000163517","ENSG00000181555","ENSG00000184897","ENSG00000178804","ENSG00000164032","ENSG00000145391","ENSG00000056050","ENSG00000268799","ENSG00000269466","ENSG00000113648","ENSG00000171720","ENSG00000164508","ENSG00000146047","ENSG00000124610","ENSG00000275714","ENSG00000278637","ENSG00000278705","ENSG00000286522","ENSG00000278463","ENSG00000276410","ENSG00000287080","ENSG00000187837","ENSG00000197061","ENSG00000187475","ENSG00000180596","ENSG00000180573","ENSG00000168298","ENSG00000158373","ENSG00000274290","ENSG00000277157","ENSG00000197409","ENSG00000196866","ENSG00000277224","ENSG00000276966","ENSG00000273802","ENSG00000277075","ENSG00000274750","ENSG00000124575","ENSG00000274618","ENSG00000275663","ENSG00000277775","ENSG00000275713","ENSG00000273983","ENSG00000278588","ENSG00000158406","ENSG00000124635","ENSG00000196787","ENSG00000276180","ENSG00000197903","ENSG00000274997","ENSG00000185130","ENSG00000196747","ENSG00000278828","ENSG00000276368","ENSG00000273703","ENSG00000197238","ENSG00000273542","ENSG00000275221","ENSG00000233822","ENSG00000276903","ENSG00000184357","ENSG00000275379","ENSG00000275126","ENSG00000197153","ENSG00000278677","ENSG00000274641","ENSG00000204371","ENSG00000196591","ENSG00000111875","ENSG00000198945","ENSG00000146414","ENSG00000048052","ENSG00000105968","ENSG00000285480","ENSG00000129691","ENSG00000181090","ENSG00000288210","ENSG00000285435","ENSG00000284841","ENSG00000285449","ENSG00000232045","ENSG00000236759","ENSG00000227333","ENSG00000238134","ENSG00000206376","ENSG00000224143","ENSG00000249467","ENSG00000187516","ENSG00000229674","ENSG00000101945","ENSG00000094631","ENSG00000147099","ENSG00000123569","ENSG00000101812","ENSG00000274183","ENSG00000277858","ENSG00000277745"]

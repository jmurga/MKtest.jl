"""
Mutable structure containing the variables required to boostrap a gene file following. All the functions are solve using the internal values of the structure. You should declare a mutable structure to the perform the analytical estimations.

# Bootstrap parameters
 - `data::String`: gene list get control genes.
 - `annotation::String`: bed file containing start and end gene coordinates.
 - `dist::Int64`: miminum distance between case and control genes.
 - `rep::Int64`: maximum number of genes by control set.
 - `tol::Float64`: 
 - `iter::Int64`: number of time defining the number of control sets. Each iter produce 100 sets.
 - `factors::String`: file containing confounding factors.
 - `output::String`

"""
@with_kw mutable struct bootstrap_parameters
    data::String = "ensembl_list.txt"
    annotation::String = "ensembl_gene_coords_v69.bed"
    dist::Int64 = 1
    rep::Int64 = 3
    tol::Float64 = 0.05
    iter::Int64 = 10
    factors::String = "confounding_factors.txt"
    output::String = "test"
end

"""
	Filtering case and intersect with control.
	The function will also filter by hla_dist and factors_raw (variables)
"""
function filter_case_control(x::Vector{String}, y::Vector{String}, distance::Vector{String})

    # setdiff!(x, hla_hist)
    # x_case::Vector{String} = intersect(x, y[:, 1])
    y_control::Vector{String} = setdiff(y[:, 1], x)

    intersect!(y_control, distance)

    return (string.(y_control))
end

"""
    Filtering case and intersect with control.
    The function will also filter by hla_dist and factors_raw (variables)
"""
function get_factors(
    case_set::Vector{String},
    control_set::Vector{String},
    factors_raw::Matrix,
)
    control_values = factors_raw[(@view factors_raw[:, 1]).∈[control_set], 1:end]
    control_dict = Dict{String,Vector{Float64}}()

    for row in eachrow(control_values)
        control_dict[row[1]] = row[2:end]
    end

    control_number::Dict{String,Int64} = Dict(control_set .=> 0)

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
    Check limits 
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

    df = Array(CSV.read(distance, DataFrame, delim = '\t', header = false))

    rm.([tmp, gene_1, gene_2])

    return df
end

function get_bootstrap(
    case_set::Vector{String},
    control_set::Vector{String};
    factors_raw::Matrix,
    max_reps::Int64,
    tolerance::Float64,
)

    # Get cont
    case_number = length(case_set)
    control_number = length(control_set)

    # Get averages
    case_avg, factors, factors_dict, used_number =
        get_factors(case_set, control_set, factors_raw)
    control_avg = vec(mean(factors, dims = 1))

    # high_low
    high_low, inf_sup = check_limits(case_avg, tolerance)
    flt_avg = case_avg .> control_avg
    high_low[flt_avg] .= "higher"
    high_low[.!flt_avg] .= "lower"

    # Increase control dataset to iter to create the control sets per se. Creating a big list speeds up the process because it avoids repeating the slow start of adding control genes to the "fake" seed.
    # gene_1, gene_2 = control_iteration(x, y)

    ######################
    # Iter combinations
    ######################
    previous_gn = 0
    same_counter = 0

    #Control sets are created by packs of 100 to speed up control creation while keeping memory usage in check.
    init_fake = fake_seed_num = 100

    # to keep the control sets close the averages for confounding factors in the group of interest, we later require that there some level of alternance with the newly added control genes swithcing between increasing the overall average and decreasing it.
    lim = size(factors, 2)

    direction = fill("start", lim)
    current_factors = case_avg

    added_nonvips = 0
    good_control = String[]

    sc_gn = length(good_control)

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

function bootstrap(param::bootstrap_parameters)
    @unpack data, annotation, dist, rep, tol, iter, factors, output = param

    @info "Opening data"

    case = vec(
        Array(CSV.read(data, header = false, delim = '\t', DataFrame, stringtype = String)),
    )
    factors_raw = Array(CSV.read(factors, header = false, DataFrame))
    factors_id = string.(@view factors_raw[:, 1])
    factors = Float64.(@view factors_raw[:, 2:end])

    @info "Estimating distance between genes"
    distance_raw = get_distance(data, annotation)

    #Filter case, control and distance
    distance = string.(distance_raw[distance_raw[:, 2].>=dist, 1])

    control = filter_case_control(case, factors_id, distance)

    case_number = length(case)
    control_number = length(control)

    @info "There are " *
          string(case_number) *
          " genes of interest and " *
          string(control_number) *
          " potential control genes at distance of at least " *
          string(dist) *
          " bases."

    if (control_number <= (1.5 * case_number))
        @warn "Warning! The number of control genes is less than 1.5 times the number of genes of interest. FDR may be high as a result."
    end

    # using ThreadsX
    n_case = fill(case, iter)
    n_control = fill(control, iter)

    @info "Running bootstrap"
    good_control = ThreadsX.map(
        (x, y) -> get_bootstrap(
            x,
            y,
            factors_raw = factors_raw,
            max_reps = rep,
            tolerance = tol,
        ),
        n_case,
        n_control,
    )

    if any(isempty.(good_control))
        @error "The bootstrap struggled too much to find matching controls. Stopping bootstrap test here. Please restart the entire pipeline for safety at input $data changing the parameters"
    else
        out_samples = ThreadsX.map(x -> get_samples(case, x), good_control)

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

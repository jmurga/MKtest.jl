"""
Mutable structure containing the variables required to boostrap a gene file following. All the functions are solve using the internal values of the structure. You should declare a mutable structure to the perform the analytical estimations.

# Bootstrap parameters
 - `data::String`: gene list get control genes.
 - `annotations::String`: bed file containing start and end gene coordinates.
 - `distance::Int64`: miminum distance between case and control genes.
 - `max_reps::Int64`: maximum number of genes by control set.
 - `tolerance::Float64`:
 - `iteration::Int64`: number of iteration defining the total number of control sets. Each iteration produce 100 control set.
 - `confounding_factors::String`: file containing confounding factors.
 - `output::String`: output file name without extension
 - `filter_hla_hist::Bool`: filter human HLA and HIST genes from Ensembl v109 annotations

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
    Filtering case and intersect with control.
    The function will also filter by hla_dist and factors_raw (variables)
"""
function get_valid(x::Vector{String15}, y::Vector{String15}, distance::Vector{String15})

    y_control = intersect(x,intersect(y, distance))

    return
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

    df = CSV.read(distance, DataFrame, delim = '\t', header = false)

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
        (x, y) -> get_bootstrap(
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

    empirical_pvalues(observed_values)

Estimating empirical p-values from a vector of summary statistics

# Arguments
 - `observed_values::Vector{Float64}` : Vector of summary statistics to convert into empirical p-values.
# Returns
 - Vector{Float64} : Empirical p-values. Large values will show small p-values.
"""
function empirical_pvalues(observed_values::Vector{Float64})
  return ordinalrank(observed_values * -1)/length(observed_values)
end

const hla_genes = ["ENSG00000233095","ENSG00000223980","ENSG00000230254","ENSG00000204642","ENSG00000204632","ENSG00000206503","ENSG00000204592","ENSG00000114455","ENSG00000235220","ENSG00000235346","ENSG00000235657","ENSG00000229252","ENSG00000228299","ENSG00000228987","ENSG00000236884","ENSG00000236418","ENSG00000233209","ENSG00000223793","ENSG00000228813","ENSG00000243496","ENSG00000226264","ENSG00000241394","ENSG00000235744","ENSG00000229685","ENSG00000237710","ENSG00000206435","ENSG00000228964","ENSG00000234794","ENSG00000227357","ENSG00000228080","ENSG00000232062","ENSG00000231939","ENSG00000231526","ENSG00000226165","ENSG00000239457","ENSG00000241296","ENSG00000243215","ENSG00000232957","ENSG00000236177","ENSG00000226826","ENSG00000237508","ENSG00000226260","ENSG00000227826","ENSG00000229074","ENSG00000225890","ENSG00000235680","ENSG00000225824","ENSG00000231834","ENSG00000231823","ENSG00000229493","ENSG00000239329","ENSG00000243189","ENSG00000231558","ENSG00000228163","ENSG00000236632","ENSG00000229295","ENSG00000237022","ENSG00000224608","ENSG00000229698","ENSG00000132297","ENSG00000206509","ENSG00000206506","ENSG00000206505","ENSG00000206493","ENSG00000206452","ENSG00000206450","ENSG00000204525","ENSG00000234745","ENSG00000232962","ENSG00000224103","ENSG00000230708","ENSG00000206308","ENSG00000196101","ENSG00000206306","ENSG00000206305","ENSG00000206302","ENSG00000225103","ENSG00000196610","ENSG00000241674","ENSG00000242685","ENSG00000206292","ENSG00000206291","ENSG00000215048","ENSG00000233841","ENSG00000223532","ENSG00000204287","ENSG00000198502","ENSG00000196126","ENSG00000196735","ENSG00000179344","ENSG00000237541","ENSG00000232629","ENSG00000241106","ENSG00000137403","ENSG00000230413","ENSG00000224320","ENSG00000197568","ENSG00000242574","ENSG00000204257","ENSG00000204252","ENSG00000231389","ENSG00000223865","ENSG00000233904","ENSG00000234487","ENSG00000237216","ENSG00000229215","ENSG00000227715","ENSG00000230726","ENSG00000231021","ENSG00000228284","ENSG00000231286","ENSG00000225201","ENSG00000206301","ENSG00000228254","ENSG00000241910","ENSG00000242386","ENSG00000239463","ENSG00000230141","ENSG00000168384","ENSG00000230763","ENSG00000230463","ENSG00000233192","ENSG00000230675","ENSG00000227993","ENSG00000243612","ENSG00000231679","ENSG00000206240","ENSG00000206237","ENSG00000204276","ENSG00000224305","ENSG00000241386","ENSG00000225691","ENSG00000242092","ENSG00000242361","ENSG00000235844","ENSG00000236693","ENSG00000232126","ENSG00000234154","ENSG00000243719"]

const hist_genes = ["ENSG00000164508","ENSG00000146047","ENSG00000124610","ENSG00000198366","ENSG00000196176","ENSG00000124529","ENSG00000124693","ENSG00000137259","ENSG00000196226","ENSG00000196532","ENSG00000187837","ENSG00000197061","ENSG00000187475","ENSG00000180596","ENSG00000180573","ENSG00000168298","ENSG00000158373","ENSG00000197697","ENSG00000188987","ENSG00000197409","ENSG00000196866","ENSG00000197846","ENSG00000198518","ENSG00000187990","ENSG00000168274","ENSG00000196966","ENSG00000124575","ENSG00000198327","ENSG00000124578","ENSG00000256316","ENSG00000197459","ENSG00000256018","ENSG00000168242","ENSG00000158406","ENSG00000124635","ENSG00000196787","ENSG00000197903","ENSG00000198339","ENSG00000184825","ENSG00000185130","ENSG00000196747","ENSG00000203813","ENSG00000182611","ENSG00000196374","ENSG00000197238","ENSG00000197914","ENSG00000184348","ENSG00000233822","ENSG00000198374","ENSG00000184357","ENSG00000182572","ENSG00000198558","ENSG00000197153","ENSG00000233224","ENSG00000196331","ENSG00000168148","ENSG00000181218","ENSG00000196890","ENSG00000203818","ENSG00000203814","ENSG00000183598","ENSG00000183941","ENSG00000203811","ENSG00000183558","ENSG00000203812","ENSG00000203852","ENSG00000182217","ENSG00000184678","ENSG00000184260","ENSG00000184270","ENSG00000197837","ENSG00000265232","ENSG00000263376","ENSG00000265133","ENSG00000265198","ENSG00000266225","ENSG00000266725","ENSG00000264719","ENSG00000263521","ENSG00000263965"]

# const hla_genes = ["ENSG00000153029","ENSG00000179583","ENSG00000221887","ENSG00000101294","ENSG00000158497","ENSG00000204642","ENSG00000204632","ENSG00000206503","ENSG00000204592","ENSG00000204525","ENSG00000234745","ENSG00000204287","ENSG00000198502","ENSG00000196126","ENSG00000196735","ENSG00000179344","ENSG00000237541","ENSG00000232629","ENSG00000241106","ENSG00000242574","ENSG00000204257","ENSG00000204252","ENSG00000231389","ENSG00000223865","ENSG00000277263","ENSG00000276051","ENSG00000229215","ENSG00000230463","ENSG00000233192","ENSG00000230675","ENSG00000243612","ENSG00000242092","ENSG00000242361","ENSG00000235844","ENSG00000236693","ENSG00000137403","ENSG00000230413","ENSG00000224320","ENSG00000233904","ENSG00000233841","ENSG00000223532","ENSG00000227993","ENSG00000231679","ENSG00000206240","ENSG00000257473","ENSG00000206237","ENSG00000224305","ENSG00000241386","ENSG00000234154","ENSG00000243719","ENSG00000232962","ENSG00000224103","ENSG00000230708","ENSG00000235220","ENSG00000235346","ENSG00000235657","ENSG00000229252","ENSG00000228299","ENSG00000228987","ENSG00000236884","ENSG00000236418","ENSG00000233209","ENSG00000223793","ENSG00000228813","ENSG00000243496","ENSG00000226264","ENSG00000241394","ENSG00000235744","ENSG00000229685","ENSG00000237710","ENSG00000237508","ENSG00000235680","ENSG00000231834","ENSG00000236632","ENSG00000237022","ENSG00000224608","ENSG00000226260","ENSG00000227826","ENSG00000229074","ENSG00000225890","ENSG00000225824","ENSG00000231823","ENSG00000229493","ENSG00000239329","ENSG00000243189","ENSG00000231558","ENSG00000228163","ENSG00000229295","ENSG00000234487","ENSG00000237216","ENSG00000227715","ENSG00000225201","ENSG00000225691","ENSG00000232126","ENSG00000230726","ENSG00000231021","ENSG00000228284","ENSG00000231286","ENSG00000206301","ENSG00000228254","ENSG00000241910","ENSG00000242386","ENSG00000239463","ENSG00000230141","ENSG00000168384","ENSG00000230763","ENSG00000206509","ENSG00000206506","ENSG00000206505","ENSG00000206493","ENSG00000206452","ENSG00000206450","ENSG00000206308","ENSG00000196101","ENSG00000206306","ENSG00000206305","ENSG00000206302","ENSG00000225103","ENSG00000196610","ENSG00000241674","ENSG00000242685","ENSG00000206292","ENSG00000206291","ENSG00000215048","ENSG00000229698","ENSG00000233095","ENSG00000223980","ENSG00000230254","ENSG00000206435","ENSG00000228964","ENSG00000234794","ENSG00000227357","ENSG00000228080","ENSG00000232062","ENSG00000231939","ENSG00000231526","ENSG00000226165","ENSG00000239457","ENSG00000241296","ENSG00000243215","ENSG00000232957","ENSG00000236177","ENSG00000226826"]

# const hist_genes = ["ENSG00000116478","ENSG00000273213","ENSG00000203814","ENSG00000183598","ENSG00000270882","ENSG00000203811","ENSG00000288825","ENSG00000288859","ENSG00000203852","ENSG00000270276","ENSG00000184678","ENSG00000184260","ENSG00000184270","ENSG00000143379","ENSG00000116539","ENSG00000117222","ENSG00000163041","ENSG00000168148","ENSG00000181218","ENSG00000196890","ENSG00000152455","ENSG00000078403","ENSG00000099284","ENSG00000149308","ENSG00000188486","ENSG00000172273","ENSG00000197837","ENSG00000246705","ENSG00000188375","ENSG00000061273","ENSG00000187166","ENSG00000139718","ENSG00000136169","ENSG00000100601","ENSG00000099381","ENSG00000177602","ENSG00000109111","ENSG00000290320","ENSG00000108840","ENSG00000132475","ENSG00000154655","ENSG00000104885","ENSG00000105011","ENSG00000162961","ENSG00000128708","ENSG00000068024","ENSG00000185513","ENSG00000274559","ENSG00000234289","ENSG00000099954","ENSG00000100084","ENSG00000189060","ENSG00000100395","ENSG00000100429","ENSG00000163517","ENSG00000181555","ENSG00000184897","ENSG00000178804","ENSG00000164032","ENSG00000145391","ENSG00000056050","ENSG00000268799","ENSG00000269466","ENSG00000113648","ENSG00000171720","ENSG00000164508","ENSG00000146047","ENSG00000124610","ENSG00000275714","ENSG00000278637","ENSG00000278705","ENSG00000286522","ENSG00000278463","ENSG00000276410","ENSG00000287080","ENSG00000187837","ENSG00000197061","ENSG00000187475","ENSG00000180596","ENSG00000180573","ENSG00000168298","ENSG00000158373","ENSG00000274290","ENSG00000277157","ENSG00000197409","ENSG00000196866","ENSG00000277224","ENSG00000276966","ENSG00000273802","ENSG00000277075","ENSG00000274750","ENSG00000124575","ENSG00000274618","ENSG00000275663","ENSG00000277775","ENSG00000275713","ENSG00000273983","ENSG00000278588","ENSG00000158406","ENSG00000124635","ENSG00000196787","ENSG00000276180","ENSG00000197903","ENSG00000274997","ENSG00000185130","ENSG00000196747","ENSG00000278828","ENSG00000276368","ENSG00000273703","ENSG00000197238","ENSG00000273542","ENSG00000275221","ENSG00000233822","ENSG00000276903","ENSG00000184357","ENSG00000275379","ENSG00000275126","ENSG00000197153","ENSG00000278677","ENSG00000274641","ENSG00000204371","ENSG00000196591","ENSG00000111875","ENSG00000198945","ENSG00000146414","ENSG00000048052","ENSG00000105968","ENSG00000285480","ENSG00000129691","ENSG00000181090","ENSG00000288210","ENSG00000285435","ENSG00000284841","ENSG00000285449","ENSG00000232045","ENSG00000236759","ENSG00000227333","ENSG00000238134","ENSG00000206376","ENSG00000224143","ENSG00000249467","ENSG00000187516","ENSG00000229674","ENSG00000101945","ENSG00000094631","ENSG00000147099","ENSG00000123569","ENSG00000101812","ENSG00000274183","ENSG00000277858","ENSG00000277745"]

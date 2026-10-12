
#################
# ABC functions #
#################

q_lower(x::Vector{Float64}) = quantile(x, [0.025])
q_upper(x::Vector{Float64}) = quantile(x, [0.975])

# A `.gz` name is written gzip-compressed.
function write_files(x::Matrix{Float64}, name::String)
    CSV.write(name, Tables.table(x), delim = '\t', header = false, compress = endswith(name, ".gz"))
end;

function write_files(x::SubArray, name::String)
    CSV.write(name, Tables.table(x), delim = '\t', header = false, compress = endswith(name, ".gz"))
end;

# Plain or gzipped text through `GZip.open`, closed afterwards.
function _with_gz(f, file::String)
    io = open(file)
    try
        return f(io)
    finally
        close(io)
    end
end

# A positive whole number given as any Real (1000, 1000.0, 1e3) as an Int.
function _as_count(x::Real, name::String)
    (isfinite(x) && isinteger(x) && x > 0) ||
        throw(ArgumentError("$name must be a positive whole number, got $x"))
    return Int(x)
end

function write_files(x::DataFrame, name::String, newfile::Bool)
    CSV.write(name, x, delim = '\t', header = false, append = newfile)
end;

# Index of the dataset a file belongs to: 10 for `summstat_10.txt.gz` or `out_10.0.post.gz`.
_dataset_index(file::String) = parse(Int, match(r"^[a-z]+_(\d+)", basename(file))[1])

# Files of `output_folder` starting with `prefix`, in dataset order. `readdir` sorts names
# as text (summstat_1, summstat_10, summstat_100, summstat_11, ...), which would return
# the posteriors in that order instead of the order of the input datasets.
function _dataset_files(output_folder::String, prefix::String)
    files = filter(x -> startswith(basename(x), prefix), readdir(output_folder, join = true))
    return sort(files, by = _dataset_index)
end

"""
	ABCreg(output_folder, S, P, tol, abcreg)

Performing ABC inference using ABCreg. Please, be sure the input folder contain the observed and summary statistics produc  ed by MKtest.summary_statistics().

# Arguments
 - `output_folder::String` : folder containing the observed data and summary estatistics. It will be used to output the posterior distributions
 - `P::Int64` : number of parameters to infer.
 - `S::Int64` : number of summary stastitics to perform the inference.
 - `tol::Real` : number of accepted values at ABC inference, a positive whole number (`1000`, `1e3`).
 - `abcreg::String` : path to ABCreg binary.
# Output
 - `posteriors:Vector{Matrix{Float64}}`: posteriors distributions estimated by ABCreg.
"""
function ABCreg(;
    output_folder::String,
    S::Int64,
    P::Int64 = 12,
    tol::Real,
    rm_summaries::Bool = false,
    abcreg::Union{Nothing,String}=nothing,
    T::Bool=false,
    L::Bool=false
)

    @assert (T == true & L != true)  |  (T != true & L == true)  "Cannot perform both transformation at the same time"

    if isnothing(abcreg)
        abcreg = abspath(joinpath(@__DIR__, "..", "scripts", "reg"))
    end
     
    # List alphas and summstat files
    a_file   = _dataset_files(output_folder, "alphas_")
    sum_file = _dataset_files(output_folder, "summstat_")

    tol = _as_count(tol, "tol")

    # Change P if omega was not estimated
    if size(_with_gz(readdlm, sum_file[1]), 2) != (P+S)
        P = P - 4
    end

    # Creating output names
    out = output_folder .* "/out_" .* string.(_dataset_index.(a_file))
    tol_abcreg = tol ./ _with_gz.(countlines, sum_file)
    @info "Running ABCreg. Please cite https://doi.org/10.1186/1471-2156-10-35"
    # Using mapi instead map to limit tasks. Bash interaction not working as expected
    ThreadsX.mapi(
        (x, y, z, t) -> run(`$abcreg -d $x -p $y -P $P -S $S -t $t -b $z`),
        a_file,
        sum_file,
        out,
        tol_abcreg,
        ntasks = Threads.nthreads(),
    )

    @info "Opening and filtering posteriors distributions"
    out = filter(x -> occursin("post", x), _dataset_files(output_folder, "out_"))
    out = filter(x -> !occursin(".1.", x), out)

    # Move from CSV.read to readdlm, bug when threads enable in node server.
    # open(x) = Array(CSV.read(x, DataFrame))
    # Control outlier inference. 2Nes non negative values
    posteriors = read_gz(out)
    # Remove if some file was not computed. The rows after it then no longer match the
    # dataset numbers.
    ok = @. !iszero(posteriors)
    all(ok) || @warn "ABCreg failed for datasets $(_dataset_index.(out[.!ok])); their posteriors are dropped"
    posteriors = posteriors[ok]

    if P != 12
        c_names = [:α_weak, :α_strong, :α, :γ₋, :β,:γ₊,:γ₊₊,:B]
    else
        c_names = [:α_weak, :α_strong, :α,:ωₐ_weak,:ωₐ_strong,:ωₐ,:ωₙₐ,:γ₋,:β,:γ₊,:γ₊₊,:B]
    end

    posteriors = map(x -> DataFrame(x,c_names),posteriors)

    # Remove summstat files
    if rm_summaries
        rm.(
            filter(
                x -> occursin("summstat", x) || occursin("alphas_", x) || occursin("out_", x),
                readdir(output_folder, join = true),
            )
        )
    end

    return posteriors
end

function get_mode(posterior::Matrix{Float64})

    # Allocating outputs
    out = zeros((1, size(posterior, 2)))

    for j::Int64 = 1:size(posterior, 2)
        y = kde(posterior[:, j])
        m = collect(y.x)[y.density.==maximum(y.density)][1]
        out[j] = m
    end

    return (out)
end

function get_mode(posterior::DataFrame)

    # Allocating outputs
    out = zeros((1, size(posterior, 2)))

    for j::Int64 = 1:size(posterior, 2)
        y = kde(posterior[:, j])
        m = collect(y.x)[y.density.==maximum(y.density)][1]
        out[j] = m
    end

    return (out)
end

# read gz since CSV.read is bugged at node server
function read_gz(out_file::String)
    try
        return readdlm(open(out_file),'\t', Float64,header = false)
    catch
        @warn "$out_file is empty"
        return zeros(1, 12)
    end
end

function read_gz(out_file::Vector{String})
    return ThreadsX.mapi(x->read_gz(x),out_file)
end


"""
	summary_abc(posteriors; stat)

Posterior statistics: Min., 2.5% Perc., Median, Mean, Mode, 97.5% Perc. and Max. of every
parameter. The mode follows the R package abc (maximum of a kernel density estimate).

# Arguments
 - `posteriors::Vector{DataFrame}` or `Vector{Matrix{Float64}}`: one posterior per dataset,
   with 8 (no ω) or 12 columns.
 - `stat::String`: statistic returned as the point estimate: "Min.", "Median", "Mean",
   "Mode" or "Max" (case-insensitive). Default "Mode".
# Output
 - `OrderedDict` with
   - `:inference`: one row per dataset with the chosen statistic.
   - `:quantiles`: the same values as strings with the 95% interval, `"0.2 [0.1-0.3]"`.
   - `:stats`: all statistics, a `Stats` column then one row per statistic and dataset.
"""
function summary_abc(posteriors::Vector{DataFrame}; stat::String = "Mode")
    @info "Computing statistics over posteriors distributions"

    labels = ["Min.", "2.5% Perc.", "Median", "Mean", "Mode", "97.5% Perc.", "Max"]
    key = uppercasefirst(lowercase(stat))
    key = key == "Min" ? "Min." : key
    key in labels || throw(ArgumentError(
        "stat must be one of Min., Median, Mean, Mode or Max; got \"$stat\""))

    c_names = propertynames(first(posteriors))

    # One 7 x P block per dataset, rows in the order of `labels`.
    blocks = ThreadsX.map(posteriors) do x
        cols = eachcol(x)
        vcat(
            minimum.(cols)',
            (only.(q_lower.(Vector{Float64}.(cols))))',
            median.(cols)',
            mean.(cols)',
            get_mode(x),
            (only.(q_upper.(Vector{Float64}.(cols))))',
            maximum.(cols)',
        )
    end

    df = vcat([hcat(DataFrame(:Stats => labels), DataFrame(b, c_names)) for b in blocks]...)

    i_stat = findfirst(==(key), labels)
    inference = DataFrame(reduce(vcat, (b[i_stat:i_stat, :] for b in blocks)), c_names)
    quantiles = DataFrame(
        reduce(vcat, (permutedims(string.(round.(b[i_stat, :], digits = 3), " [",
                                          round.(b[2, :], digits = 3), "-",
                                          round.(b[6, :], digits = 3), "]"))
                      for b in blocks)),
        c_names,
    )

    length(posteriors) == 1 && @info df

    return OrderedDict(:inference => inference, :quantiles => quantiles, :stats => df)
end

function summary_abc(posteriors::Vector{Matrix{Float64}}; stat::String = "Mode")
    c_names(x) = size(x, 2) == 12 ?
        [:α_weak, :α_strong, :α, :ωₐ_weak, :ωₐ_strong, :ωₐ, :ωₙₐ, :γ₋, :β, :γ₊, :γ₊₊, :B] :
        [:α_weak, :α_strong, :α, :γ₋, :β, :γ₊, :γ₊₊, :B]
    return summary_abc([DataFrame(x, c_names(x)) for x in posteriors]; stat = stat)
end

#############


function abc(; output_folder::String,
    S::Int64,
    P::Int64 = 12,
    tol::Real,
    transformation::String = "none",
    kernel::String = "epanechnikov",
    rm_summaries::Bool = false
)

    # Change BLAS threads to 1 to proper parallel.
    nthreads_og = get_num_threads()
    set_num_threads(1)

    @assert lowercase(kernel) ∈ ["gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine"] "Kernel is incorrectly defined. Use gaussian, epanechnikov, rectangular, triangular, biweight or cosine"
    @assert lowercase(transformation) ∈ ["none", "log", "tan"] "Apply one of the following transformations: none, log, tan"

    # List/open alphas and summstat files
    targets = _dataset_files(output_folder, "alphas_")
    param_summaries = _dataset_files(output_folder, "summstat_")

    tol = _as_count(tol, "tol")

    # Change P if omega was not estimated
    if size(CSV.read(param_summaries[1],matrix,header=false), 2) != (P + S)
        P -= 4
    end

    @info "Running ABC"
    # Using mapi instead of map to limit tasks. Bash interaction not working as expected
    posteriors, posteriors_adjusted = unzip(ThreadsX.map((x, y) -> abc_loclinear(x,y, P=P, tol=tol, transformation=transformation, kernel=kernel), targets, param_summaries))

    @info "Filtering posterior distributions"
    # Remove if some file wasn't computed
    posteriors = filter(!iszero, posteriors)
    posteriors_adjusted = filter(!iszero, posteriors_adjusted)

    # Remove summstat files
    if rm_summaries
        rm.(filter(x -> occursin("summstat", x) || occursin("alphas_", x), readdir(output_folder, join = true)))
    end

    set_num_threads(nthreads_og)

    if P != 12
        c_names = [:α_weak, :α_strong, :α, :γ₋, :β,:γ₊,:γ₊₊,:B]
    else
        c_names = [:α_weak, :α_strong, :α,:ωₐ_weak,:ωₐ_strong,:ωₐ,:ωₙₐ,:γ₋,:β,:γ₊,:γ₊₊,:B]
    end

    posteriors = map(x -> DataFrame(x,c_names),posteriors)
    posteriors_adjusted = map(x -> DataFrame(x,c_names),posteriors_adjusted)

    return OrderedDict("rejection"=>posteriors, "loclinear"=>posteriors_adjusted)
end

function normalise(x::Union{SubArray,Float64},y::SubArray)
    if mad(y) != 0
        return x/mad(y)
    else
        return x
    end
end

function rejection(target::Vector{Float64},sumstat::Matrix{Float64},param::Matrix{Float64},tol::Float64,kernel::String)

    # Scale everything and euclidean distance
    sumstat_scaled = similar(sumstat)
    target_scaled  = similar(target)
    dist           = zeros(size(sumstat,1))

    for j in eachindex(target)
        sumstat_scaled[:,j] = normalise(view(sumstat,:,j),view(sumstat,:,j))
        target_scaled[j]    = normalise(target[j],view(sumstat,:,j))
        dist .+= (view(sumstat_scaled,:,j) .- target_scaled[j]).^2
    end
    dist = @. sqrt(dist)

    # Sort and get minimum distance to return values inside tolerance range
    n_accept = ceil(Int,length(dist)*tol)
    n_limit  = sort(dist)[n_accept]
    # Ensure get only n_limit, if more than one
    n_idx = findall(dist .<= n_limit)[1:n_accept]

    ## Weigthened distances
    wts = if(kernel == "epanechnikov")
        @.  1 - (dist[n_idx]/n_limit)^2
    elseif(kernel == "rectangular")
        @.  dist[n_idx]/n_limit
    elseif(kernel == "gaussian")
        @.  1/sqrt(2*pi)*exp(-0.5*(dist/(ds/2))^2)
    elseif(kernel == "triangular")
        @.  1 - abs(dist[n_idx]/n_limit)
    elseif(kernel == "biweight")
        @. (1 - (dist[n_idx]/n_limit)^2)^2
    elseif(kernel == "cosine")
        @. cos(pi/2*dist[n_idx]/n_limit)
    end


    return(target_scaled,sumstat[n_idx,:],sumstat_scaled[n_idx,:],param[n_idx,:],dist[n_idx],wts)
end

function regression(sumstat::Matrix{Float64},param::Matrix{Float64},weights::Vector{Float64})

    # Allocating results
    lm_residuals = zeros(size(sumstat,1),size(param,2))
    lm_coefficients = zeros(size(param,2),size(sumstat,2))

    # Loclinear regression
    for (i,v) in enumerate(eachcol(param))
        model = lm(sumstat,v,wts=weights)
        lm_coefficients[i,:] .= vec(coef(model))
        lm_residuals[:,i] .= vec(residuals(model))
    end

    return(lm_coefficients,lm_residuals)
end

function tangent_transfromation(x::Float64, minval::Float64, maxval::Float64)
    # Apply the tangent transformation of Hamilton et al. 2005 PNAS 7476
    return -log(1.0 / tan(((x - (minval - 1e-04)) / ((maxval + 1e-04) - (minval - 1e-04))) * (π / 2)))
end

function undo_tangent_transfromation(x::Float64, minval::Float64, maxval::Float64)
    # Undo the Hamilton tangent transformation
    return (minval - 1e-4) + (2 / π) * ((maxval + 1e-4) - (minval - 1e-4)) * atan(exp(y))
end

function abc_loclinear(target_file::String,param_summaries_file::String;P::Int64,tol::Int64,transformation::String="none",kernel::String= "epanechnikov")

    # Reading into matrix using Tables.matrix
    target = vec(CSV.read(target_file,matrix,ntasks=1,header=false))
    summaries = CSV.read(param_summaries_file,matrix,ntasks=1,header=false)
    # target = vec(readdlm(target_file,'\t',Float64,header=false))
    # summaries = readdlm(param_summaries_file,'\t',Float64,header=false))

    param   = summaries[:, 1:P]
    sumstat = summaries[:, (P+1):end]

    # transformation="none";kernel="epanechnikov";tol=0.025

    @assert lowercase(kernel) ∈ ["gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine"] "Kernel is incorrectly defined. Use gaussian, epanechnikov, rectangular, triangular, biweight or cosine"
    @assert lowercase(transformation) ∈ ["none", "log","tan"] "Apply one of the following transformations: none, log, tan"

    @assert length(target) == size(sumstat,2) "Number of summary statistics in target has to be the same as in sumstat."

    num_params = size(param,2)
    num_stats = size(sumstat,2)

    mins = minimum(param,dims=1)
    maxs = maximum(param,dims=1)

    ## Scale and euclidean distance
    ## #########################
    target_scaled, sumstat_accepted, sumstat_scaled, param_accepted, dist_accepted, wts = rejection(vec(target),sumstat,param,tol/size(summaries,1),kernel)
    num_accepted = size(sumstat_accepted,1)

    ## Transform parameters
    ## ######################
    if any(param .<= 0 .&& transformation == "log")
        @warn "log transform: values out of bounds - correcting..."
    end

    for (i,v) ∈ enumerate(eachcol(param_accepted))
        if transformation == "none"
            continue
        elseif transformation == "log"
            if(minimum(v) <= 0)
                v[findall(iszero, v)] .= minimum(filter(!iszero,v))
                @show minimum(filter(!iszero,v))
            end
            @. v .= log(v)
        elseif transformation == "tan"
            @. v = tangent_transfromation(v,mins[i],maxs[i])
        end
    end

    sumstat_intercept = hcat(ones(size(sumstat_scaled,1)),sumstat_scaled)

    # Linear regression
    lm_coefficients, lm_residuals = regression(sumstat_intercept,param_accepted,wts)


    pred = lm_coefficients * vcat(1,target_scaled)
    pred = repeat(pred',num_accepted)

    # rsdl = param_accepted .- (lm_coefficients * sumstat_intercept')'
    # rsdl = @view lm_residuals

    rsdl_mean = mapslices(mean,lm_residuals,dims=1)
    rsdl_corrected = lm_residuals .- rsdl_mean

    pred_corrected = pred .+ rsdl_mean

    f(x::Vector,wts::Vector{Float64}=wts) = sum((x).^2 .* wts)/sum(wts)
    σ = mapslices(f,rsdl_corrected,dims=1)
    aic = num_accepted * sum(log.(σ)) + 2*(num_stats+1)*num_params
    bic = num_accepted * sum(log.(σ)) + log(sum(num_accepted))*(num_stats+1)*num_params

    # Heteroscedasticity correction
    rsdl_log = @. log((lm_residuals)^2)
    lm_coefficients, lm_residuals = regression(sumstat_intercept,rsdl_log,wts)

    pred_sd = lm_coefficients * vcat(1,target_scaled)
    pred_sd = @. sqrt(exp(pred_sd))
    pred_sd = repeat(pred_sd',num_accepted)
    pred_si = (lm_coefficients * sumstat_intercept')'
    pred_si = @. sqrt(exp(pred_si))

    param_adjusted = @. pred + (pred_sd*rsdl_corrected) / pred_si
    rsdl_adjusted = @. (pred_sd*rsdl_corrected) / pred_si

    # Back transform parameter values
    for i=1:num_params
        if(transformation == "log")
            param_accepted[:,i] .= exp(param_accepted[:,i])
            param_adjusted[:,i] .= exp(param_adjusted[:,i])
        elseif(transformation == "tan")
            param_accepted[:,i] .= undo_tangent_transfromation(param_accepted[:,i],mins[i],maxs[i])
            param_adjusted[:,i] .= undo_tangent_transfromation(param_adjusted[:,i],mins[i],maxs[i])
        end
    end

    return param_accepted,param_adjusted
end

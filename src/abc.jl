
#################
# ABC functions #
#################

q_lower(x::Vector{Float64}) = quantile(x, [0.025])
q_upper(x::Vector{Float64}) = quantile(x, [0.975])

function write_files(x::Matrix{Float64}, name::String)
    CSV.write(name, Tables.table(x), delim = '\t', header = false)
end;

function write_files(x::SubArray, name::String)
    CSV.write(name, Tables.table(x), delim = '\t', header = false)
end;

function write_files(x::DataFrame, name::String, newfile::Bool)
    CSV.write(name, x, delim = '\t', header = false, append = newfile)
end;

"""
	ABCreg(output_folder, S, P, tol, abcreg)

Performing ABC inference using ABCreg. Please, be sure the input folder contain the observed and summary statistics produc  ed by MKtest.summary_statistics().

# Arguments
 - `output_folder::String` : folder containing the observed data and summary estatistics. It will be used to output the posterior distributions
 - `P::Int64` : number of parameters to infer.
 - `S::Int64` : number of summary stastitics to perform the inference.
 - `tol::Float64` : tolerance value. It defines the number of accepted values at ABC inference.
 - `abcreg::String` : path to ABCreg binary.
# Output
 - `posteriors:Vector{Matrix{Float64}}`: posteriors distributions estimated by ABCreg.
"""
function ABCreg(;
    output_folder::String,
    S::Int64,
    P::Int64 = 12,
    tol::Int64,
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
    a_file   = filter(x -> occursin("alphas", x), readdir(output_folder, join = true))
    sum_file = filter(x -> occursin("summstat", x), readdir(output_folder, join = true))

    # Change P if omega was not estimated
    if size(readdlm(sum_file[1]),2) != (P+S)
        P = P - 4
    end

    # Creating output names
    out = output_folder .* "/out_" .* sort!(string.(1:size(a_file, 1)))
    tol_abcreg = tol./countlines.(sum_file)
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
    out = filter(x -> occursin("post", x), readdir(output_folder, join = true))
    out = filter(x -> !occursin(".1.", x), out)

    # Move from CSV.read to readdlm, bug when threads enable in node server.
    # open(x) = Array(CSV.read(x, DataFrame))
    # Control outlier inference. 2Nes non negative values
    posteriors = read_gz(out)
    # Remove is some file wasnt computed
    posteriors = posteriors[@. !iszero(posteriors)]

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
	summary_abc(posteriors,stat)
Posterior distributions statistics. The function estimates Min., 0.5% Perc., Median, Mean, Mode, "99.5% Perc., Max values. The mode is estimated following R package abc mode function. The argument `stat` define the output estimation.

# Arguments
 - `posterior::Matrix{Float64}` : posterior distribution.
 - `stat::String`: output estimation.
# Output
 - `DataFrame`: posterior statistics.
 - `DataFrame`: chosen statistic inference.
"""
function summary_abc(posteriors::Vector{Matrix{Float64}}; stat::String = "Mode")
    @info "Computing statistics over posteriors distributions"

    p_min    = vcat(ThreadsX.map(x -> minimum(x, dims = 1), posteriors)...)
    p_max    = vcat(ThreadsX.map(x -> maximum(x, dims = 1), posteriors)...)
    p_mean   = vcat(ThreadsX.map(x -> mean(x, dims = 1), posteriors)...)
    p_median = vcat(ThreadsX.map(x -> median(x, dims = 1), posteriors)...)
    p_mode   = vcat(ThreadsX.map(get_mode, posteriors)...)
    p_lower  = vcat(ThreadsX.map(x -> mapslices(q_lower, x, dims = 1), posteriors)...)
    p_upper  = vcat(ThreadsX.map(x -> mapslices(q_upper, x, dims = 1), posteriors)...)

    if size(p_mean,2) != 12
        c_names = [:α_weak, :α_strong, :α, :γ₋, :β,:γ₊,:γ₊₊,:B]
    else
        c_names = [:α_weak, :α_strong, :α,:ωₐ_weak,:ωₐ_strong,:ωₐ,:ωₙₐ,:γ₋,:β,:γ₊,:γ₊₊,:B]
    end
    if length(posteriors) == 1

        quantile_string = String[]
        for i in zip(p_lower,p_upper)
            push!(quantile_string," [" * join(string.(round.(i,digits=3)),"-") * "]")
        end
        quantile_string = DataFrame(x1=vcat(quantile_string))

        stats = DataFrame(
            :Stats => repeat(
                ["Min.", "2.5% Perc.", "Median", "Mean", "Mode", "97.5% Perc.", "Max"],
                inner = length(posteriors),
            ),
        )

        tmp = DataFrame(
            vcat(p_min, p_lower, p_median, p_mean, p_mode, p_upper, p_max),
	    c_names
        )

        df = hcat(stats, tmp)

        @info df

        df_quantiles = permutedims(string.(permutedims(string.(round.(df[df.Stats.==uppercasefirst(stat), 2:end],digits=3))),quantile_string))
        rename!(df_quantiles,c_names)

        # return df[df.Stats.==uppercasefirst(stat), 2:end], df_quantiles, df
        return OrderedDict(:inference=>df[df.Stats.==uppercasefirst(stat), 2:end],:quantiles=>df_quantiles,:stats=>df)

    else
        df = DataFrame[]
        for i::Int64 = 1:length(posteriors)
            stats = DataFrame(
                :Stats => [
                    "Min.",
                    "2.5% Perc.",
                    "Median",
                    "Mean",
                    "Mode",
                    "97.5% Perc.",
                    "Max",
                ],
            )

            tmp = DataFrame(
                hcat(
                    p_min[i, :],
                    p_lower[i, :],
                    p_median[i, :],
                    p_mean[i, :],
                    p_mode[i, :],
                    p_upper[i, :],
                    p_max[i, :],
                )',
                c_names
            )
            push!(df, hcat(stats, tmp))
        end

        df_vcat = vcat(df...)

        # return  df_vcat[df_vcat.Stats.==uppercasefirst(stat), 2:end], df
        return OrderedDict(:inference=>df[df.Stats.==uppercasefirst(stat), 2:end],:quantiles=> df_quantiles,:stats=>df)

    end
end

function summary_abc(posteriors::Vector{DataFrame}; stat::String = "Mode")
    @info "Computing statistics over posteriors distributions"

    p_min    = vcat(ThreadsX.map(x -> minimum.(eachcol(x))', posteriors)...)
    p_max    = vcat(ThreadsX.map(x -> maximum.(eachcol(x))', posteriors)...)
    p_mean   = vcat(ThreadsX.map(x -> mean.(eachcol(x))', posteriors)...)
    p_median = vcat(ThreadsX.map(x -> median.(eachcol(x))', posteriors)...)
    p_mode   = vcat(ThreadsX.map(get_mode, posteriors)...)
    p_lower  = vcat(ThreadsX.map(x -> hcat(q_lower.(eachcol(x))...), posteriors)...)
    p_upper  = vcat(ThreadsX.map(x -> hcat(q_upper.(eachcol(x))...), posteriors)...)

    if size(p_mean,2) != 12
        c_names = [:α_weak, :α_strong, :α, :γ₋, :β,:γ₊,:γ₊₊,:B]
    else
        c_names = [:α_weak, :α_strong, :α,:ωₐ_weak,:ωₐ_strong,:ωₐ,:ωₙₐ,:γ₋,:β,:γ₊,:γ₊₊,:B]
    end

    if length(posteriors) == 1

        quantile_string = String[]
        for i in zip(p_lower,p_upper)
            push!(quantile_string," [" * join(string.(round.(i,digits=3)),"-") * "]")
        end
        quantile_string = DataFrame(x1=vcat(quantile_string))

        stats = DataFrame(
            :Stats => repeat(
                ["Min.", "5% Perc.", "Median", "Mean", "Mode", "95% Perc.", "Max"],
                inner = length(posteriors),
            ),
        )

        tmp = DataFrame(
            vcat(p_min, p_lower, p_median, p_mean, p_mode, p_upper, p_max),
        c_names
        )

        df = hcat(stats, tmp)

        # @info df

        df_quantiles = permutedims(string.(permutedims(string.(round.(df[df.Stats.==uppercasefirst(stat), 2:end],digits=3))),quantile_string))
        rename!(df_quantiles,c_names)

        return OrderedDict(:inference=>df[df.Stats.==uppercasefirst(stat), 2:end],:quantiles=>df_quantiles,:stats=>df)
        # return df[df.Stats.==uppercasefirst(stat), 2:end], df_quantiles, df

    else
        df = DataFrame[]
        for i::Int64 = 1:length(posteriors)
            stats = DataFrame(
                :Stats => [
                    "Min.",
                    "5% Perc.",
                    "Median",
                    "Mean",
                    "Mode",
                    "95% Perc.",
                    "Max",
                ],
            )

            tmp = DataFrame(
                hcat(
                    p_min[i, :],
                    p_lower[i, :],
                    p_median[i, :],
                    p_mean[i, :],
                    p_mode[i, :],
                    p_upper[i, :],
                    p_max[i, :],
                )',c_names
            )
            push!(df, hcat(stats, tmp))
        end

        df_vcat = vcat(df...)

        return OrderedDict(:inference=>df_vcat[df_vcat.Stats.==uppercasefirst(stat), 2:end],:stats=>df_vcat)
        # return df_vcat[df_vcat.Stats.==uppercasefirst(stat), 2:end], df
    end
end

#############


function abc(; output_folder::String,
    S::Int64,
    P::Int64 = 12,
    tol::Int64,
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
    targets = filter(x -> occursin("alphas", x), readdir(output_folder, join = true))
    param_summaries = filter(x -> occursin("summstat", x), readdir(output_folder, join = true))

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

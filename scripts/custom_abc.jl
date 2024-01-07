using StatsBase,ThreadsX,JLD2,SharedArrays

import LinearAlgebra.BLAS: get_num_threads, set_num_threads

function sampling_summaries(
    models::SharedMatrix,
    fs::Vector{Float64},
    d::Vector{Float64},
    neut::SharedMatrix,
    sel::SharedMatrix,
    dsdn::SharedMatrix,
)
    ds = @view dsdn[:, 1]
    dn = @view dsdn[:, 2]
    dweak = @view dsdn[:, 3]
    dstrong = @view dsdn[:, 4]
    gn = abs.(@view models[:, end])
    sh = round.(view(models, :, size(models, 2) - 1), digits = 5)
    gL = ceil.(abs.(@view models[:, end-3]))
    gH = ceil.(abs.(@view models[:, end-2]))
    B  = @view models[:,1]

    ## Outputs
    α, ω, expected_dn, expected_ds = MKtest.poisson_fixation(d, ds, dn, dweak, dstrong)

    expected_pn, expected_ps = MKtest.poisson_polymorphism(fs, permutedims(neut), permutedims(sel))

    ## Alpha from expected values. Used as summary statistics
    α_summaries = @. round(1 - ((expected_ds / expected_dn) * (expected_pn / expected_ps)'), digits = 5)
    if any(isnan.(ω))
        expected_values = hcat(round.(α, digits = 5), gn, sh, gH, gL, B, α_summaries)
    else
        expected_values = hcat(round.(α, digits = 5),round.(ω, digits = 5), gn, sh, gL, gH, B, α_summaries)
    end

    expected_values = MKtest.filter_expected(expected_values)

    return(expected_values)
end

function summary_statistics(
    param::MKtest.parameters,
    sfs::Vector,
    divergence::Vector;
    h5_file::String,
    summstat_size::Int64,
    output_folder::String,
    alpha::Union{Nothing,Vector{Float64}} = nothing,
    B_bins::Union{Nothing,Vector{Float64}} = nothing
)

    ## Opening files
    MKtest.assertion_params(param)

    α, sfs_p, divergence_p = MKtest.data_to_poisson(sfs, divergence, param.dac)

    # Stop execution if any of the SFS contains 0 values
    # Low polymorphism distribution will cause ABC errors due to Inf values on α estimations
    if any(0 .∈ sfs_p) | any(0 .∈ getindex.(divergence_p,1))

        throw(
            ArgumentError(
                "Your SFS contains 0 values at the selected DACs or the divergence is 0. Please consider to bin the SFS and re-estimate the rates using the selected bin as sample the new sample size.",
            ),
        )
    end

    # Open rates
    @info "Opening and random sampling $summstat_size expected SFS and fixation rates"
    string_cutoff =
        "cutoff=[" * string(param.cutoff[1]) * "," * string(param.cutoff[end]) * "]"
    h = jldopen(h5_file)
    tmp = h[string(param.N)*"/"*string(param.n)*"/"*string_cutoff]

    # Sample index to subset models
    idx = trues(size(tmp["models"],1))
    idx_alpha = trues(length(idx))
    idx_B = trues(length(idx))

    if !isnothing(alpha)
        idx_alpha = tmp["models"].al_tot .>= alpha[1] .&& tmp["models"].al_tot .<= alpha[end]
    end

    if !isnothing(B_bins)
        idx_B = tmp["models"].B .>= B_bins[1].*1000 .&& tmp["models"].B .<= B_bins[end].*1000
    end

    idx = findall(idx .&& idx_alpha .&& idx_B)

    @assert length(idx) >= summstat_size "Check the filters or increase the number of rates solutions"

    sampled_idx = sample(idx,summstat_size;replace=true)

    @assert length(sampled_idx) >= summstat_size "Check the filters or increase the number of rates solutions"

    # Convert random models to solve to a SharedMatrix. Only reading shouldn't generate race-condition
    models = SharedMatrix(Array(view(tmp["models"], sampled_idx, :)))
    dsdn   = SharedMatrix(tmp["dsdn"][sampled_idx, :])

    # Filtering polymorphic rate by dac
    neut = SharedMatrix(zeros(summstat_size,length(param.dac)))
    sel = SharedMatrix(zeros(summstat_size,length(param.dac)))

    for (i,v) in enumerate(param.dac)
        neut[:,i] .= view(tmp["neut"][v],sampled_idx,:)
        sel[:,i] .= view(tmp["sel"][v],sampled_idx,:)
    end

    # Making summaries

    @info "Sampling and computing summary statistics"

    # 20 threads: case + control ~ 25GB
    expected_values = ThreadsX.map(
        (x, y) -> sampling_summaries(models, x, y, neut, sel, dsdn),
        sfs_p,
        divergence_p
    );

    return(α,expected_values)
end

function abc(targets,summstat;
    S::Int64,
    P::Int64 = 12,
    tol::Float64,
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
    # targets = filter(x -> occursin("alphas", x), readdir(output_folder, join = true))
    # summstat = filter(x -> occursin("summstat", x), readdir(output_folder, join = true))

    # Change P if omega was not estimated
    if size(summstat[1], 2) != (P + S)
        P -= 4
    end

    @info "Running ABC"
    # Using mapi instead of map to limit tasks. Bash interaction not working as expected
    posteriors, posteriors_adjusted = unzip(ThreadsX.map((x, y) -> abc_loclinear(x,y, P=P, tol=tol, transformation=transformation, kernel=kernel), targets, summstat))

    @info "Filtering posterior distributions"
    # Remove if some file wasn't computed
    posteriors = filter(!iszero, posteriors)
    posteriors_adjusted = filter(!iszero, posteriors_adjusted)

    # Remove summstat files
    # if rm_summaries
    #     rm.(filter(x -> occursin("summstat", x) || occursin("alphas_", x), readdir(output_folder, join = true)))
    # end

    set_num_threads(nthreads_og)

    return posteriors, posteriors_adjusted
end


function abc_loclinear(target::Matrix{Float64},summaries::Matrix{Float64};P::Int64,tol::Float64,transformation::String="none",kernel::String= "epanechnikov")

    # Reading into matrix using Tables.matrix
    # target = vec(CSV.read(target_file,matrix,ntasks=1,header=false))
    # summaries = CSV.read(param_summaries_file,matrix,ntasks=1,header=false)
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
    target_scaled, sumstat_accepted, sumstat_scaled, param_accepted, dist_accepted, wts = MKtest.rejection(vec(target),sumstat,param,tol,kernel)
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
    lm_coefficients, lm_residuals = MKtest.regression(sumstat_intercept,param_accepted,wts)


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
    lm_coefficients, lm_residuals = MKtest.regression(sumstat_intercept,rsdl_log,wts)

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

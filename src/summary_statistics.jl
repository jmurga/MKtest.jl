################################
###    Summary statistics    ###
################################
"""
    poisson_fixation(observed_values,λds, λdn)

Divergence sampling from Poisson distribution. The expected neutral and selected fixations are subset through their relative expected rates ([`fix_neut`](@ref), [`fix_neg_b`](@ref), [`fix_pos_sim`](@ref)). Empirical values are used are used to simulate the locus *L* along a branch of time *T* from which the expected *Ds* and *Dn* raw count estimated given the mutation rate (``\\mu``). Random number generation is used to subset samples arbitrarily given the success rate ``\\lambda`` in the distribution.

```math
\\mathbb{E}[D_N] = X \\in Poisson\\left(\\lambda = D \\times \\left[\\frac{\\mathbb{E}[D_+] + \\mathbb{E}[D_-]}{\\mathbb{E}[D_+] + \\mathbb{E}[D_-] + \\mathbb{E}[D_0]}\\right]\\right)
```
```math
\\mathbb{E}[D_S] = X \\in Poisson\\left(\\lambda = D \\times \\left[\\frac{\\mathbb{E}[D_0]}{\\mathbb{E}[D_+] + \\mathbb{E}[D_-] + \\mathbb{E}[D_0]}\\right]\\right)
```
# Arguments
 - `observed_values::Array`: Array containing the total observed divergence.
 - ` λds::Float64`: expected neutral fixations rate.
 - ` λdn::Float64`: expected selected fixations rate.
# Returns
 - `Array{Int64,1}` containing the expected count of neutral and selected fixations.

"""
function poisson_fixation(
    observed_values::Vector,
    λds::SubArray,
    λdn::SubArray,
    λweak::SubArray,
    λstrong::SubArray,
)
    observed_d, ln, ls = observed_values

    ds                 = @. (λds / (λds + λdn)) * (observed_d)
    dn                 = @. (λdn / (λds + λdn)) * (observed_d)
    dweak              = @. (λweak / (λds + λdn)) * (observed_d)
    dstrong            = @. (λstrong / (λds + λdn)) * (observed_d)

    sampled_ds         = @. pois_rand(ds)
    sampled_dn         = @. pois_rand(dn)
    sampled_weak       = @. pois_rand(dweak)
    sampled_strong     = @. pois_rand(dstrong)

    α                  = @. [sampled_weak / sampled_dn sampled_strong / sampled_dn (
        sampled_weak + sampled_strong
    ) / sampled_dn]

    dₙ                 = @. (sampled_dn/ln)
    dₛ                 = @. (sampled_ds/ls)
    d₋                 = @. (sampled_dn-sampled_strong-sampled_weak)/ln
    D₋                 = @. sampled_dn-sampled_strong-sampled_weak
    D₊                 = @. sampled_strong+sampled_weak

    # ω                = @. (sampled_dn/ln)/(sampled_ds/ls)
    # ωₐ               = @. α *  ω
    # ωₙ               = @. (D₋/ln) / dₛ
    # ω_params         = hcat(ωₐ,ωₙ)

    ω                  = @. dₙ / dₛ
    ωₙ                 = @. (D₋/ln) / dₛ
    ωₐ_weak            = @. (sampled_weak/ln) / dₛ
    ωₐ_strong          = @. (sampled_strong/ln) / dₛ
    ωₐ                 = @. ((sampled_weak+sampled_strong)/ln) / dₛ
    ω_params           = hcat(ωₐ_weak,ωₐ_strong,ωₐ,ωₙ)

    out                = α, ω_params, sampled_dn, sampled_ds
    # out              = sampled_dn, sampled_ds, dₙ./dₛ
    return  out
end

"""
    poisson_polymorphism(observed_values,λps,λpn)

Polymorphism sampling from Poisson distributions. The total expected neutral and selected polimorphism are subset through the relative expected rates at the frequency spectrum ([`fix_neut`](@ref), [`sfs_neut`](@ref),). Empirical sfs are used to simulate the locus *L* along a branch of time *T* from which the expected *Ps* and *Pn* raw count are estimated given the mutation rate (``\\mu``). Random number generation is used to subset samples arbitrarily from the whole sfs given each frequency success rate ``\\lambda`` in the distribution.

The success rate managing the Poisson distribution by the observed count each frequency.  We considered both sampling variance and process variance is affecting the number of variable alleles we sample from SFS. This variance arises from the random mutation-fixation process along the branch. To incorporate this variance we do one sample per frequency-bin and use the total sampled variation and the SFS for the summary statistics.

```math
\\mathbb{E}[P_N]   = \\sum_{x=0}^{x=1} X \\in Poisson\\left(\\lambda = SFS_{(x)} \\times \\left[\\frac{\\mathbb{E}[P_{+(x)}] + \\mathbb{E}[P_{-(x)}]}{\\mathbb{E}[P_{+(x)}] + \\mathbb{E}[P_{-(x)}] + \\mathbb{E}[P_{0(x)}]}\\right]\\right)
```

```math
\\mathbb{E}[P_S] = \\sum_{x=0}^{x=1} X \\in Poisson\\left(\\lambda = SFS_{(x)} \\times \\left[\\frac{\\mathbb{E}[P_{0(x)}]}{\\mathbb{E}[P_{+(x)}] + \\mathbb{E}[P_{-(x)}] + \\mathbb{E}[P_{0(x)}]}\\right]\\right)
```

# Arguments
 - `observed_values::Array{Int64,1}`: Array containing the total observed divergence.
 - ` λps::Array{Float64,1} `: expected neutral site frequency spectrum rate.
 - ` λpn::Array{Float64,1} `: expected selected site frequency spectrum rate.
# Returns
 - `Array{Int64,2}` containing the expected total count of neutral and selected polymorphism.

"""
function poisson_polymorphism(
    observed_values::Vector,
    λps::Matrix{Float64},
    λpn::Matrix{Float64},
)

    # r = ifelse(sum(observed_values) < 1e4,fill(1000,length(observed_values)),fill(1,length(observed_values)))
    # Neutral λ;
    λ1 = @. (λps / (λps + λpn)) * (observed_values)
    # Selected λ;
    λ2 = @. (λpn / (λps + λpn)) * (observed_values)

    # Relative rates output NaN values due to 0 divisons.
    replace!(λ1, NaN => 1)
    replace!(λ2, NaN => 1)

    sampled_ps = @. pois_rand(λ1)
    sampled_pn = @. pois_rand(λ2)

    return (sampled_pn, sampled_ps)
end

"""
    sampled_alpha(observed_values,λds, λdn)

Ouput the expected values from the Poisson sampling process. Please check [`poisson_fixation`](@ref) and [`poisson_polymorphism`](@ref) to understand the samplingn process. α(x) is estimated through the expected values of Dn, Ds, Pn and Ps.

# Arguments
 - `param::parameters`: Array containing the total observed divergence.
 - `d::Array`: observed divergence.
 - `afs::Array`: observed polymorphism.
 - ` λdiv::Array{Float64,2}`: expected fixations rate.
 - ` λdiv::Array{Float64,2}`: expected site frequency spectrum rates.
# Returns
α_summaries,exp_dn,exp_ds,exp_pn,exp_ps,ssAlpha
 - `Array{Int64,2}` containing α(x) values.
 - `Array{Int64,1}` expected non-synonymous divergence.
 - `Array{Int64,1}` expected synonymous divergence.
 - `Array{Int64,1}` expected non-synonymous polymorphism.
 - `Array{Int64,1}` expected synonymous polymorphism.
 - `Array{Int64,1}` expected synonymous polymorphism.
 - `Array{Int64,1}` expected synonymous polymorphism.
 - `Array{Int64,2}` containing α(x) binned values.

    sampling_summaries(gammaL,gammaH,ppos_l,ppos_h,observedData,nopos)

"""
function sampling_summaries(
    models::SharedMatrix,
    fs::Vector,
    d::Vector,
    neut::SharedMatrix,
    sel::SharedMatrix,
    dsdn::SharedMatrix,
    output::String,
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
    α, ω, expected_dn, expected_ds = poisson_fixation(d, ds, dn, dweak, dstrong)
    # expected_dn, expected_ds, ω = poisson_fixation(d, ds, dn, dweak, dstrong)
    expected_pn, expected_ps = poisson_polymorphism(fs, permutedims(neut), permutedims(sel))
    # expected_pn, expected_ps = poisson_polymorphism(fs,d[2:3],permutedims(neut), permutedims(sel))

    ## Alpha from expected values. Used as summary statistics
    α_summaries = @. round(1 - ((expected_ds / expected_dn) * (expected_pn / expected_ps)'), digits = 5)
    if any(isnan.(ω))
        expected_values = hcat(round.(α, digits = 5), gn, sh, gH, gL, B, α_summaries)
    else
        # α = hcat(view(models,:,2),view(models,:,3).-view(models,:,2),view(models,:,3))
        # ωₐ = α .* ω
        # ωₙ = ω .- view(ωₐ,:,3)
        # ω = hcat(ωₐ,ωₙ)
        expected_values = hcat(round.(α, digits = 5),round.(ω, digits = 5), gn, sh, gL, gH, B, α_summaries)
    end

    expected_values = filter_expected(expected_values)

    write_files(expected_values, output)
end

"""
    summary_statistics(param,h5_file,sfs,divergence,output_folder,summstat_size)

Estimate summary statistics using observed data and analytical rates. *output_folder* will be used to output summary statistics

# Arguments
 - `param::parameters` : Mutable structure containing the models
 - `h5_file::String` : HDF5 containing solved models, fixation and polymorphic rates
 - `sfs::Vector`: SFS data parsed from parse_sfs().
 - `divergence::Vector`: divergence data from parse_sfs().
 - `output_folder::String` : path to save summary statistics and observed data.
 - `summstat_size::Int64` : number of summary statistics to sample.
"""
function summary_statistics(
    param::parameters,
    sfs::Vector,
    divergence::Vector;
    h5_file::String,
    summstat_size::Int64,
    output_folder::String,
    alpha::Union{Nothing,Vector{Float64}} = nothing,
    B_bins::Union{Nothing,Vector{Float64}} = nothing
)

    @assert 1 ∉ param.dac "Please exclude singletons from the inference"

    ## Opening files
    assertion_params(param)
    mkpath(output_folder)

    α, sfs_p, divergence_p = data_to_poisson(sfs, divergence, param.dac)

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

    #=# Compute the weights for each value in the vector
    ws = (1/(alpha[end]-alpha[1])) * (length(idx)/summstat_size) * ones(length(idx))

    # Sample N values from the vector with weights. Uniform prior alpha
    sampled_idx = sample(idx, Weights(ws), summstat_size;replace=true)=#

    @assert length(sampled_idx) >= summstat_size "Check the filters or increase the number of rates solutions"

    # Convert random models to solve to a SharedMatrix. Only reading shouldn't generate race-condition
    models = SharedMatrix(Array(view(tmp["models"], sampled_idx, :)))
    dsdn   = SharedMatrix(tmp["dsdn"][sampled_idx, :])

    # Filtering polymorphic rate by dac
    #=n = hcat(map(x -> view(tmp["neut"][x], :), param.dac)...)
    s = hcat(map(x -> view(tmp["sel"][x], :), param.dac)...)

    neut = SharedMatrix(n[sampled_idx, :])
    sel = SharedMatrix(s[sampled_idx, :])=#

    neut = SharedMatrix(zeros(summstat_size,length(param.dac)))
    sel = SharedMatrix(zeros(summstat_size,length(param.dac)))

    for (i,v) in enumerate(param.dac)
        neut[:,i] .= view(tmp["neut"][v],sampled_idx,:)
        sel[:,i] .= view(tmp["sel"][v],sampled_idx,:)
    end


    # Making summaries
    summ_output = output_folder .* "/summstat_" .* string.(1:size(sfs_p, 1)) .* ".txt"

    @info "Sampling and writting summary statistics at $output_folder"

    # 20 threads: case + control ~ 25GB
    expected_values = ThreadsX.map(
        (x, y, z) -> sampling_summaries(models, x, y, neut, sel, dsdn, z),
        sfs_p,
        divergence_p,
        summ_output
    );

    α_output = output_folder * "/alphas_" .* string.(1:size(sfs_p, 1)) .* ".txt"

    @. write_files(α, α_output)
end


function filter_expected(x::Matrix{Float64})
    replace!(x, -Inf => NaN)
    x = @view x[vec(.!any(isnan.(x), dims = 2)), :]

    return (Matrix(x))
end

function pol_correction!(
    sfs_all::Vector{Matrix{Float64}},
    sfs_in::Vector{Matrix{Float64}};
    column::Vector{Int} = [2],
)
    pn_all = map(x -> sum(x[:, column]), sfs_all)
    sfs_pn = map(x -> sum(x[:, column]), sfs_in)
    ratio_all_sub = map((x, y) -> x / y, pn_all, sfs_pn)

    map((x, y) -> y[:, column] .= x * y[:, column], ratio_all_sub, sfs_in)
end




# if any(data_filter)
#     sfs_p = sfs_p[.!data_filter]
#     divergence_p = divergence_p[.!data_filter]
#     α=α[.!data_filter]
#     @warn "Files distribution were reduce from " * string(length(data_filter)) * " values to " * string(sum(.!data_filter)) *", to exclude SFS or divergence distributions containing 0 values"

#     if isempty(sfs_p) | isempty(divergence_p)
#         throw(
#             ArgumentError(
#             "Your SFS contains 0 values at the selected DACs or the divergence is 0. Please consider to bin the SFS and re-estimate the rates using the selected bin as sample the new sample size.",
#         ),
#     )
#     end
# else all(data_filter)

#     throw(
#         ArgumentError(
#             "Your SFS contains 0 values at the selected DACs or the divergence is 0. Please consider to bin the SFS and re-estimate the rates using the selected bin as sample the new sample size.",
#         ),
#     )
# end

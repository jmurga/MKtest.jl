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
    observed_values::Vector{Int64},
    λds::SubArray,
    λdn::SubArray,
    λweak::SubArray,
    λstrong::SubArray,
)
    ds = @. λds / (λds + λdn) * observed_values
    dn = @. λdn / (λds + λdn) * observed_values
    dweak = @. λweak / (λds + λdn) * observed_values
    dstrong = @. λstrong / (λds + λdn) * observed_values

    sampled_ds = pois_rand.(ds)
    sampled_dn = pois_rand.(dn)
    sampled_weak = pois_rand.(dweak)
    sampled_strong = pois_rand.(dstrong)

    alphas = @. [sampled_weak / sampled_dn sampled_strong / sampled_dn (
        sampled_weak + sampled_strong
    ) / sampled_dn]

    out = alphas, sampled_dn, sampled_ds
    return out
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
    observed_values::Vector{Float64},
    λps::Matrix{Float64},
    λpn::Matrix{Float64},
)

    # Neutral λ;
    λ1 = @. λps / (λps + λpn) * observed_values
    # Selected λ;
    λ2 = @. λpn / (λps + λpn) * observed_values

    # Replace negative and 0 rates values. Strange behaviour depending on θ and Γ values
    # Relative rates output NaN values due to 0 divisons.

    replace!(λ1, NaN => 1)
    replace!(λ2, NaN => 1)

    sampled_ps = pois_rand.(λ1)
    sampled_pn = pois_rand.(λ2)

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
    fs::Vector{Float64},
    d::Vector{Int64},
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

    ## Outputs
    alphas, exp_dn, exp_ds = poisson_fixation(d, ds, dn, dweak, dstrong)
    exp_pn, exp_ps = poisson_polymorphism(fs, permutedims(neut), permutedims(sel))

    ## Alpha from expected values. Used as summary statistics
    α_summaries = @. round(1 - ((exp_ds / exp_dn) * (exp_pn / exp_ps)'), digits = 5)
    expected_values = hcat(round.(alphas, digits = 5), gn, sh, α_summaries)
    expected_values = filter_expected(expected_values)

    write_files(expected_values, output)
end

"""
	summary_statistics(param,h5_file,sfs,divergence,analysis_folder,summstat_size)

Estimate summary statistics using observed data and analytical rates. *analysis_folder* will be used to output summary statistics

# Arguments
 - `param::parameters` : Mutable structure containing the models
 - `h5_file::String` : HDF5 containing solved models, fixation and polymorphic rates
 - `sfs::Vector`: SFS data parsed from parse_sfs().
 - `divergence::Vector`: divergence data from parse_sfs().
 - `analysis_folder::String` : path to save summary statistics and observed data.
 - `summstat_size::Int64` : number of summary statistics to sample.
"""
function summary_statistics(
    param::parameters;
    h5_file::String,
    sfs::Vector,
    divergence::Vector,
    analysis_folder::String,
    summstat_size::Int64,
)

    ## Opening files
    # if length(sfs) > 1
    # 	throw(ArgumentError("You have more than one SFS and divergence file. Please be sure you have on set of files to bootstrap manually your data."))
    # end

    assertion_params(param)

    α, sfs_p, divergence_p = data_to_poisson(sfs, divergence, param.dac)

    if any(0 .∈ sfs_p) | any(0 .∈ divergence_p)
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
    idx = sample(1:size(tmp["models"], 1), summstat_size, replace = false)
    # Convert random models to solve to a SharedMatrix. Only reading shouldn't generate race-condition
    models = SharedMatrix(Array(view(tmp["models"], idx, :)))
    dsdn = SharedMatrix(tmp["dsdn"][idx, :])

    # Filtering polymorphic rate by dac
    n = hcat(map(x -> view(tmp["neut"][x], :), param.dac)...)
    s = hcat(map(x -> view(tmp["sel"][x], :), param.dac)...)

    neut = SharedMatrix(n[idx, :])
    sel = SharedMatrix(s[idx, :])

    # Making summaries
    summ_output = analysis_folder .* "/summstat_" .* string.(1:size(sfs_p, 1)) .* ".txt"

    @info "Sampling and writting summary statistics at $analysis_folder"

    # 20 threads: case + control ~ 25GB
    expected_values = ThreadsX.map(
        (x, y, z) -> sampling_summaries(models, x, y, neut, sel, dsdn, z),
        sfs_p,
        divergence_p,
        summ_output,
    )

    α_output = analysis_folder * "/alphas_" .* string.(1:size(sfs_p, 1)) .* ".txt"

    @. write_files(α, α_output)
end

function filter_expected(x::Matrix{Float64})
    replace!(x, -Inf => NaN)
    x = @view x[vec(.!any(isnan.(x), dims = 2)), :]
    x = @view x[(x[:, 3].<1), :]
    x = @view x[(x[:, 1].>0), :]

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

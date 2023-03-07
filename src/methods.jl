"""
	α(x)(sfs,divergence)

Function to estimate alpha_x.

# Arguments
 - `sfs::Matrix`: SFS data Matrix parsed from parse_sfs().
 - `divergence::Matrix`: divergence data Matrix from parse_sfs().
 - `cumulative::Bool`: perform the estimation using the cumulative SFS.
# Returns
 - `Vector{Float64}: α(x) estimations.
"""
function α_x(sfs::Matrix, divergence::Matrix; cumulative::Bool = true)
    if cumulative
        sfs = cumulative_sfs(sfs)
    end

    alpha_x = @. 1 - (sfs[:, 2] / sfs[:, 3] * divergence[2] / divergence[1])

    return (alpha_x)
end

"""
    α(x)(sfs,divergence)

Function to estimate alpha_x.

# Arguments
 - `sfs::Vector`: Multiple SFS data Matrix parsed from parse_sfs().
 - `divergence::Vector`: Multiple divergence data Matrix from parse_sfs().
 - `cumulative::Bool`: perform the estimation using the cumulative SFS.
# Returns
 - `Vector{Float64}: α(x) estimations.
"""
function α_x(sfs::Vector, divergence::Vector; cumulative::Bool = true)

    output = Vector{Vector}(undef, length(sfs))
    @unpack nn, dac = param
    Threads.@threads for i in eachindex(sfs)
        if cumulative
            sfs = cumulative_sfs(sfs)
        end

        alpha_x = @. 1 - (sfs[:, 2] / sfs[:, 3] * divergence[2] / divergence[1])
        output[i] = alpha_x
    end
    return (alpha_x)
end

"""
	aMK(param,α;dac_rm,na_rm)

Function to estimate the asymptotic value of α(x).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model
 - `α::Matrix`: α(x) values.
 - `dac_rm::Bool`: Remove any value not selected at param.dac.
 - `na_rm::Bool`: Remove any NA value at α(x)

# Returns
 - `Tuple{Float64, Vector{Float64}, Vector{Float64}}`: Tuple containing aMK estimation, CI estimation, model output.
"""
function aMK(
    param::parameters,
    α::Vector{Float64};
    dac_rm::Bool = false,
    na_rm::Bool = false,
)
    @unpack nn, dac = param

    f = collect(1:length(α)) / nn

    if dac_rm
        f = dac / nn
        α = α[dac]
    end

    if na_rm
        flt = .!isnan.(α) .&& .!isinf.(α) .&& α .!= 1
        α = α[flt]
        f = f[flt]
    end

    # Model
    model(x, p) = @. p[1] + p[2] * exp(-x * p[3])

    # Fit values

    fitted1 = curve_fit(
        model,
        f,
        α,
        [-1.0, -1.0, 1.0];
        lower = [-1.0, -1.0, 1.0],
        upper = [1.0, 1.0, 10.0],
    )
    fitted2 = curve_fit(model, f, α, fitted1.param)
    asymp = model(1, fitted2.param)

    ci = try
        [confidence_interval(fitted2)[1][1], confidence_interval(fitted2)[1][2]]
    catch
        [0.0, 0.0]
    end

    return (asymp, ci, fitted2.param)
end


"""
    aMK(param,α;dac_rm,na_rm)

Function to estimate the asymptotic value of α(x).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model
 - `α::Vector`: Multiple α(x) vectors.
 - `dac_rm::Bool`: Remove any value not selected at param.dac.
 - `na_rm::Bool`: Remove any NA value at α(x)

# Returns
 - `Vector{Tuple}`: Tuple containing aMK estimation, CI estimation, model output.
"""
function aMK(param::parameters, α::Vector; dac_rm::Bool = false, na_rm::Bool = false)

    output = Vector{Tuple}(undef, length(α))

    @unpack nn, dac = param

    f = Vector{Vector{Float64}}(undef, length(α))
    for i ∈ eachindex(α)
        tmp = collect(1:length(α[i])) / nn

        if dac_rm
            tmp = dac / nn
            α[i] = α[i][dac]
        end

        if na_rm
            flt = .!isnan.(α[i]) .&& .!isinf.(α[i]) .&& α[i] .!= 1
            α[i] = α[i][flt]
            tmp = tmp[flt]
        end
        f[i] = tmp
    end

    # Model
    model(x, p) = @. p[1] + p[2] * exp(-x * p[3])

    # Fit values
    fitted1 = ThreadsX.map((x,y) -> curve_fit(model,y,x,[-1.0, -1.0, 1.0];lower = [-1.0, -1.0, 1.0],upper = [1.0, 1.0, 10.0]),α,f)
    fitted2 = ThreadsX.map((x,y,z) -> curve_fit(model, y,x, z.param),α,f,fitted1)
    asymp = ThreadsX.map(x -> model(1, x.param),fitted2)

    ci = ThreadsX.map(x -> confidence_interval(x)[1],fitted2)
    model_params = map(x-> x.param,fitted2)

    return (asymp,ci,model_params)
end

"""
	imputedMK(param,sfs,divergence;m,cutoff)

Function to estimate the imputedMK (https://doi.org/10.1093/g3journal/jkac206).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model.
 - `sfs::Matrix{Float64}`: SFS data parsed from parse_sfs().
 - `divergence::Matrix{Int64}`: divergence data from parse_sfs().
 - `cutoff{Float64}`: frequency cutoff to perform imputedMK.
 - `m::Matrix{Int64}`: total number of sites parsed from parse_sfs()
# Returns
 - `Dict: Dictionary containing imputedMK estimation.
"""

function imputedMK(
    param::parameters,
    sfs::Vector,
    divergence::Vector,
    cutoff::Float64 = 0.15,
    m::T = nothing,
) where {T<:Union{Nothing,Array}}
    output = Vector{OrderedDict}(undef, length(sfs))
    tmp = OrderedDict{String,Float64}()

    @unpack isolines = param
    Threads.@threads for i ∈ eachindex(sfs)

        pn = sum(sfs[i][:, 2])
        ps = sum(sfs[i][:, 3])
        dn = divergence[i][1]
        ds = divergence[i][2]

        deleterious = 0
        ### Estimating slightly deleterious with pn/ps ratio
        sfs_tmp = deepcopy(sfs[i])
        s_size = if isolines
            param.n
        else
            param.nn
        end

        sfs_tmp[:, 1] = sfs_tmp[:, 1] ./ s_size

        flt_low = (sfs_tmp[:, 1] .<= cutoff)
        pn_low = sum(sfs_tmp[flt_low, 2])
        ps_low = sum(sfs_tmp[flt_low, 3])

        flt_inter = (sfs_tmp[:, 1] .>= cutoff) .& (sfs_tmp[:, 1] .<= 1)
        pn_inter = sum(sfs_tmp[flt_inter, 2])
        ps_inter = sum(sfs_tmp[flt_inter, 3])

        ratio_ps = ps_low / ps_inter
        deleterious = pn_low - (pn_inter * ratio_ps)

        if (deleterious > pn) || (deleterious < 0)
            deleterious = 0
            pn_neutral = round(pn - deleterious, digits = 3)
        else
            pn_neutral = round(pn - deleterious, digits = 3)
        end

        tmp["alpha"] = round(1 - ((pn_neutral / ps) * (ds / dn)), digits = 3)

        #  method = :minlike same results R, python two.sides
        tmp["pvalue"] =
            pvalue(FisherExactTest(Int(ps), Int(ceil(pn_neutral)), Int(ds), Int(dn)))

        if (!isnothing(m))
            mn = m[1]
            ms = m[2]
            # ## Estimation of b: weakly deleterious
            tmp["b"] = (deleterious / ps) * (ms / mn)

            ## Estimation of f: neutral sites
            tmp["f"] = (ms * pn_neutral) / (mn * ps)

            ## Estimation of d, strongly deleterious sites
            tmp["d"] = 1 - (tmp["f"] + tmp["b"])

            ka = dn / mn
            ks = ds / ms
            tmp["omega"] = ka / ks

            # Omega A and Omega D
            tmp["omegaA"] = tmp["omega"] * tmp["alpha"]
            tmp["omegaD"] = tmp["omega"] - tmp["omegaA"]
        end
        output[i] = tmp
    end
    return output
end
"""
	fwwMK(param,sfs,divergence;m,cutoff)

Function to estimate the fwwMK (https://doi.org/10.1038/4151024a).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model.
 - `sfs::Matrix{Float64}`: SFS data parsed from parse_sfs().
 - `divergence::Matrix{Int64}`: divergence data from parse_sfs().
 - `cutoff{Float64}`: frequency cutoff to perform imputedMK.
 - `m::Matrix{Int64}`: total number of sites parsed from parse_sfs()
# Returns
 - `Dict: Dictionary containing imputedMK estimation.
"""
function fwwMK(
    param::parameters,
    sfs::Vector,
    divergence::Vector;
    cutoff::Float64 = 0.15,
    m::T = nothing,
) where {T<:Union{Nothing,Array}}
    output = Vector{OrderedDict}(undef, length(sfs))
    tmp = OrderedDict{String,Float64}()

    @unpack isolines = param;
    Threads.@threads for i ∈ eachindex(sfs)
        ps = sum(sfs[i][:, 2])
        pn = sum(sfs[i][:, 3])
        dn = divergence[i][1]
        ds = divergence[i][2]

        deleterious = 0
        ### Estimating slightly deleterious with pn/ps ratio
        s_size = if isolines
            param.n
        else
            param.nn
        end
        sfs_tmp = deepcopy(sfs[i])
        sfs_tmp[:, 1] = sfs_tmp[:, 1] ./ s_size
        flt_inter = (sfs_tmp[:, 1] .>= cutoff) .& (sfs_tmp[:, 1] .<= 1)
        pn_inter = sum(sfs_tmp[flt_inter, 2])
        ps_inter = sum(sfs_tmp[flt_inter, 3])


        tmp["alpha"] = round(1 - ((pn_inter / ps_inter) * (ds / dn)), digits = 3)
        #  method = :minlike same results R, python two.sides
        tmp["pvalue"] =
            pvalue(FisherExactTest(Int(ps_inter), Int(ceil(pn_inter)), Int(ds), Int(dn)))

        if (!isnothing(m))
            mn = m[1]
            ms = m[2]
            ka = dn / mn
            ks = ds / ms
            tmp["omega"] = ka / ks

            # Omega A and Omega D
            tmp["omegaA"] = tmp["omega"] * tmp["alpha"]
            tmp["omegaD"] = tmp["omega"] - tmp["omegaA"]
        end
        output[i] = tmp
    end
    return output
end

"""
	standardMK(sfs,divergence;m)

Function to estimate the original α value.

# Arguments
 - `sfs::Matrix{Float64}`: SFS array
 - `divergence::Array`: divegence count array
 - `m::Union{Nothing,Array}`: non-synonymous and synonymous sites# Returns

# Output
 - `Dict: Dictionary containing results
"""
function standardMK(
    sfs::Vector,
    divergence::Vector;
    m::T = nothing,
) where {T<:Union{Nothing,Array}}
    output = Vector{OrderedDict}(undef, length(sfs))
    tmp = OrderedDict{String,Float64}()

    Threads.@threads for i ∈ eachindex(sfs)
        pn = sum(sfs[i][:, 2])
        ps = sum(sfs[i][:, 3])
        dn = divergence[i][1]
        ds = divergence[i][2]

        tmp["alpha"] = round(1 - ((pn / ps) * (ds / dn)), digits = 3)
        #  method = :mnnlike same results R, python two.sides
        tmp["pvalue"] = pvalue(FisherExactTest(Int(ps), Int(ceil(pn)), Int(ds), Int(dn)))

        if (!isnothing(m))
            mn = m[1]
            ms = m[2]

            ka = dn / mn
            ks = ds / ms
            tmp["omega"] = ka / ks

            # Omega A and Omega D
            tmp["omegaA"] = tmp["omega"] * tmp["alpha"]
            tmp["omegaD"] = tmp["omega"] - tmp["omegaA"]
        end
        output[i] = tmp
    end
    return output
end

"""
    grapes(sfs,divergence,m,model,folder,bins)

Run grapes using SFS and divergence data. The function will install grapes using Conda from genomedk channel. 

# Arguments
 - `sfs::Vector{Matrix{Float64}}`: SFS data parsed from parse_sfs()
 - `divergence::Vector{Matrix{Int64}}`: divergence data parsed from parse_sfs()
 - `m::Vector{Matrix{Int64}}`: total number of sites parsed from parse_sfs()
 - `model::String`: grapes model.
 - `folder::String`: output folder.
 - `bins::Int64`:  bin size to reduce the SFS. Grapes became unstable when inputing large SFS.

# Output
 - `DataFrame`: grapes model estimation.
"""

function grapes(
    sfs::Vector,
    divergence::Vector,
    m::Vector,
    model::String,
    folder::String,
    bins::Int64;
    nearly_neutral::Int64 = 5,
    FWW_threshold::Float64 = 0.15,
    nb_rand_start::Int64 = 0,
    anc_to_rec_Ne_ratio::Float64 = 1.0,
    no_div_data::Bool = false,
    no_div_param::Bool = false,
    no_syn_orient_error::Bool = false,
    fold::Bool = false,
    fixed_param::String = "",
)

    grapes_bin = CondaPkg.which("grapes")

    if (isnothing(grapes_bin))
        CondaPkg.add("grapes-static", channel = "genomedk")
        grapes_bin = CondaPkg.which("grapes")
    end

    @assert model ∈ [
        "GammaZero",
        "GammaExpo",
        "GammaGamma",
        "DisplGamma",
        "ScaledBeta",
        "FGMBesselK",
        "all",
    ] "Please select a valid model: GammaZero GammaExpo GammaGamma DisplGamma ScaledBeta FGMBesselK all"

    sfs = reduce_sfs.(sfs, bins)

    pn = map(x -> permutedims(x[:, 2]), sfs)
    ps = map(x -> permutedims(x[:, 3]), sfs)

    dn = map(x -> x[1], divergence)
    ds = map(x -> x[2], divergence)

    mn = map(x -> x[1], m)
    ms = map(x -> x[2], m)

    idx = string.(collect(1:length(sfs)))

    # Temporal function to broadcast pol and div data
    function f(pn, ps, dn, ds, mn, ms, w)
        DataFrame(hcat("dofe_" * string(w), bins, mn, pn, ms, ps, mn, dn, ms, ds...), :auto)
    end

    @info "Converting SFS to dofe file"
    dofe = @. f(pn, ps, dn, ds, mn, ms, idx)
    h = fill(DataFrame(["" ""; "#unfolded" ""], :auto), length(sfs))
    output_dofe = @. folder * "/dofe_" * idx * ".txt"
    output_grapes = @. folder * "/dofe_" * idx * "." * model

    @. write_files(h, output_dofe, false)
    @. write_files(dofe, output_dofe, true)

    @info "Running Grapes"
    if model == "GammaZero" || no_div_param
        no_div_param = false
    end

    bool_options = ["-no_div_data", "-no_div_param", "-no_syn_orient_error", "-fold"]
    bool_options = join(
        bool_options[any(
            vcat(no_div_data, no_div_param, no_syn_orient_error, fold),
            dims = 2,
        )],
        " ",
    )

    if !isempty(fixed_param)
        fixed_param = "-fixed_param " * fixed_param
    end


    r(
        d,
        o,
        md = model,
        gr = grapes_bin,
        nt = nearly_neutral,
        fww = FWW_threshold,
        nb = nb_rand_start,
        ne = anc_to_rec_Ne_ratio,
        bopt = bool_options,
        fp = fixed_param,
    ) = run(
        `$gr -in $d -out $o -model $md -nearly_neutral $nt -FWW_threshold $fww -nb_rand_start $nb -anc_to_rec_Ne_ratio $ne $bopt $fp`,
    )

    @suppress_out begin
        ThreadsX.mapi(r, output_dofe, output_grapes, ntasks = Threads.nthreads())
    end

    if model == "all"
        df = @suppress begin
            CSV.read.(output_grapes, DataFrame)
        end
    else
        df = @suppress begin
            CSV.read.(output_grapes, DataFrame, footerskip = 1, skipto = 3)
        end
    end
    return (vcat(df...))
end

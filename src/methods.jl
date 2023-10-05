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

    α_i = deepcopy(α)

    if dac_rm
        f = dac / nn
        α_i = α_i[dac]
    end

    if na_rm
        flt = .!isnan.(α_i) .&& .!isinf.(α_i) .&& α_i .!= 1
        α_i = α_i[flt]
        f = f[flt]
    end

    # Model
    model(x, p) = @. p[1] + p[2] * exp(-x * p[3])

    # Fit values

    fitted1 = curve_fit(
        model,
        f,
        α_i,
        [-1.0, -1.0, 1.0];
        lower = [-1.0, -1.0, 1.0],
        upper = [1.0, 1.0, 10.0],
    )
    fitted2 = curve_fit(model, f, α_i, fitted1.param)
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

    α_i = deepcopy(α)

    f = Vector{Vector{Float64}}(undef, length(α))
    for i ∈ eachindex(α_i)
        tmp = collect(1:length(α_i[i])) / nn

        if dac_rm
            tmp = dac / nn
            α_i[i] = α_i[i][dac]
        end

        if na_rm
            flt = .!isnan.(α_i[i]) .&& .!isinf.(α_i[i]) .&& α_i[i] .!= 1
            α_i[i] = α_i[i][flt]
            tmp = tmp[flt]
        end
        f[i] = tmp
    end

    # Model
    model(x, p) = @. p[1] + p[2] * exp(-x * p[3])

    # Fit values
    amk_v = Vector{Float64}(undef, length(α))
    model_params = Vector{Vector{Float64}}(undef, length(α))
    ci_v = Vector{Vector{Float64}}(undef, length(α))

    Threads.@threads for i ∈ eachindex(α_i)
        fitted1 = curve_fit(model,f[i],α_i[i],[-1.0, -1.0, 1.0];lower = [-1.0, -1.0, 1.0],upper = [1.0, 1.0, 10.0])
        fitted2 = curve_fit(model, α_i[i],f[i], fitted1.param)
        asymp = model(1, fitted2.param)
        amk_v[i] = asymp

        model_params[i] = fitted2.param

        ci_v[i] = try
            [confidence_interval(fitted2)[1][1], confidence_interval(fitted2)[1][2]]
        catch
            [0.0, 0.0]
        end

    end

    return (amk_v,ci_v,model_params)
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
    divergence::Vector;
    cutoff::Float64 = 0.15,
    m::T = nothing,
) where {T<:Union{Nothing,Array}}
    output = Vector{OrderedDict}(undef, length(sfs))

    @unpack isolines = param
    Threads.@threads for i ∈ eachindex(sfs)
        tmp = OrderedDict{String,Float64}()

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

        flt_inter = (sfs_tmp[:, 1] .> cutoff) .& (sfs_tmp[:, 1] .<= 1)
        pn_inter = sum(sfs_tmp[flt_inter, 2])
        ps_inter = sum(sfs_tmp[flt_inter, 3])

        if (pn_inter == 0 ||  ps_inter == 0)
            @warn "There is no polymorphic data above the cutoff on dataset $i"
            tmp["pvalue"] = NaN
            tmp["alpha"] = NaN
            output[i] = tmp
        else
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
                ln = m[1]
                ls = m[2]
                # ## Estimation of b: weakly deleterious
                tmp["b"] = (deleterious / ps) * (ls / ln)

                ## Estimation of f: neutral sites
                tmp["f"] = (ls * pn_neutral) / (ln * ps)

                ## Estimation of d, strongly deleterious sites
                tmp["d"] = 1 - (tmp["f"] + tmp["b"])

                ka = dn / ln
                ks = ds / ls
                tmp["omega"] = ka / ks

                # Omega A and Omega D
                tmp["omegaA"] = tmp["omega"] * tmp["alpha"]
                tmp["omegaD"] = tmp["omega"] - tmp["omegaA"]
            end
            output[i] = tmp
        end
    end
    return output
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
    sfs::Matrix,
    divergence::Matrix;
    cutoff::Float64 = 0.15,
    m::T = nothing,
) where {T<:Union{Nothing,Array}}

    @unpack isolines = param

    output = OrderedDict{String,Float64}()

    pn = sum(sfs[:, 2])
    ps = sum(sfs[:, 3])
    dn = divergence[1]
    ds = divergence[2]

    deleterious = 0
    ### Estimating slightly deleterious with pn/ps ratio

    s_size = if isolines
        param.n
    else
        param.nn
    end

    sfs_tmp= deepcopy(sfs)
    sfs_tmp[:, 1] = sfs_tmp[:, 1] ./ s_size

    flt_low = (sfs_tmp[:, 1] .<= cutoff)
    pn_low = sum(sfs_tmp[flt_low, 2])
    ps_low = sum(sfs_tmp[flt_low, 3])

    flt_inter = (sfs_tmp[:, 1] .>= cutoff) .& (sfs_tmp[:, 1] .<= 1)
    pn_inter = sum(sfs_tmp[flt_inter, 2])
    ps_inter = sum(sfs_tmp[flt_inter, 3])

    @assert pn_inter != 0 ||  ps_inter != 0 "There is no polymorphic data above the cutoff"

    ratio_ps = ps_low / ps_inter
    deleterious = pn_low - (pn_inter * ratio_ps)

    if (deleterious > pn) || (deleterious < 0)
        deleterious = 0
        pn_neutral = round(pn - deleterious, digits = 3)
    else
        pn_neutral = round(pn - deleterious, digits = 3)
    end

    output["alpha"] = round(1 - ((pn_neutral / ps) * (ds / dn)), digits = 3)

    #  method = :minlike same results R, python two.sides
    output["pvalue"] =
        pvalue(FisherExactTest(Int(ps), Int(ceil(pn_neutral)), Int(ds), Int(dn)))

    if (!isnothing(m))
        ln = m[1]
        ls = m[2]
        # ## Estimation of b: weakly deleterious
        output["b"] = (deleterious / ps) * (ls / ln)

        ## Estimation of f: neutral sites
        output["f"] = (ls * pn_neutral) / (ln * ps)

        ## Estimation of d, strongly deleterious sites
        # output["d"] = 1 - (output["f"] + output["b"])

        ka = dn / ln
        ks = ds / ls
        output["omega"] = ka / ks

        # Omega A and Omega D
        output["omegaA"] = output["omega"] * output["alpha"]
        output["omegaD"] = output["omega"] - output["omegaA"]
    end

    return output
end

"""
	fwwMK(param,sfs,divergence;m,cutoff)

Function to estimate the fwwMK (https://doi.org/10.1038/4151024a).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model.
 - `sfs::Vector`: SFS data parsed from parse_sfs().
 - `divergence::Vector`: divergence data from parse_sfs().
 - `cutoff{Float64}`: frequency cutoff to perform fwwMK.
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

    @unpack isolines = param;
    Threads.@threads for i ∈ eachindex(sfs)
        tmp = OrderedDict{String,Float64}()
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


        if (pn_inter == 0 ||  ps_inter == 0)
            @warn "There is no polymorphic data above the cutoff on dataset $i"
            tmp["pvalue"] = NaN
            tmp["alpha"] = NaN
            output[i] = tmp
        else
            tmp["alpha"] = round(1 - ((pn_inter / ps_inter) * (ds / dn)), digits = 3)
            #  method = :minlike same results R, python two.sides
            tmp["pvalue"] =
                pvalue(FisherExactTest(Int(ps_inter), Int(ceil(pn_inter)), Int(ds), Int(dn)))

            if (!isnothing(m))
                ln = m[1]
                ls = m[2]
                ka = dn / ln
                ks = ds / ls
                tmp["omega"] = ka / ks

                # Omega A and Omega D
                tmp["omegaA"] = tmp["omega"] * tmp["alpha"]
                tmp["omegaD"] = tmp["omega"] - tmp["omegaA"]
            end
            output[i] = tmp
        end
    end
    return output
end

"""
    fwwMK(param,sfs,divergence;m,cutoff)

Function to estimate the fwwMK (https://doi.org/10.1038/4151024a).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model.
 - `sfs::Matrix`: SFS data parsed from parse_sfs().
 - `divergence::Matrix`: divergence data from parse_sfs().
 - `cutoff{Float64}`: frequency cutoff to perform fwwMK.
 - `m::Matrix{Int64}`: total number of sites parsed from parse_sfs()
# Returns
 - `Dict: Dictionary containing imputedMK estimation.
"""
function fwwMK(
    param::parameters,
    sfs::Matrix,
    divergence::Matrix;
    cutoff::Float64 = 0.15,
    m::T = nothing,
) where {T<:Union{Nothing,Array}}
    output = OrderedDict{String,Float64}()

    @unpack isolines = param;

    ps = sum(sfs[:, 2])
    pn = sum(sfs[:, 3])
    dn = divergence[1]
    ds = divergence[2]

    deleterious = 0
    ### Estimating slightly deleterious with pn/ps ratio
    s_size = if isolines
        param.n
    else
        param.nn
    end

    sfs[:, 1] = sfs[:, 1] ./ s_size
    flt_inter = (sfs[:, 1] .>= cutoff) .& (sfs[:, 1] .<= 1)
    pn_inter = sum(sfs[flt_inter, 2])
    ps_inter = sum(sfs[flt_inter, 3])

    @assert pn_inter != 0 ||  ps_inter != 0 "There is no polymorphic data above the cutoff"


    output["alpha"] = round(1 - ((pn_inter / ps_inter) * (ds / dn)), digits = 3)
    #  method = :minlike same results R, python two.sides
    output["pvalue"] =
        pvalue(FisherExactTest(Int(ps_inter), Int(ceil(pn_inter)), Int(ds), Int(dn)))

    if (!isnothing(m))
        ln = m[1]
        ls = m[2]
        ka = dn / ln
        ks = ds / ls
        output["omega"] = ka / ks

        # Omega A and Omega D
        output["omegaA"] = output["omega"] * output["alpha"]
        output["omegaD"] = output["omega"] - output["omegaA"]
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

    Threads.@threads for i ∈ eachindex(sfs)
        tmp = OrderedDict{String,Float64}()
        pn = sum(sfs[i][:, 2])
        ps = sum(sfs[i][:, 3])
        dn = divergence[i][1]
        ds = divergence[i][2]

        tmp["alpha"] = round(1 - ((pn / ps) * (ds / dn)), digits = 3)
        #  method = :lnnlike same results R, python two.sides
        tmp["pvalue"] = pvalue(FisherExactTest(Int(ps), Int(ceil(pn)), Int(ds), Int(dn)))

        if (!isnothing(m))
            ln = m[1]
            ls = m[2]

            ka = dn / ln
            ks = ds / ls
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
    model::String,
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

    if any(0 .∈ divergence)
        throw(
            ArgumentError(
                "Your SFS contains 0 values at the selected DACs or the divergence is 0. Please consider to bin the SFS and re-estimate the rates using the selected bin as sample the new sample size.",
            ),
        )
    end
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

    sfs_r = reduce_sfs.(sfs, bins)

    pn = map(x -> permutedims(x[:, 2]), sfs_r)
    ps = map(x -> permutedims(x[:, 3]), sfs_r)

    dn = map(x -> [x[1]], divergence)
    ds = map(x -> [x[2]], divergence)
    ln = map(x -> [x[3]], divergence)
    ls = map(x -> [x[4]], divergence)

    idx = string.(collect(1:length(sfs_r)))

    @info "Converting SFS to dofe file"
    dofe = @. sfs_to_dofe(pn, ps, dn, ds, ln, ls, idx,bins)
    h = fill(DataFrame(["" ""; "#unfolded" ""], :auto), length(sfs_r))

    output_dofe   = @. tempname() * idx
    output_grapes = @. output_dofe * "." * model

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
        ThreadsX.mapi((x,y) -> r(x,y), output_dofe, output_grapes; ntasks = Threads.nthreads())
    end

    df = map((x,y) -> read_clean_grapes(x,y,model),output_dofe,output_grapes)
    df = df[.!isnothing.(df)]

    return (vcat(df...))
end

function sfs_to_dofe(pn, ps, dn, ds, ln, ls, w, bins)
    return(DataFrame(hcat("dofe_" * string(w), bins, ln, pn, ls, ps, ln, dn, ls, ds...), :auto))
end

function read_clean_grapes(dofe::String,grapes_file::String,model::String)
    # Try to open file if grapes executed properly
    output = @suppress_err begin
        try
            if model == "all"
                df = CSV.read(grapes_file, DataFrame)
            else
                df = CSV.read(grapes_file, DataFrame, footerskip = 1, skipto = 3)
            end

            rm(dofe)
            rm(grapes_file)
            rm(grapes_file * ".profile")
            rm(grapes_file * ".messages")

            df

        catch
            @warn "$dofe did not complete the Grapes estimation. Please check your data."
            nothing
        end
    end

    return(output)
end

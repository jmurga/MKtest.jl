"""
	α(x)(sfs,divergence)

Function to estimate alpha_x.

# Arguments
 - `sfs::Matrix{Float64}`: SFS data parsed from parse_sfs().
 - `divergence::Matrix{Int64}`: divergence data from parse_sfs().
 - `cumulative::Bool`: perform the estimation using the cumulative SFS.
# Returns
 - `Vector{Float64}: α(x) estimations.
"""
function α_x(sfs::Matrix{Float64}, divergence::Matrix{Int64}; cumulative::Bool = true)
    if cumulative
        sfs = cumulative_sfs(sfs)
    end

    alpha_x = @. 1 - (sfs[:, 2] / sfs[:, 3] * divergence[2] / divergence[1])

    return (alpha_x)
end

"""
	aMK(param,α;dac_rm,na_rm)

Function to estimate the asymptotic value of α(x).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model
 - `α::Matrix{Float64}`: α(x) values.
 - `dac_rm::Bool`: Remove any value not selected at param.dac.
 - `na_rm::Bool`: Remove any NA value at α(x)

# Returns
 - `Tuple{Float64, Vector{Float64}, Vector{Float64}}`: Tuple containing aMK estimation, CI estimation, model output.
"""
function aMK(param::parameters,
             α::Vector{Float64};
             dac_rm::Bool = false,
             na_rm::Bool = false)
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

    fitted1 = curve_fit(model,
                        f,
                        α,
                        [-1.0, -1.0, 1.0];
                        lower = [-1.0, -1.0, 1.0],
                        upper = [1.0, 1.0, 10.0])
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
function imputedMK(param::parameters, sfs::Matrix{Float64}, divergence::Matrix{Int64};
                   cutoff::Float64 = 0.15,
                   m::T = nothing) where {T <: Union{Nothing, Array}}
    output = OrderedDict{String, Float64}()

    pn = sum(sfs[:, 2])
    ps = sum(sfs[:, 3])
    dn = divergence[1]
    ds = divergence[2]

    deleterious = 0
    ### Estimating slightly deleterious with pn/ps ratio
    sfs_tmp = deepcopy(sfs)
    sfs_tmp[:, 1] = sfs_tmp[:, 1] ./ param.nn
    println(sfs_tmp[1, 1])
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

    output["alpha"] = round(1 - ((pn_neutral / ps) * (ds / dn)), digits = 3)

    #  method = :minlike same results R, python two.sides
    output["pvalue"] = pvalue(FisherExactTest(Int(ps), Int(ceil(pn_neutral)), Int(ds),
                                              Int(dn)))

    if (!isnothing(m))
        mn = m[1]
        ms = m[2]
        # ## Estimation of b: weakly deleterious
        output["b"] = (deleterious / ps) * (ms / mn)

        ## Estimation of f: neutral sites
        output["f"] = (ms * pn_neutral) / (mn * ps)

        ## Estimation of d, strongly deleterious sites
        output["d"] = 1 - (output["f"] + output["b"])

        ka = dn / mn
        ks = ds / ms
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
 - `sfs::Matrix{Float64}`: SFS data parsed from parse_sfs().
 - `divergence::Matrix{Int64}`: divergence data from parse_sfs().
 - `cutoff{Float64}`: frequency cutoff to perform imputedMK.
 - `m::Matrix{Int64}`: total number of sites parsed from parse_sfs()
# Returns
 - `Dict: Dictionary containing imputedMK estimation.
"""
function fwwMK(param::parameters, sfs::Matrix{Float64}, divergence::Matrix{Int64};
               cutoff::Float64 = 0.15, m::T = nothing) where {T <: Union{Nothing, Array}}
    output = OrderedDict{String, Float64}()

    ps = sum(sfs[:, 2])
    pn = sum(sfs[:, 3])
    dn = divergence[1]
    ds = divergence[2]

    deleterious = 0
    ### Estimating slightly deleterious with pn/ps ratio
    sfs_tmp = deepcopy(sfs)
    sfs_tmp[:, 1] = sfs_tmp[:, 1] ./ param.nn
    flt_inter = (sfs_tmp[:, 1] .>= cutoff) .& (sfs_tmp[:, 1] .<= 1)
    pn_inter = sum(sfs_tmp[flt_inter, 2])
    ps_inter = sum(sfs_tmp[flt_inter, 3])

    # output['alpha'] = 1 - (((pn - deleterious) / ps) * (ds / dn))
    output["alpha"] = round(1 - ((pn_inter / ps_inter) * (ds / dn)), digits = 3)
    #  method = :minlike same results R, python two.sides
    output["pvalue"] = pvalue(FisherExactTest(Int(ps_inter), Int(ceil(pn_inter)), Int(ds),
                                              Int(dn)))

    if (!isnothing(m))
        mn = m[1]
        ms = m[2]
        ka = dn / mn
        ks = ds / ms
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
function standardMK(sfs::Matrix{Float64}, divergence::Matrix{Int64};
                    m::T = nothing) where {T <: Union{Nothing, Array}}
    output = OrderedDict{String, Float64}()

    pn = sum(sfs[:, 2])
    ps = sum(sfs[:, 3])
    dn = divergence[1]
    ds = divergence[2]

    output["alpha"] = round(1 - ((pn / ps) * (ds / dn)), digits = 3)
    #  method = :mnnlike same results R, python two.sides
    output["pvalue"] = pvalue(FisherExactTest(Int(ps), Int(ceil(pn)), Int(ds), Int(dn)))

    if (!isnothing(m))
        mn = m[1]
        ms = m[2]

        ka = dn / mn
        ks = ds / ms
        output["omega"] = ka / ks

        # Omega A and Omega D
        output["omegaA"] = output["omega"] * output["alpha"]
        output["omegaD"] = output["omega"] - output["omegaA"]
    end

    return output
end

"""
    grapes(sfs,divergence,m,model,folder,bins)

Run grapes using SFS and divergence data. The function will install grapes using Conda from genomedk channel. 

# Arguments
 - `sfs::Vector{Matrix{Float64}}`: SFS data parsed from MKtest.parse_sfs()
 - `divergence::Vector{Matrix{Int64}}`: divergence data parsed from MKtest.parse_sfs()
 - `m::Vector{Matrix{Int64}}`: total number of sites parsed from MKtest.parse_sfs()
 - `model::String`: grapes model.
 - `folder::String`: output folder.
 - `bins::Int64`:  bin size to reduce the SFS. Grapes became unstable when inputing large SFS.

# Output
 - `DataFrame`: grapes model estimation.
"""
function grapes(sfs::Vector{Matrix{Float64}},
                divergence::Vector{Matrix{Int64}},
                m::Vector{Matrix{Int64}},
                model::String,
                folder::String,
                bins::Int64)

    grapes_bin = CondaPkg.which("grapes")

    if (isnothing(grapes_bin))
        CondaPkg.add("grapes-static", channel = "genomedk")
    end

    @assert model ∈ ["GammaZero", "GammaExpo", "DisplGamma", "ScaledBeta", "FGMBesselK"] "Please select a valid model: GammaZero GammaExpo DisplGamma ScaledBeta FGMBesselK"

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
    r(d, o, m = model, gr = grapes_bin) = run(`$gr -in $d -out $o -model $m`)

    @suppress_out begin ThreadsX.map(r, output_dofe, output_grapes) end

    df = @suppress begin CSV.read.(output_grapes, DataFrame, footerskip = 1, skipto = 3) end

    return (vcat(df...))
end

"""
    α_x(sfs,divergence)

Function to estimate α(x).

# Arguments
 - `sfs::Matrix`: SFS data Matrix parsed from parse_sfs().
 - `divergence::Matrix`: divergence data Matrix from parse_sfs().
 - `cutoff::Vector{Float64}`: frequency cutoff to perform α(x) estimation.
 - `cumulative::Bool`: perform the estimation using the cumulative SFS.

# Returns
 - `Vector{Float64}: α(x) estimations.
"""
function α_x(sfs::Matrix, divergence::Matrix; cutoff::Vector{Float64}=[0.0,1.0],cumulative::Bool = true)

    sfs_tmp = copy(sfs)

    n = size(sfs,1) + 1
    sfs_tmp[:,1] .= collect(1:n-1) / n
    if cumulative
        sfs_tmp = sfs_tmp[(sfs_tmp[:,1] .>= cutoff[1]) .& (sfs_tmp[:,1] .<= cutoff[2]),:]
        sfs_tmp = cumulative_sfs(sfs_tmp)
    else
        sfs_tmp = sfs_tmp[(sfs_tmp[:,1] .>= cutoff[1]) .& (sfs_tmp[:,1] .<= cutoff[2]),:]
    end

    alpha_x = @. 1 - (sfs_tmp[:, 2] / sfs_tmp[:, 3] * divergence[2] / divergence[1])

    return (alpha_x)
end

"""
    α_x(sfs,divergence)

Function to estimate α(x).

# Arguments
 - `sfs::Vector`: Multiple SFS data Matrix parsed from parse_sfs().
 - `divergence::Vector`: Multiple divergence data Matrix from parse_sfs().
 - `cutoff::Vector{Float64}`: frequency cutoff to perform α(x) estimation.
 - `cumulative::Bool`: perform the estimation using the cumulative SFS.
# Returns
 - `Vector{Float64}: α(x) estimations.
"""
function α_x(sfs::Vector, divergence::Vector; cutoff::Vector{Float64}=[0.0,1.0],cumulative::Bool = true)

    output = Vector{Vector}(undef, length(sfs))
    sfs_tmp = deepcopy(sfs)

    Threads.@threads for i in eachindex(sfs)
        @show
        n = size(sfs_tmp[i],1) + 1
        sfs_tmp[i][:,1] .= collect(1:n-1)/n

        if cumulative
            sfs_tmp[i] = sfs_tmp[i][(sfs_tmp[i][:,1] .>= cutoff[1]) .& (sfs_tmp[i][:,1] .<= cutoff[2]),:]
            sfs_tmp[i] .= cumulative_sfs(sfs_tmp[i])
        else
            sfs_tmp[i] = sfs_tmp[i][(sfs_tmp[i][:,1] .>= cutoff[1]) .& (sfs_tmp[i][:,1] .<= cutoff[2]),:]
        end

        output[i] = @. 1 - (sfs_tmp[i][:, 2] / sfs_tmp[i][:, 3] * divergence[i][2] / divergence[i][1])
    end
    return (output)
end

"""
    aMK(param,α;dac_rm,na_rm)

Function to estimate the asymptotic value of α(x).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model.
 - `sfs::Matrix`: SFS data parsed from parse_sfs().
 - `divergence::Matrix`: divergence data from parse_sfs().
 - `threshold{Float64}`: frequency cutoff to perform imputedMK.
 - `cumulative::Bool`: Perform α(x) estimation using the cumulative SFS
 - `dac_rm::Bool`: Remove any value not selected at param.dac.
 - `na_rm::Bool`: Remove any NA value at α(x)

# Returns
 - `DataFrame`
"""
function aMK(
    param::parameters,
    sfs::Matrix,
    divergence::Matrix;
    cumulative::Bool=true,
    dac_rm::Bool = false,
    na_rm::Bool = false,
)


    @unpack nn, dac, cutoff = param

    f = collect(1:size(sfs,1)) / nn
    f = f[(f.>=cutoff[1]) .& (f.<=cutoff[2])]

    α_i = α_x(sfs,divergence,cutoff=cutoff,cumulative=cumulative)

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

    a,b,c = fitted2.param
    return DataFrame(:amk=>asymp,:ci_low=>ci[1],:ci_high=>ci[2], :a=>a, :b=>b, :c=>c)
end


"""
    aMK(param,α;dac_rm,na_rm)

Function to estimate the asymptotic value of α(x).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model.
 - `sfs::Vector`: SFS data parsed from parse_sfs().
 - `divergence::Vector`: divergence data from parse_sfs().
 - `threshold{Float64}`: frequency cutoff to perform imputedMK.
 - `cumulative::Bool`: Perform α(x) estimation using the cumulative SFS
 - `dac_rm::Bool`: Remove any value not selected at param.dac.
 - `na_rm::Bool`: Remove any NA value at α(x)

# Returns
 - `DataFrame`
"""
function aMK(
    param::parameters,
    sfs::Vector,
    divergence::Vector;
    cumulative::Bool=true,
    dac_rm::Bool = false,
    na_rm::Bool = false
)

    output = Vector{Tuple}(undef, length(sfs))
    f = Vector{Vector{Float64}}(undef, length(sfs))

    @unpack nn, dac, cutoff = param

    α_i = α_x(sfs,divergence,cutoff=cutoff,cumulative=cumulative)

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
    amk_v = Vector{Float64}(undef, length(sfs))
    model_params = Vector{Vector{Float64}}(undef, length(sfs))
    ci_v = Vector{Tuple{Float64,Float64}}(undef,length(sfs))
    fitted1 = ThreadsX.mapi((x,y) -> curve_fit(model,x,y,[-1.0, -1.0, 1.0];lower = [-1.0, -1.0, 1.0],upper = [1.0, 1.0, 10.0]),f,α_i)
    fitted2 = ThreadsX.mapi((x,y,z) -> curve_fit(model,x,y,z.param),f,α_i,fitted1)

    for (i,v) ∈ enumerate(fitted2)
        amk_v[i] = model(1, v.param)
        model_params[i] = v.param
        try
            ci_v[i] = confidence_interval(v)[1]
        catch
            [0.0,0.0]
        end
    end

    return DataFrame(hcat(amk_v,hcat(collect.(ci_v)...)',hcat(model_params...)'),[:amk,:ci_low,:ci_high,:a,:b,:c])
end

"""
    imputedMK(param,sfs,divergence;m,cutoff)

Function to estimate the imputedMK (https://doi.org/10.1093/g3journal/jkac206).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model.
 - `sfs::Vector`: SFS data parsed from parse_sfs().
 - `divergence::Vector`: divergence data from parse_sfs().
 - `threshold{Float64}`: frequency cutoff to perform imputedMK.

# Returns
 - `DataFrame`
"""

function imputedMK(
    param::parameters,
    sfs::Vector,
    divergence::Vector;
    threshold::Float64 = 0.15
)

    alpha   = Vector{Float64}(undef, length(sfs))
    p_value = Vector{Float64}(undef, length(sfs))

    b       = Vector{Float64}(undef, length(sfs))
    f       = Vector{Float64}(undef, length(sfs))
    d       = Vector{Float64}(undef, length(sfs))
    omega   = Vector{Float64}(undef, length(sfs))
    omega_a = Vector{Float64}(undef, length(sfs))
    omega_d = Vector{Float64}(undef, length(sfs))

    @unpack isolines, n, nn, cutoff = param
    sfs_tmp = deepcopy(sfs)

    for i ∈ eachindex(sfs)
        if(sfs_tmp[i][1,1] >= 1)
            s_size = ifelse(isolines,n,nn)
            freqs = collect(1:s_size-1) ./ s_size
            sfs_tmp[i][:,1] .= freqs[(freqs .>= cutoff[1]) .& (freqs .<= cutoff[2])]
        end

        pn = sum(sfs_tmp[i][:, 2])
        ps = sum(sfs_tmp[i][:, 3])
        dn = divergence[i][1]
        ds = divergence[i][2]

        deleterious = 0

        flt_low = (sfs_tmp[i][:, 1] .<= threshold)
        pn_low = sum(sfs_tmp[i][flt_low, 2])
        ps_low = sum(sfs_tmp[i][flt_low, 3])

        flt_inter = (sfs_tmp[i][:, 1] .> threshold) .& (sfs_tmp[i][:, 1] .<= 1)
        pn_inter = sum(sfs_tmp[i][flt_inter, 2])
        ps_inter = sum(sfs_tmp[i][flt_inter, 3])

        ratio_ps = ps_low / ps_inter
        deleterious = pn_low - (pn_inter * ratio_ps)

        if (deleterious > pn) || (deleterious < 0)
            deleterious = 0
            pn_neutral = round(pn - deleterious, digits = 3)
        else
            pn_neutral = round(pn - deleterious, digits = 3)
        end

        alpha[i] = round(1 - ((pn_neutral / ps) * (ds / dn)), digits = 3)

        #  method = :minlike same results R, python two.sides
        p_value[i] =
            pvalue(FisherExactTest(Int(ceil(ps)), Int(ceil(pn_neutral)), Int(ds), Int(dn)))

        if size(divergence[i],2) > 2
            ln = divergence[i][3]
            ls = divergence[i][4]

            # ## Estimation of b: weakly deleterious
            b[i] = (deleterious / ps) * (ls / ln)

            ## Estimation of f: neutral sites
            f[i] = (ls * pn_neutral) / (ln * ps)

            ## Estimation of d, strongly deleterious sites
            d[i] = 1 - (f[i] + b[i])

            ka = dn / ln
            ks = ds / ls
            omega[i] = ka / ks

            # omega_a and omega_d
            omega_a[i] = omega[i] * alpha[i]
            omega_d[i] = omega[i] - omega_a[i]
        end
    end

    return DataFrame(:alpha => alpha, :p_value => p_value, :b=>b, :f=>f, :d=>d, :omega=>omega, :omega_a=>omega_a, :omega_d=>omega_d)
end

"""
    imputedMK(param,sfs,divergence;m,cutoff)

Function to estimate the imputedMK (https://doi.org/10.1093/g3journal/jkac206).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model.
 - `sfs::Matrix{Float64}`: SFS data parsed from parse_sfs().
 - `divergence::Matrix{Int64}`: divergence data from parse_sfs().
 - `threshold{Float64}`: frequency cutoff to perform imputedMK.

# Returns
 - `DataFrame`
"""

function imputedMK(
    param::parameters,
    sfs::Matrix,
    divergence::Matrix;
    threshold::Float64 = 0.15,
)

    @unpack isolines,n,nn,cutoff= param

    sfs_tmp = deepcopy(sfs)
    if(sfs_tmp[1,1] >= 1)
        s_size = ifelse(isolines,n,nn)
        freqs = collect(1:s_size-1) ./ s_size
        sfs_tmp[:,1] .= freqs[(freqs .>= cutoff[1]) .& (freqs .<= cutoff[2])]
    end

    pn = sum(sfs_tmp[:, 2])
    ps = sum(sfs_tmp[:, 3])
    dn = divergence[1]
    ds = divergence[2]

    deleterious = 0
    ### Estimating slightly deleterious with pn/ps ratio
    flt_low = (sfs_tmp[:, 1] .<= threshold)
    pn_low = sum(sfs_tmp[flt_low, 2])
    ps_low = sum(sfs_tmp[flt_low, 3])

    flt_inter = (sfs_tmp[:, 1] .>= threshold) .& (sfs_tmp[:, 1] .<= 1)
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

    alpha = round(1 - ((pn_neutral / ps) * (ds / dn)), digits = 3)

    #  method = :minlike same results R, python two.sides
    p_value =
        pvalue(FisherExactTest(Int(ceil(ps)), Int(ceil(pn_neutral)), Int(ds), Int(dn)))

    b = f = d = omega = omega_a = omega_d = NaN;
    if size(divergence,2) > 2
        ln = divergence[3]
        ls = divergence[4]
        # ## Estimation of b: weakly deleterious
        b = (deleterious / ps) * (ls / ln)

        ## Estimation of f: neutral sites
        f = (ls * pn_neutral) / (ln * ps)

        ## Estimation of d, strongly deleterious sites
        d = 1 - (f + b)

        ka = dn / ln
        ks = ds / ls
        omega = ka / ks

        # omega_a and omega_d
        omega_a = omega * alpha
        omega_d = omega - omega_a
    end

    return DataFrame(:alpha => alpha, :p_value => p_value, :b=>b, :f=>f, :d=>d, :omega=>omega, :omega_a=>omega_a, :omega_d=>omega_d)
end

"""
    fwwMK(param,sfs,divergence;m,cutoff)

Function to estimate the fwwMK (https://doi.org/10.1038/4151024a).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model.
 - `sfs::Vector`: SFS data parsed from parse_sfs().
 - `divergence::Vector`: divergence data from parse_sfs().
 - `threshold{Float64}`: frequency cutoff to perform fwwMK.
# Returns
 - `DataFrame`
"""
function fwwMK(
    param::parameters,
    sfs::Vector,
    divergence::Vector;
    threshold::Float64 = 0.15,
)

    alpha   = Vector{Float64}(undef, length(sfs))
    p_value = Vector{Float64}(undef, length(sfs))

    omega   = Vector{Float64}(undef, length(sfs))
    omega_a = Vector{Float64}(undef, length(sfs))
    omega_d = Vector{Float64}(undef, length(sfs))

    @unpack isolines,n,nn,cutoff = param;
    sfs_tmp = deepcopy(sfs)

    for i ∈ eachindex(sfs_tmp)
        if(sfs_tmp[i][1,1] >= 1)
            s_size = ifelse(isolines,n,nn)
            freqs = collect(1:s_size-1) ./ s_size
            sfs_tmp[i][:,1] .= freqs[(freqs .>= cutoff[1]) .& (freqs .<= cutoff[2])]
        end

        ps = sum(sfs_tmp[i][:, 2])
        pn = sum(sfs_tmp[i][:, 3])
        dn = divergence[i][1]
        ds = divergence[i][2]

        ### Estimating slightly deleterious with pn/ps ratio
        flt_inter = (sfs_tmp[i][:, 1] .>= threshold) .& (sfs_tmp[i][:, 1] .<= 1)
        pn_inter = sum(sfs_tmp[i][flt_inter, 2])
        ps_inter = sum(sfs_tmp[i][flt_inter, 3])


        alpha[i] = round(1 - ((pn_inter / ps_inter) * (ds / dn)), digits = 3)
        #  method = :minlike same results R, python two.sides
        p_value[i] =
            pvalue(FisherExactTest(Int(ceil(ps_inter)), Int(ceil(pn_inter)), Int(ds), Int(dn)))

        if size(divergence[i],2) > 2
            ln = divergence[i][3]
            ls = divergence[i][4]
            ka = dn / ln
            ks = ds / ls

            omega[i] = ka / ks
            # omega_a and omega_d
            omega_a[i] = omega[i] * alpha[i]
            omega_d[i] = omega[i] - omega_a[i]
        end
    end

    return DataFrame(:alpha => alpha, :p_value => p_value, :omega => omega, :omega_a => omega_a, :omega_d => omega_d)
end

"""
    fwwMK(param,sfs,divergence;m,cutoff)

Function to estimate the fwwMK (https://doi.org/10.1038/4151024a).

# Arguments
 - `param::parameters`: mutable structure containing the variables required to solve the model.
 - `sfs::Matrix`: SFS data parsed from parse_sfs().
 - `divergence::Matrix`: divergence data from parse_sfs().
 - `threshold{Float64}`: frequency cutoff to perform fwwMK.

# Returns
 - `DataFrame`
"""
function fwwMK(
    param::parameters,
    sfs::Matrix,
    divergence::Matrix;
    threshold::Float64 = 0.15,
)

    @unpack isolines,n,nn,cutoff= param;

    # DAC to freqs
    if(sfs_tmp[1,1] >= 1)
        s_size = ifelse(isolines,n,nn)
        freqs = collect(1:s_size-1) ./ s_size
        sfs_tmp[:,1] .= freqs[(freqs .>= cutoff[1]) .& (freqs .<= cutoff[2])]
    end

    ps = sum(sfs_tmp[:, 2])
    pn = sum(sfs_tmp[:, 3])
    dn = divergence[1]
    ds = divergence[2]

    deleterious = 0
    ### Estimating slightly deleterious with pn/ps ratio
    s_size = ifelse(isolines,n,nn)

    sfs_tmp[:, 1] = sfs_tmp[:, 1] ./ s_size
    flt_inter = (sfs_tmp[:, 1] .>= threshold) .& (sfs_tmp[:, 1] .<= 1)
    pn_inter = sum(sfs_tmp[flt_inter, 2])
    ps_inter = sum(sfs_tmp[flt_inter, 3])


    alpha = round(1 - ((pn_inter / ps_inter) * (ds / dn)), digits = 3)
    #  method = :minlike same results R, python two.sides
    p_value =
        pvalue(FisherExactTest(Int(ps_inter), Int(ceil(pn_inter)), Int(ds), Int(dn)))

    omega = omega_a = omega_d = 0
    if size(divergence,2) > 2
        ln = divergence[3]
        ls = divergence[4]
        ka = dn / ln
        ks = ds / ls

        omega = ka / ks

        # omega_a and omega_d
        omega_a = omega * alpha
        omega_d = omega - omega_a
    end

    return  DataFrame(:alpha => alpha, :p_value => p_value, :omega=>omega, :omega_a=>omega_a, :omega_d=>omega_d)
end

"""
	standardMK(sfs,divergence;m)

Function to estimate the original α value.

# Arguments
 - `sfs::Matrix{Float64}`: SFS array
 - `divergence::Array`: divegence count array

# Output
 - `DataFrame`
"""
function standardMK(
    sfs::Vector,
    divergence::Vector;
)

    alpha   = Vector{Float64}(undef, length(sfs))
    p_value = Vector{Float64}(undef, length(sfs))

    omega   = Vector{Float64}(undef, length(sfs))
    omega_a = Vector{Float64}(undef, length(sfs))
    omega_d = Vector{Float64}(undef, length(sfs))
    sfs_tmp = deepcopy(sfs)

    for i ∈ eachindex(sfs)
        # DAC to freqs
        if(sfs[i][1,1] > 1)
            nn = size(view(sfs[i],:,1),1)
            sfs[i][:,1] .= collect(1:nn) ./ (nn + 1)
        end
        pn = sum(sfs_tmp[i][:, 2])
        ps = sum(sfs_tmp[i][:, 3])
        dn = divergence[i][1]
        ds = divergence[i][2]

        alpha[i] = round(1 - ((pn / ps) * (ds / dn)), digits = 3)
        #  method = :lnnlike same results R, python two.sides
        p_value = pvalue(FisherExactTest(Int(ceil(ps)), Int(ceil(pn)), Int(ds), Int(dn)))

        if size(divergence[i],2) > 2
            ln = divergence[i][3]
            ls = divergence[i][4]

            ka = dn / ln
            ks = ds / ls
            omega[i] = ka / ks

            # omega_a and omega_d
            omega_a[i] = omega[i] * alpha[i]
            omega_d[i] = omega[i] - omega_a[i]
        end
    end
    return DataFrame(:alpha => alpha, :p_value => p_value, :omega=>omega, :omega_a=>omega_a, :omega_d=>omega_d)
end

"""
    standardMK(sfs,divergence;m)

Function to estimate the original α value.

# Arguments
 - `sfs::Matrix{Float64}`: SFS array
 - `divergence::Array`: divegence count array

# Output
 - `DataFrame`
"""
function standardMK(
    sfs::Matrix,
    divergence::Matrix;
)

    # DAC to freqs
    if(sfs[1,1] >= 1)
        nn = size(view(sfs,:,1),1)
        sfs[:,1] .= collect(1:nn) ./ (nn + 1)
    end

    pn = sum(sfs[:, 2])
    ps = sum(sfs[:, 3])
    dn = divergence[1]
    ds = divergence[2]

    alpha = round(1 - ((pn / ps) * (ds / dn)), digits = 3)
    #  method = :lnnlike same results R, python two.sides
    p_value = pvalue(FisherExactTest(Int(ps), Int(ceil(pn)), Int(ds), Int(dn)))

    if size(divergence,2) > 2
        ln = divergence[3]
        ls = divergence[4]

        ka = dn / ln
        ks = ds / ls
        omega = ka / ks

        # omega_a and omega_d
        omega_a = omega * alpha
        omega_d = omega - omega_a
    end
    return DataFrame(:alpha => alpha, :p_value => p_value, :omega=>omega, :omega_a=>omega_a, :omega_d=>omega_d)
end

"""
    grapes(sfs,divergence,model)

Run grapes using SFS and divergence data. The function will install grapes using Conda from genomedk channel. 

# Arguments
 - `sfs::Vector{Matrix{Float64}}`: SFS data parsed from parse_sfs()
 - `divergence::Vector{Matrix{Int64}}`: divergence data parsed from parse_sfs()
 - `model::String`: grapes model (GammaZero, GammaExpo, DisplGamma, ScaledBeta, FGMBesselK and all)
 - `n::Int64`: smaller sample size to project the SFS
 - `nearly_neutral::Int64`: Sadp threshold of above which a mutation is considered adaptive
 - `FWW_threshold::Float64`: minimal allele frequency in FWW α
 - `nb_rand_start::Int64`: number of random starting values in model optimization (default=0); setting positive values will slow down the program but decrease the probability of being trapped in local optima.
 - `anc_to_rec_Ne_ratio::Float64`: divergence/polymorphism Ne ratio
 - `no_div_data::Bool`: only use divergence data to estimate DFE
 - `no_div_param::Bool`: implements so-called [-A] version in Galtier (2016); also called a DFE by Tataru et al. (2017), and indicated with * in Rousselle et al. (2018); irrelevant if model = GammaZero; (default=false)
 - `no_syn_orient_error::Bool`: force equal synonymous and non-synonymous mis-orientation rate(default=false)
 - `fold::Bool`: fold the SFS
 - `fixed_param::String`: this option should be used if one does not want to optimize every parameter, but rather have some parameters fixed to predefined values; parameter names and predefined values are passed via a control file; see example at the bottom of this file (default=none).
# Output
 - `DataFrame`: grapes model estimation.
"""
function grapes(
    sfs::Union{Vector,Matrix},
    divergence::Union{Vector,Matrix},
    model::String;
    n::Int64=0,
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

    # Change to Vectors of Matrix if needed
    if sfs isa Matrix
        sfs = [sfs]
        divergence = [divergence]
    end

    # Reduce to n-1 to input in Grapes as SFS = n - 1 entries
    if n != 0
        sfs_r = project.(sfs, n)
    else
        sfs_r = sfs
        n     = size.(sfs,1) .+ 1
    end

    pn = map(x -> permutedims(x[:, 2]), sfs_r)
    ps = map(x -> permutedims(x[:, 3]), sfs_r)

    dn = map(x -> [x[1]], divergence)
    ds = map(x -> [x[2]], divergence)
    ln = map(x -> [x[3]], divergence)
    ls = map(x -> [x[4]], divergence)

    idx = string.(collect(1:length(sfs_r)))

    @info "Converting SFS to dofe file"
    dofe = @. sfs_to_dofe(pn, ps, dn, ds, ln, ls, idx, n)
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

function sfs_to_dofe(pn, ps, dn, ds, ln, ls, w, n)
    return(DataFrame(hcat("dofe_" * string(w), n, ln, pn, ls, ps, ln, dn, ls, ds...), :auto))
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

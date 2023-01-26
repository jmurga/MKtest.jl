function simulate_models!(param::parameters,
                          binom::Dict{Float64, SparseMatrixCSC{Float64, Int64}},
                          var_params::Vector{Float64})
    (gH, gL, gam_flanking, gam_dfe, shape, alpha, alpha_low) = var_params

    param.shape = shape
    param.gam_dfe = gam_dfe
    param.scale = abs(shape / gam_dfe)
    param.gH = gH
    param.gL = gL
    param.gam_flanking = gam_flanking
    param.gam_dfe = gam_dfe
    param.al_tot = alpha
    param.al_low = alpha * alpha_low

    m, r_ps, r_pn, r_f = solve_model!(param, binom)

    return (m, r_ps, r_pn, r_f)
end

function solve_model!(param::parameters,
                      binom::Dict{Float64, SparseMatrixCSC{Float64, Int64}})
    @unpack B, dac, B_bins = param

    assertion_params(param)

    m = zeros(size(B_bins, 1), 8)
    r_ps = zeros(size(B_bins, 1), size(dac, 1))
    r_pn = zeros(size(B_bins, 1), size(dac, 1))
    r_f = zeros(size(B_bins, 1), 4)

    # Solving θ on non-coding region and probabilites to get α value without BGS
    param.B = 0.999
    set_θ!(param)
    try
        set_ppos!(param)
        for j::Int64 in 1:length(B_bins)
            # Set B value
            param.B = param.B_bins[j]
            # Solve θ non-coding for the B value.
            set_θ!(param)
            # Solve model for the B value
            x, y, z, w = try
                solve_rates(param, binom[param.B])
            catch
                zeros(size(m, 2)),
                zeros(size(r_ps, 2)),
                zeros(size(r_pn, 2)),
                zeros(size(r_pf, 2))
            end
            @inbounds m[j, :] = x
            @inbounds r_ps[j, :] = y
            @inbounds r_pn[j, :] = z
            @inbounds r_f[j, :] = w
        end
        return (m, r_ps, r_pn, r_f)
    catch
        return (m, r_ps, r_pn, r_f)
    end

    return (m, r_ps, r_pn, r_f)
end

function solve_rates(param::parameters, binom::SparseMatrixCSC{Float64, Int64})

    ################################################
    # Subset rates accounting for positive alleles #
    ################################################
    @unpack B,
    ppos_h,
    ppos_l,
    gL,
    gH,
    B,
    al_low,
    al_tot,
    gam_flanking,
    gL,
    gH,
    shape,
    scale,
    dac = param

    # Fixation
    f_n = B * fix_neut(param)
    f_neg = B * fix_neg(param, 0.5 * ppos_h + 0.5 * ppos_l)
    f_pos_l = fix_pos_sim(param, gL, 0.5 * ppos_l)
    f_pos_h = fix_pos_sim(param, gH, 0.5 * ppos_h)

    ds = f_n
    dn = f_neg + f_pos_l + f_pos_h

    ## Polymorphism
    neut::Vector{Float64} = sfs_neut(param, binom)
    sel_h::Vector{Float64} = if isinf(exp(gH * 2))
        sfs_pos_float(param, gH, ppos_h, binom)
    else
        sfs_pos(param, gH, ppos_h, binom)
    end

    sel_l::Vector{Float64} = sfs_pos(param, gL, ppos_l, binom)
    sel_neg::Vector{Float64} = sfs_neg(param, ppos_l + ppos_h, binom)

    split_columns(matrix::Matrix{Float64}) = (view(matrix, :, i) for i in 1:size(matrix, 2))
    tmp = cumulative_sfs(hcat(neut, sel_h, sel_l, sel_neg), false)

    neut, sel_h, sel_l, sel_neg = split_columns(tmp)
    sel = (sel_h + sel_l) + sel_neg

    ##########
    # Output #
    ##########
    analytical_m::Matrix{Float64} = hcat(B, al_low, al_tot, gam_flanking, gL, gH, shape,
                                         -(shape / scale))
    analytical_ps::Matrix{Float64} = permutedims(neut[dac])
    analytical_pn::Matrix{Float64} = permutedims(sel[dac])
    analytical_f::Matrix{Float64} = hcat(ds, dn, f_pos_l, f_pos_h)

    return (analytical_m, analytical_ps, analytical_pn, analytical_f)
end

"""
	rates(param,gH,gL,gam_flanking,gam_dfe,alpha,iterations,output)

Function to solve randomly *N* scenarios. The function will create *N* models, defined by ```MKtest.parameters()```, to solve analytically fixation rates and the expected SFS for each model. The rates will be used to compute summary statistics required at ABC inference. The function output a HDF5 file containing the solved models, the selected DAC and the analytical solutions. 

# Arguments
 - `param::parameters`: mutable structure containing the model.
 - `gH::Array{Int64,1}`: Range of strong selection coefficients.
 - `gL::Union{Array{Int64,1},Nothing}`: Range of weak selection coefficients.
 - `gam_flanking::Array{Int64,1}`: Range of deleterious selection coefficients at the flanking region.
 - `gam_dfe::Array{Int64,1}`: Range of deleterious selection coefficients at the coding region.
 - `alpha::Vector{Float64}`: Range of α value to solve.
 - `iterations::Int64`: Number of solutions.
 - `output::String`: File to output HDF5 file.
# Returns
 - `DataFrame`: models solved.
 - `Output`: HDF5 file containing models solved and rates.
"""
function rates(param::parameters;
               gH::Vector{Int64},
               gL::Vector{Int64},
               gam_flanking::Vector{Int64},
               gam_dfe::Vector{Int64},
               alpha::Vector{Float64} = [0.1, 0.9],
               iterations::Int64,
               output::String)

    # param = parameters(N=1000,n=661);gH=[200,2000];gL=[1,10];gam_flanking=[-1000.,-100];gam_dfe=[-1000.,-100];iterations = 100;alpha=[0.1,0.9];B=[0.1,0.999];iterations=100000

    @info "Solving binomial convolution to downsample the SFS"
    binom = binom_op(param)

    assertion_params(param)

    @unpack shape = param

    @info "Creating priors distributions"

    # Parameters ranges
    u_gh = gH[1]:gH[end]
    u_gl = gL[1]:gL[end]
    u_gam_flanking = gam_flanking[1]:gam_flanking[end]
    u_gam_dfe = gam_dfe[1]:gam_dfe[end]
    afac = -2:0.05:2
    u_tot = alpha[1]:0.01:alpha[end]
    u_low = 0.0:0.05:0.9

    # Priors
    priors = Vector{Float64}[]
    for i::Int64 in 1:iterations
        push!(priors,
              [
                  rand(u_gh),
                  rand(u_gl),
                  rand(u_gam_flanking),
                  rand(u_gam_dfe),
                  shape * 2^rand(afac),
                  rand(u_tot),
                  rand(u_low),
              ])
    end

    # Solving models in multi-threading
    @info "Solving models in multiple threads"
    m, r_ps, r_pn, r_f = unzip(ThreadsX.map(x -> simulate_models!(deepcopy(param), binom, x),
                                            priors))

    @info "Saving solved models in $output"
    # Reducing models array
    df = vcat(m...)
    # Filtering cases where deleterious DFE cannot be solved.
    idx = sum.(eachrow(df)) .!= 0

    # Saving models and rates
    models = @view DataFrame(df,
                             [
                                 :B,
                                 :al_low,
                                 :al_tot,
                                 :gam_flanking,
                                 :gL,
                                 :gH,
                                 :shape,
                                 :gam_dfe,
                             ])[idx,
                                :]
    neut = @view vcat(r_ps...)[idx, :]
    sel = @view vcat(r_pn...)[idx, :]
    dsdn = @view vcat(r_f...)[idx, :]

    # Saving multiple summary statistics
    n = OrderedDict{Int, Array}()
    s = OrderedDict{Int, Array}()
    for i in eachindex(param.dac)
        n[param.dac[i]] = neut[:, i]
        s[param.dac[i]] = sel[:, i]
    end

    # Writting HDF5 file
    string_cutoff = "cutoff=[" * string(param.cutoff[1]) * "," * string(param.cutoff[end]) *
                    "]"
    JLD2.jldopen(output, "a+") do file
        file[string(param.N) * "/" * string(param.n) * "/" * string_cutoff * "/models"] = models
        file[string(param.N) * "/" * string(param.n) * "/" * string_cutoff * "/neut"] = n
        file[string(param.N) * "/" * string(param.n) * "/" * string_cutoff * "/sel"] = s
        file[string(param.N) * "/" * string(param.n) * "/" * string_cutoff * "/dsdn"] = dsdn
        file[string(param.N) * "/" * string(param.n) * "/" * string_cutoff * "/dac"] = param.dac
    end

    return(Dict(:models => models,:neut => n,:sel => s,:dsdn => dsdn,:dax => param.dac))
end

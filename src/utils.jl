function cumulative_vector!(v::Union{SubArray,Vector{T}}) where {T<:Number}
    @inbounds for i = length(v)-1:-1:1
        v[i] += v[i+1]
    end
    return v
end

"""
	cumulative_sfs(sfs_tmp)

Changing SFS considering all values above a frequency *x*. The original aMK approach takes Pn(x) and Ps(x) as the number of polymorphic sites at frequency *x* rather than above *x*, but this approach scales poorly as sample size increases. We define the polymorphic spectrum as stated above since these quantities trivially have the same asymptote but are less affected by changing sample size.

# Arguments
 - `sfs_tmp::Union{Array,SubArray}`: SFS vector.
 - `freqs::Bool`: true/false wheter the input SFS vector containing a first columns with the frequencies.
# Output 
 - `Vector{Float64}`: SFS vector.
"""
function cumulative_sfs(sfs_tmp::Union{Array, SubArray}, freqs::Bool = true)
    n, m = size(sfs_tmp)
    out = similar(sfs_tmp)

    idx = ifelse(freqs, 2, 1)

    out[1, idx:end] = sum(sfs_tmp[:, idx:end], dims=1)

    @inbounds @simd for i in 2:n
        app = view(out,i-1, idx:size(out,2)) .- view(sfs_tmp,i-1, idx:size(out,2))

        if sum(app) > 0.0
            out[i, idx:end] .= app
        else
            out[i, idx:end] .= 0.0
        end
    end

    if freqs
        out[:, 1] .= sfs_tmp[:, 1]
    end

    return out
end


"""
	reduce_sfs(sfs_tmp,bins)

Function to reduce the SFS into N bins.

# Arguments
 - `sfs_tmp::Vector{Float64}`: SFS vector.
 - `bins::Int64`: bin size to collapse the SFS.
# Output 
 - `Vector{Float64}`: Binned SFS vector.
"""
function reduce_sfs(sfs_tmp::Array, bins::Int64)
    freq = collect(0:(size(sfs_tmp, 1)-1)) / size(sfs_tmp, 1)
    h1 = StatsBase.fit(StatsBase.Histogram, freq, 0:(1/(bins-1)):1)
    xmap1 = StatsBase.binindex.(Ref(h1), freq)

    tmp = hcat(sfs_tmp, xmap1)
    out = Array{Float64}(undef, bins - 1, size(sfs_tmp, 2))
    vIter = convert(Array, unique(xmap1)')
    @simd for i in eachindex(vIter)
        @inbounds out[i, 2:end] = sum(tmp[tmp[:, end].==i, 2:(end-1)], dims = 1)
    end

    out[:, 1] = collect(1:(bins-1)) ./ bins

    return (out)
end

"""
Log of N choose k.

# Arguments
 - `N::Int64`
 - `k::Union{Vector{Int64},Int64}`
# Output
 - `Vector{Float64}`
"""
function _lncomb(N::Int64,k::Union{Vector{Int64},Int64})

    return @. loggamma(N+1) - loggamma(k+1) - loggamma(N-k+1)
end

"""
Coefficients for projection from a different sfs size.

# Arguments
 - `proj_to::Int64`: Numper of samples to project down to.
 - `proj_from::Int64`: Numper of samples to project from.
 - `hits::Int64`: Number of derived alleles projecting from.
# Output
 - `Vector{Float64}`: projection weights.
"""
function _cached_projection(proj_to::Int64, proj_from::Int64, hits::Int64)

    proj_hits = collect(0:proj_to)

    # For large sample sizes, we need to do the calculation in logs, and it is accurate enough for small sizes as well.
    lncontrib = _lncomb(proj_to,proj_hits)
    lncontrib += _lncomb(proj_from-proj_to,hits.-proj_hits)
    lncontrib = lncontrib .- _lncomb(proj_from, hits)
    contrib = @. exp(lncontrib)

    return contrib
end

"""
Project SFS to smaller sample size following DADI and moments implementation

# Arguments
 - `sfs::Matrix{Float64}`: SFS parsed from parse_sfs
 - `n::Int64`: new sample sample size
# Output:
 - `Vector{Float64}`: projected SFS.
"""
function project(sfs::Matrix{Float64},n::Int64)
    @assert size(sfs,1) > n "n is larger than the SFS. Select another sample size to project"

    n *= 2;

    proj_from     = size(sfs,1) + 1
    p_fs          = zeros(n + 1,2)
    # Assuming SFS contains freqs and only segregating sites, so not 0 and 1 frequency mutations
    sfs_fixations = vcat([0 0],view(sfs,:,2:3),[0 0])

    for hits = 1:proj_from
        # Adjust the slice in the array we're projecting from.
        # These are the least and most possible hits we could have in the projected fs
        least, most = max(n - (proj_from - hits), 1), min(hits,n);
        proj_slice = to_slice = least:most+1;
        # The projection weights.
        proj = _cached_projection(n, proj_from, hits)
        @. p_fs[to_slice,1]  += sfs_fixations[hits+1,1] * proj[proj_slice]
        @. p_fs[to_slice,2]  += sfs_fixations[hits+1,2] * proj[proj_slice]

    end

    p_sfs = hcat(collect(1:n-1),trunc.(@view p_fs[2:end-1,:]))
    # p_sfs = hcat(collect(1:n-1),@view p_fs[2:end-1,:])
    return p_sfs
end

function project(sfs::Vector{Matrix{Float64}},n::Int64)
    return ThreadsX.mapi(x -> project(x,n),sfs)
end

function fold(sfs::Matrix{Float64})

    n, m = size(sfs); n+=1;
    u_sfs = vcat(zeros(m-1)',view(sfs,:,2:m),zeros(m-1)')

    mask = 0:n .> (n / 2)
    # Remove non-mutation values too
    mask[1] = true;

    reversed = mask .* u_sfs
    reversed[.!mask,:] .= 0

    f_sfs = view(u_sfs + reverse(reversed,dims=1),.!mask,:)


    return(hcat(view(sfs,1:size(f_sfs,1)),f_sfs))

end

function fold(sfs::Vector{Matrix{Float64}})
    return ThreadsX.mapi(fold,sfs)
end

"""
Get standard Kingman SFS.

# Arguments
 - n: sample size
# Output
 - Expected SFS
"""
function standard_kingman(n::Int64)
    return 1 ./ collect(1:n)
end

"""
Calculate site-wise theta using Watterson's estimator.
"""
function θ(sfs::Matrix,divergence::Matrix)
    n_polymorphic = sum(@view sfs[:,3])
    n_sites = n_polymorphic + divergence[4]
    return n_polymorphic / sum(1 ./ collect(1:size(sfs,1))) / n_sites
end


function rᵢ!(sfs::Matrix,divergence::Matrix)

    # sfs_adjusted = copy(sfs)

    # normalized counts of the standard coalescent
    counts_kingman = standard_kingman(size(sfs,1)) * θ(sfs,divergence) * (sum(sfs[:,3])+divergence[4])

    # These counts transform the standard Kingman case to the observed
    # neutral SFS when multiplied and thus account for demography as we assume
    # the distortion in counts_neut to be due to demography only.
    # r = counts_kingman ./ sfs_adjusted[:,3]
    r = counts_kingman ./ sfs[:,3]

    # adjust for demography and polarization error
    # sfs_adjusted[:,2] .= r .*  @view sfs[:,2]
    sfs[:,2] .= r .*  @view sfs[:,2]

    return sfs
    # return sfs_adjusted
end

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

function write_files(x::DataFrame, name::String, newfile::Bool)
    CSV.write(name, x, delim = '\t', header = false, append = newfile)
end;

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
        n *= 2
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

    output_dofe   = @. tempname() * idx
    output_grapes = @. output_dofe * "." * model

    # GRAPES' header: a free-text title line, then "#unfolded". Plain text: written as a
    # DataFrame through CSV, the empty cells came out as `""` on both lines.
    for f in output_dofe
        Base.open(f, "w") do io
            print(io, "MKtest dofe\n#unfolded\n")
        end
    end
    @. write_files(dofe, output_dofe, true)

    @info "Running Grapes"
    # GammaZero has no beneficial class, so -no_div_param means nothing there. The test used
    # to read `model == "GammaZero" || no_div_param`, which reset the flag whenever it was
    # set and made -no_div_param impossible to pass.
    if model == "GammaZero"
        no_div_param = false
    end

    # Vectors splat into separate Cmd arguments; an empty vector adds none. The options
    # used to be joined into one string, which a Cmd passes as ONE argument
    # ("-no_div_data -fold"), and an empty one as an empty argument.
    bool_options = ["-no_div_data", "-no_div_param", "-no_syn_orient_error", "-fold"]
    bool_options = bool_options[[no_div_data, no_div_param, no_syn_orient_error, fold]]

    fixed_param = isempty(fixed_param) ? String[] : ["-fixed_param", fixed_param]

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

                idx_1 = findall(occursin.(model,names(df)))
                idx_2 = findfirst(occursin.("theta",names(df)))
                df = hcat(df[:,1:20],df[:,idx_1],df[:,idx_2:end])
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

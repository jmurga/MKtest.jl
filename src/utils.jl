function cumulative_vector!(v::Union{SubArray,Vector{T}}) where {T<:Number}
    @inbounds for i = length(v)-1:-1:1
        v[i] += v[i+1]
    end
    return v
end

"""
    cumulative_sfs(sfs, freqs = true)

Tail sums of an SFS: row i holds the counts of rows `i, i+1, ..., end`, so `Pn(x)` and
`Ps(x)` become the numbers of polymorphisms at frequency `>= x`.

The original asymptotic MK test uses the counts at x itself; the tail sums have the same
asymptote and are less sensitive to the sample size. A row whose tail sum (over all
count columns) is not positive is set to 0.

# Arguments
 - `sfs::Union{Array,SubArray}`: an SFS matrix, `[DAC Pn Ps]` from [`parse_sfs`](@ref).
 - `freqs::Bool`: column 1 holds the DAC (or frequency) and is copied unchanged;
   `false` cumulates every column.

# Returns
 - A matrix of the same size and element type as `sfs`.
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
    reduce_sfs(sfs, bins; n = nothing)

Collapse an SFS into `bins` derived allele frequency bins of equal width.

Row `i` has frequency `x = DAC / n`; it goes to bin `k` when `(k - 1) / bins < x <= k / bins`,
and the counts of columns 2 to end are summed within each bin. Bins that receive no row
hold 0.

# Arguments
 - `sfs::AbstractMatrix`: an SFS `[DAC Pn Ps ...]` (any number of count columns).
 - `bins::Int`: number of frequency bins.

# Keywords
 - `n::Union{Nothing,Int}`: haploid sample size. By default `size(sfs, 1) + 1`, which needs
   the full SFS (DAC `1..n-1`). Pass it for any SFS trimmed by `parse_sfs(...; cutoff)`:
   a trim at the low end throws without it, but a top-only trim (DAC `1..k`) cannot be
   told apart from a full SFS of `k + 1` chromosomes and would give wrong frequencies.

# Returns
 - `Matrix{Float64}` of `bins` rows: column 1 is the upper edge `k / bins` of each bin (a
   frequency, which `fww_mk`, `imputed_mk` and `standard_mk` accept as such), the other
   columns the summed counts.

# Example
```julia
sfs = hcat(1:9, ones(9), 2 .* ones(9))   # n = 10 chromosomes
MKtest.reduce_sfs(sfs, 2)                # [0.5 5 10; 1.0 4 8]
```
"""
function reduce_sfs(sfs::AbstractMatrix, bins::Int; n::Union{Nothing,Int}=nothing)
    bins > 0 || throw(ArgumentError("bins must be positive, got $bins"))
    freq = view(sfs, :, 1) ./ _alpha_x_n(sfs, n)
    out = zeros(bins, size(sfs, 2))
    out[:, 1] .= (1:bins) ./ bins
    for i in axes(sfs, 1)
        # the bin of x: ceil(x bins), with x = k / bins exactly kept in bin k
        k = clamp(ceil(Int, round(freq[i] * bins; digits = 10)), 1, bins)
        @views out[k, 2:end] .+= sfs[i, 2:end]
    end
    return out
end

"""
    _lncomb(N, k)

Log of the binomial coefficient `N choose k`, through `loggamma`; broadcast over a
vector `k`.
"""
function _lncomb(N::Int64,k::Union{Vector{Int64},Int64})

    return @. loggamma(N+1) - loggamma(k+1) - loggamma(N-k+1)
end

"""
    _cached_projection(proj_to, proj_from, hits)

Hypergeometric weights of `0:proj_to` derived alleles in a subsample of `proj_to`
chromosomes drawn from `proj_from` chromosomes carrying `hits` derived alleles,
computed in logs. Nothing is cached, despite the name.
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
    project(sfs, n)

Project an SFS down to `n` diploid individuals (`2n` chromosomes), by hypergeometric
sampling as in dadi and moments.

The projected counts are truncated to integers. Only columns 2 and 3 are projected, and
the input is taken to be a full spectrum of `size(sfs, 1) + 1` chromosomes: an SFS
trimmed by `parse_sfs(...; cutoff)` gives a wrong projection. For data read with
`ploidy = 1`, `n` still means `2n` chromosomes.

# Arguments
 - `sfs::Matrix{Float64}` (or a `Vector` of them, projected in parallel threads): the
   full raw SFS `[DAC Pn Ps]` from [`parse_sfs`](@ref).
 - `n::Int64`: number of diploid individuals to project to. Throws an `AssertionError`
   unless `n < size(sfs, 1)`.

# Returns
 - `Matrix{Float64}` `[DAC Pn Ps]` with `2n - 1` rows (a `Vector` of them for a
   `Vector` input).
"""
function project(sfs::Matrix{Float64},n::Int64)
    @assert 2n <= size(sfs,1) + 1 "n is larger than the SFS. Select another sample size to project"

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

"""
    fold(sfs)

Fold an SFS into the minor allele spectrum: the counts at DAC `i` and `m - i` are added,
for `m = size(sfs, 1) + 1` chromosomes.

The input must be a full, untrimmed spectrum. All columns after the first are folded; the
middle bin (`i = m / 2`, `m` even) is counted once.

# Arguments
 - `sfs::Matrix{Float64}` (or a `Vector` of them, folded in parallel threads): raw SFS
   `[DAC Pn Ps]`.

# Returns
 - `Matrix{Float64}` with `floor(m / 2)` rows, `[DAC Pn Ps]` for DAC `1:floor(m / 2)`
   (a `Vector` of them for a `Vector` input).
"""
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
    standard_kingman(n)

Shape of the neutral SFS under the standard coalescent, `1 ./ (1:n)` (unscaled by θ).
"""
function standard_kingman(n::Int64)
    return 1 ./ collect(1:n)
end

"""
    θ(sfs, divergence)

Per-site Watterson θ from the synonymous polymorphisms: `Ps / a / (Ps + Ls)`, with
`a = Σ 1/i` over the `size(sfs, 1)` SFS rows and `Ls = divergence[4]`.

`divergence` must be a `Matrix`; the `Vector` returned by [`parse_sfs`](@ref) has no
method.
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

"""
    write_files(x::Union{Matrix{Float64},SubArray}, name)
    write_files(x::DataFrame, name, newfile)

Write `x` to `name` as tab-separated text without a header.

A matrix is written gzip-compressed when `name` ends in `.gz`, replacing any existing
file. A `DataFrame` is appended to `name` when `newfile = true` (despite the name), and
replaces it when `false`; it is never compressed.
"""
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

"""
    read_gz(file)
    read_gz(files::Vector{String})

Read a tab-separated numeric file, plain or gzipped, into a `Matrix{Float64}` (no
header). A file that cannot be read gives a warning and `zeros(1, 12)`. The `Vector`
form reads the files in parallel threads and returns a `Vector` of matrices.
"""
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
    grapes(sfs, divergence, model; n = 0, nearly_neutral = 5, FWW_threshold = 0.15,
           nb_rand_start = 0, anc_to_rec_Ne_ratio = 1.0, no_div_data = false,
           no_div_param = false, no_syn_orient_error = false, fold = false,
           fixed_param = "")

Run the GRAPES binary (Galtier 2016) on each dataset and return its estimates.

The binary comes from the `grapes-static` package of the `genomedk` conda channel,
installed through CondaPkg on first use. Each dataset is written to a temporary `.dofe`
file and the runs go in parallel threads; the files of a successful run are removed.
For a Julia implementation of the same model see `dfe_ml`.

# Arguments
 - `sfs`: raw SFS `[DAC Pn Ps]`, a `Matrix` or a `Vector` of them. It must be a full
   spectrum (GRAPES reads `size(sfs, 1) + 1` chromosomes, or `2n` after projection).
 - `divergence`: `[Dn, Ds, Ln, Ls]`, a `Vector` or a `Vector` of them. Ln and Ls are
   required.
 - `model::String`: `"GammaZero"`, `"GammaExpo"`, `"GammaGamma"`, `"DisplGamma"`,
   `"ScaledBeta"`, `"FGMBesselK"` or `"all"`; anything else throws an `AssertionError`.

# Keywords (GRAPES options)
 - `n::Int64`: if positive, [`project`](@ref) the SFS to `n` diploid individuals first.
 - `nearly_neutral::Int64`: `-nearly_neutral`, the S threshold above which a mutation
   counts as adaptive.
 - `FWW_threshold::Float64`: `-FWW_threshold`, the lowest frequency kept in GRAPES' FWW
   alpha.
 - `nb_rand_start::Int64`: `-nb_rand_start`, random starting points of the optimizer.
 - `anc_to_rec_Ne_ratio::Float64`: `-anc_to_rec_Ne_ratio`, the ratio of the divergence to
   the polymorphism Ne.
 - `no_div_data::Bool`: `-no_div_data`, estimate the DFE from the polymorphism alone.
 - `no_div_param::Bool`: `-no_div_param`, the [-A] version of Galtier (2016); ignored
   for GammaZero.
 - `no_syn_orient_error::Bool`: `-no_syn_orient_error`, the same misorientation rate for
   synonymous and nonsynonymous sites.
 - `fold::Bool`: `-fold`, fit the folded SFS.
 - `fixed_param::String`: `-fixed_param`, the path of a GRAPES control file of
   parameters fixed to given values; empty for none.

A `0` anywhere in `divergence` (including Ln = Ls = 0 from a table without site counts)
throws an `ArgumentError`.

# Returns
 - `DataFrame`: the GRAPES output tables of the datasets that ran, stacked. For a single
   model the columns are GRAPES' first 20, those of `model`, and those from `theta` on;
   `"all"` keeps every column. A dataset that fails gives a warning and no row,
   so rows then no longer follow the dataset order.
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

    # GRAPES' header: a blank line, then "#unfolded". Plain text: written as a
    # DataFrame through CSV, the empty cells came out as `""` on both lines.
    for f in output_dofe
        Base.open(f, "w") do io
            print(io, "\n#unfolded\n")
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

"""
    sfs_to_dofe(pn, ps, dn, ds, ln, ls, w, n)

One GRAPES `.dofe` data line as a one-row `DataFrame`: the name `"dofe_<w>"`, the number
of chromosomes `n`, `Ln`, the Pn spectrum, `Ls`, the Ps spectrum, then `Ln, Dn, Ls, Ds`
(the same site counts for polymorphism and divergence). `pn` and `ps` are `1 x (n - 1)`
row vectors.
"""
function sfs_to_dofe(pn, ps, dn, ds, ln, ls, w, n)
    return(DataFrame(hcat("dofe_" * string(w), n, ln, pn, ls, ps, ln, dn, ls, ds...), :auto))
end

"""
    read_clean_grapes(dofe, grapes_file, model)

Read the GRAPES output `grapes_file` into a `DataFrame` and delete it, its `.profile`
and `.messages` files and the input `dofe`. For a single `model` the first 20 columns,
the columns of `model` and those from `theta` on are kept; `"all"` keeps the whole table.
If any step fails (no output, or a file missing) it warns and returns `nothing`, leaving
the files in place.
"""
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

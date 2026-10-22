
#################
# ABC functions #
#################

q_lower(x::Vector{Float64}) = quantile(x, [0.025])
q_upper(x::Vector{Float64}) = quantile(x, [0.975])


# A positive whole number given as any Real (1000, 1000.0, 1e3) as an Int.
function _as_count(x::Real, name::String)
    (isfinite(x) && isinteger(x) && x > 0) ||
        throw(ArgumentError("$name must be a positive whole number, got $x"))
    return Int(x)
end


# Index of the dataset a file belongs to: 10 for `summstat_10.txt.gz` or `out_10.0.post.gz`.
_dataset_index(file::String) = parse(Int, match(r"^[a-z]+_(\d+)", basename(file))[1])

# Files of `output_folder` starting with `prefix`, in dataset order. `readdir` sorts names
# as text (summstat_1, summstat_10, summstat_100, summstat_11, ...), which would return
# the posteriors in that order instead of the order of the input datasets.
function _dataset_files(output_folder::String, prefix::String)
    files = filter(x -> startswith(basename(x), prefix), readdir(output_folder, join = true))
    return sort(files, by = _dataset_index)
end

# Columns of a (possibly gzipped) table, from its first line only.
_count_cols(file::String) = _with_gz(io -> count(==('\t'), readline(io)) + 1, file)

# Posterior column names for P = 8 (no ω block) or 12 parameters.
_abc_param_names(P::Integer) = P == 12 ?
    [:α_weak, :α_strong, :α, :ωₐ_weak, :ωₐ_strong, :ωₐ, :ωₙₐ, :γ_deleterious, :β, :γ_weak, :γ_strong, :B] :
    [:α_weak, :α_strong, :α, :γ_deleterious, :β, :γ_weak, :γ_strong, :B]

# A dataset without site counts (Ln or Ls = 0) has no ω block in its summaries, so its
# ABC runs on 8 parameters. Its ω columns are then reported as NaN, after the α columns,
# so every posterior has them. Each parameter is regressed on its own over the same
# accepted rows, so the other columns are what the ABC gives on the table without the
# block. (A draw with Ds = 0 drops its row only; see `_assemble_summaries!`.)
_with_ω_columns(post::AbstractMatrix, P::Integer) =
    P == 12 ? post : hcat(post[:, 1:3], fill(NaN, size(post, 1), 4), post[:, 4:end])

function _report_ω_missing(Ps)
    k = count(<(12), Ps)
    k > 0 && @info "ω not estimated for $k of $(length(Ps)) datasets (Ln or Ls = 0: " *
                   "ω is undefined); their ω columns are NaN"
end

"""
    ABCreg(; output_folder, S, P = 12, tol, rm_summaries = false, abcreg = nothing,
           T = false, L = false) -> Vector{DataFrame}

Run the ABCreg binary (rejection plus local-linear regression ABC,
[Thornton 2009](https://doi.org/10.1186/1471-2156-10-35)) on every dataset of
`output_folder`, as written by [`summary_statistics`](@ref), and read the adjusted
posteriors.

Each pair `alphas_i.txt.gz` / `summstat_i.txt.gz` is run on its own thread, accepting the
`tol` closest table rows; ABCreg writes `out_i.0.post.gz`. A dataset whose summaries have
no ω columns (no site counts, Ln or Ls = 0) runs with `P - 4` and reports its ω columns as NaN.

# Keywords
 - `output_folder::String`: folder with the observed α(x) and the summaries.
 - `S::Int64`: number of summary statistics, `length(param.dac)`.
 - `P::Int64 = 12`: number of parameters (leading columns of the summaries).
 - `tol::Real`: number of accepted rows, a positive whole number (`1000`, `1e3`).
 - `rm_summaries::Bool = false`: delete the `summstat_*`, `alphas_*` and `out_*` files
   afterwards.
 - `abcreg::Union{Nothing,String} = nothing`: path to the binary; `nothing` uses the
   copy shipped in `scripts/reg`.
 - `T`, `L::Bool = false`: not passed to ABCreg; only checked against each other.

# Returns
One `DataFrame` per dataset, in dataset order, with columns
`α_weak, α_strong, α, ωₐ_weak, ωₐ_strong, ωₐ, ωₙₐ, γ_deleterious, β, γ_weak, γ_strong, B` (12); the
four ω columns are NaN for a dataset whose summaries have none. `γ_deleterious` is the mean
deleterious coefficient (absolute value), `β` the DFE shape, `γ_weak`/`γ_strong` the weakly and
strongly beneficial coefficients, and `B` is B x 1000.
Datasets whose ABCreg output is missing or empty are dropped with a warning, so the
rows after them no longer match the dataset numbers.
"""
function ABCreg(;
    output_folder::String,
    S::Int64,
    P::Int64 = 12,
    tol::Real,
    rm_summaries::Bool = false,
    abcreg::Union{Nothing,String}=nothing,
    T::Bool=false,
    L::Bool=false
)

    @assert !(T && L) "Cannot perform both transformation at the same time"

    if isnothing(abcreg)
        abcreg = abspath(joinpath(@__DIR__, "..", "..", "scripts", "reg"))
    end
     
    # List alphas and summstat files
    a_file   = _dataset_files(output_folder, "alphas_")
    sum_file = _dataset_files(output_folder, "summstat_")

    tol = _as_count(tol, "tol")

    # Per file: P, or P - 4 when that dataset's summaries have no ω block (see
    # `_with_ω_columns`). Taking it from the first file alone made the binary read a mixed
    # set's summary columns as parameters.
    Ps = ThreadsX.map(f -> _count_cols(f) == P + S ? P : P - 4, sum_file)

    # Creating output names
    out = output_folder .* "/out_" .* string.(_dataset_index.(a_file))
    # One task per file: counting 1000 files serially delays the binary.
    tol_abcreg = tol ./ ThreadsX.map(f -> _with_gz(countlines, f), sum_file)
    @info "Running ABCreg. Please cite https://doi.org/10.1186/1471-2156-10-35"
    # Using mapi instead map to limit tasks. Bash interaction not working as expected
    ThreadsX.mapi(
        (x, y, z, t, Pi) -> run(`$abcreg -d $x -p $y -P $Pi -S $S -t $t -b $z`),
        a_file,
        sum_file,
        out,
        tol_abcreg,
        Ps,
        ntasks = Threads.nthreads(),
    )

    @info "Opening and filtering posteriors distributions"
    # Exactly the file each dataset's run wrote, so each posterior stays paired with its P
    # (a glob over `out_*post*` lost the pairing when a run wrote nothing).
    out = out .* ".0.post.gz"

    # Move from CSV.read to readdlm, bug when threads enable in node server.
    # open(x) = Array(CSV.read(x, DataFrame))
    # Control outlier inference. 2Nes non negative values
    posteriors = read_gz(out)
    # Remove if some file was not computed. The rows after it then no longer match the
    # dataset numbers.
    ok = @. !iszero(posteriors)
    all(ok) || @warn "ABCreg failed for datasets $(_dataset_index.(out[.!ok])); their posteriors are dropped"
    posteriors = [_with_ω_columns(x, Pi) for (x, Pi) in zip(posteriors[ok], Ps[ok])]
    _report_ω_missing(Ps[ok])

    c_names = _abc_param_names(size(first(posteriors), 2))
    posteriors = map(x -> DataFrame(x,c_names),posteriors)

    # Remove summstat files
    if rm_summaries
        rm.(
            filter(
                x -> occursin("summstat", x) || occursin("alphas_", x) || occursin("out_", x),
                readdir(output_folder, join = true),
            )
        )
    end

    return posteriors
end

function get_mode(posterior::Matrix{Float64})

    # Allocating outputs
    out = zeros((1, size(posterior, 2)))

    for j::Int64 = 1:size(posterior, 2)
        y = kde(posterior[:, j])
        m = collect(y.x)[y.density.==maximum(y.density)][1]
        out[j] = m
    end

    return (out)
end

function get_mode(posterior::DataFrame)

    # Allocating outputs
    out = zeros((1, size(posterior, 2)))

    for j::Int64 = 1:size(posterior, 2)
        y = kde(posterior[:, j])
        m = collect(y.x)[y.density.==maximum(y.density)][1]
        out[j] = m
    end

    return (out)
end



"""
    AbstractPointEstimate

How [`summary_abc`](@ref) (and so [`abcmk`](@ref)) reduces each posterior column to the
point estimate in `:inference`: [`PosteriorMedian`](@ref) (default), [`PosteriorMean`](@ref),
[`PosteriorMode`](@ref) or [`PosteriorQuantile`](@ref). A new estimator is a subtype plus
one method, `MKtest.point_estimate(::MyEstimate, v::Vector{Float64}) -> Float64`.
"""
abstract type AbstractPointEstimate end

"""
    PosteriorMedian()

Posterior median of each column: the ABC-MK point estimate (default).
"""
struct PosteriorMedian <: AbstractPointEstimate end

"""
    PosteriorMean()

Posterior mean of each column. Linear, so it keeps `α = α_weak + α_strong` (to the 5-digit
rounding of the summaries) without [`Composed`](@ref).
"""
struct PosteriorMean <: AbstractPointEstimate end

"""
    PosteriorMode()

Posterior mode of each column, as the R package abc: the maximum of a kernel density
estimate (KernelDensity.jl defaults). The ABC-MK v1.0.0 point estimate.
"""
struct PosteriorMode <: AbstractPointEstimate end

"""
    PosteriorQuantile(p)

Posterior `p` quantile of each column, `0 <= p <= 1` (Statistics' default definition, the
one the `:stats` percentile rows use): `PosteriorQuantile(0)` is the minimum,
`PosteriorQuantile(1)` the maximum, `PosteriorQuantile(0.05)` a one-sided 95% lower bound.
Under [`Composed`](@ref) the derived strong component, e.g.
`q_0.975(α) - q_0.975(α_weak)`, is not a quantile of the α_strong posterior and carries no
coverage probability.
"""
struct PosteriorQuantile <: AbstractPointEstimate
    p::Float64
    function PosteriorQuantile(p::Real)
        0 <= p <= 1 || throw(ArgumentError("PosteriorQuantile: p = $p is outside [0, 1]"))
        return new(Float64(p))
    end
end

"""
    point_estimate(est::AbstractPointEstimate, v::Vector{Float64}) -> Float64

The point estimate `est` of one posterior column `v`.
"""
function point_estimate end
point_estimate(::PosteriorMedian, v::Vector{Float64}) = median(v)
point_estimate(::PosteriorMean, v::Vector{Float64}) = mean(v)
point_estimate(::PosteriorMode, v::Vector{Float64}) = _kde_mode(v)
point_estimate(e::PosteriorQuantile, v::Vector{Float64}) = quantile(v, e.p)

# Maximum of a kernel density estimate, the first grid point if several tie (as R abc).
function _kde_mode(v::Vector{Float64})
    y = kde(v)
    return collect(y.x)[y.density.==maximum(y.density)][1]
end

"""
    AbstractComponents

How [`summary_abc`](@ref) (and so [`abcmk`](@ref)) reports the weak and strong components
of α and ωₐ in the point estimate: [`Marginal`](@ref) (default) or [`Composed`](@ref).
"""
abstract type AbstractComponents end

"""
    Marginal()

Each column summarised on its own (default). Every posterior draw satisfies
`α = α_weak + α_strong`, but the median, the mode and a quantile of the three marginals
need not add up: on 825-gene sets the median gap is 0.01-0.02 and the mode gap 0.03-0.05.
"""
struct Marginal <: AbstractComponents end

"""
    Composed()

The strong components of the point estimate are derived from the other two,
`α_strong = α - α_weak` and `ωₐ_strong = ωₐ - ωₐ_weak` (ωₐ only when the posterior has the
ω block), so `α == α_weak + α_strong` holds exactly in `:inference`. The 95% interval in
`:quantiles` and every row of `:stats` stay marginal.
"""
struct Composed <: AbstractComponents end

"""
    compose!(components::AbstractComponents, point::Vector{Float64}, names) -> point

Apply the components rule to one dataset's point estimates, `point[j]` for column
`names[j]`.
"""
function compose! end
compose!(::Marginal, point::Vector{Float64}, names) = point
function compose!(::Composed, point::Vector{Float64}, names)
    for (w, s, t) in ((:α_weak, :α_strong, :α), (:ωₐ_weak, :ωₐ_strong, :ωₐ))
        jw, js, jt = findfirst(==(w), names), findfirst(==(s), names), findfirst(==(t), names)
        any(isnothing, (jw, js, jt)) && continue
        point[js] = point[jt] - point[jw]
    end
    return point
end

# NaN for a column with NaN: an ω the ABC did not estimate (see `_with_ω_columns`).
_point(est::AbstractPointEstimate, c::AbstractVector) =
    any(isnan, c) ? NaN : point_estimate(est, Vector{Float64}(c))

# The seven statistics of one posterior column, in `summary_abc`'s order; NaN as `_point`.
function _posterior_stats(c::AbstractVector)
    any(isnan, c) && return fill(NaN, 7)
    v = Vector{Float64}(c)
    return [minimum(v), only(q_lower(v)), median(v), mean(v), _kde_mode(v), only(q_upper(v)),
            maximum(v)]
end

"""
    summary_abc(posteriors; stat = PosteriorMedian(), components = Marginal()) -> OrderedDict

Summarise each posterior: minimum, 2.5% quantile, median, mean, mode, 97.5% quantile and
maximum of every column, and one point estimate per column and dataset. The mode is the
maximum of a kernel density estimate, as in the R package abc.

# Arguments
 - `posteriors`: `Vector{DataFrame}` (as [`ABCreg`](@ref) returns) or
   `Vector{Matrix{Float64}}` with 8 or 12 columns, one posterior per dataset.

# Keywords
 - `stat::AbstractPointEstimate = PosteriorMedian()`: the point estimate,
   [`PosteriorMedian`](@ref)`()` (the ABC-MK point estimate), [`PosteriorMean`](@ref)`()`,
   [`PosteriorMode`](@ref)`()` (the v1.0.0 point estimate) or
   [`PosteriorQuantile`](@ref)`(p)`. A string is a `TypeError`.
 - `components::AbstractComponents = Marginal()`: [`Marginal`](@ref)`()` summarises every
   column on its own; [`Composed`](@ref)`()` reports `α_strong = α - α_weak` and
   `ωₐ_strong = ωₐ - ωₐ_weak` of the point estimates, so they compose exactly.

# Returns
An `OrderedDict` with
 - `:inference`: `DataFrame`, one row per dataset with the point estimate of every
   parameter;
 - `:quantiles`: the point estimates as strings with the marginal 95% interval,
   `"0.2 [0.1-0.3]"`, rounded to 3 digits;
 - `:stats`: `DataFrame` with a `Stats` column and, per dataset, 7 rows (one per marginal
   statistic). Neither `stat` nor `components` changes it.
"""
summary_abc(posteriors::Vector{DataFrame}; stat::AbstractPointEstimate = PosteriorMedian(),
            components::AbstractComponents = Marginal()) =
    _summary_abc(posteriors, stat, components)

function summary_abc(posteriors::Vector{Matrix{Float64}}; kwargs...)
    c_names(x) = size(x, 2) == 12 ?
        [:α_weak, :α_strong, :α, :ωₐ_weak, :ωₐ_strong, :ωₐ, :ωₙₐ, :γ_deleterious, :β, :γ_weak, :γ_strong, :B] :
        [:α_weak, :α_strong, :α, :γ_deleterious, :β, :γ_weak, :γ_strong, :B]
    return summary_abc([DataFrame(x, c_names(x)) for x in posteriors]; kwargs...)
end

function _summary_abc(posteriors::Vector{DataFrame}, est::AbstractPointEstimate,
                      components::AbstractComponents)
    @info "Computing statistics over posteriors distributions"

    labels = ["Min.", "2.5% Perc.", "Median", "Mean", "Mode", "97.5% Perc.", "Max"]
    c_names = propertynames(first(posteriors))

    # One 7 x P block per dataset, rows in the order of `labels`, and one point estimate
    # per column with the components rule applied.
    blocks = ThreadsX.map(x -> reduce(hcat, _posterior_stats.(eachcol(x))), posteriors)
    points = ThreadsX.map(x -> compose!(components, [_point(est, c) for c in eachcol(x)], c_names),
                          posteriors)

    df = vcat([hcat(DataFrame(:Stats => labels), DataFrame(b, c_names)) for b in blocks]...)

    inference = DataFrame(reduce(vcat, (permutedims(p) for p in points)), c_names)
    quantiles = DataFrame(
        reduce(vcat, (permutedims(string.(round.(p, digits = 3), " [",
                                          round.(b[2, :], digits = 3), "-",
                                          round.(b[6, :], digits = 3), "]"))
                      for (p, b) in zip(points, blocks))),
        c_names,
    )

    length(posteriors) == 1 && @info df

    return OrderedDict(:inference => inference, :quantiles => quantiles, :stats => df)
end

#############


"""
    abc(; output_folder, S, P = 12, tol, transformation = "none",
        kernel = "epanechnikov", rm_summaries = false) -> OrderedDict

Pure-Julia ABC on every dataset of `output_folder`, as written by
[`summary_statistics`](@ref): rejection followed by local-linear regression adjustment
with heteroscedastic correction, a port of the R package abc (`method = "loclinear"`).
The alternative to [`ABCreg`](@ref) in `abcmk(...; abc_method = "abc")`.

Each pair `alphas_i.txt.gz` / `summstat_i.txt.gz` goes through `abc_loclinear`,
one dataset per thread (BLAS runs on one thread meanwhile). A dataset whose summaries have
no ω columns (no site counts, Ln or Ls = 0) runs with `P - 4` and reports its ω columns as NaN.

# Keywords
 - `output_folder::String`: folder with the observed α(x) and the summaries.
 - `S::Int64`: number of summary statistics, `length(param.dac)`.
 - `P::Int64 = 12`: number of parameters (leading columns of the summaries).
 - `tol::Real`: number of accepted rows, a positive whole number.
 - `transformation::String = "none"`: parameter transformation before the regression,
   `"none"`, `"log"` or `"tan"` (see `abc_loclinear`).
 - `kernel::String = "epanechnikov"`: weighting kernel, `"epanechnikov"`, `"gaussian"`,
   `"rectangular"`, `"triangular"`, `"biweight"` or `"cosine"`.
 - `rm_summaries::Bool = false`: delete the `summstat_*` and `alphas_*` files afterwards.

# Returns
`OrderedDict("rejection" => posteriors, "loclinear" => adjusted)`: two
`Vector{DataFrame}`, one per dataset, with the accepted and the regression-adjusted
values. Columns
`α_weak, α_strong, α, ωₐ_weak, ωₐ_strong, ωₐ, ωₙₐ, γ_deleterious, β, γ_weak, γ_strong, B` (12); the
four ω columns are NaN for a dataset whose summaries have none. `γ_deleterious` is the mean
deleterious coefficient (absolute value), `β` the DFE shape, `γ_weak`/`γ_strong` the weakly and
strongly beneficial coefficients, and `B` is B x 1000.
"""
function abc(; output_folder::String,
    S::Int64,
    P::Int64 = 12,
    tol::Real,
    transformation::String = "none",
    kernel::String = "epanechnikov",
    rm_summaries::Bool = false
)

    kernel, transformation = lowercase(kernel), lowercase(transformation)
    @assert lowercase(kernel) ∈ ["gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine"] "Kernel is incorrectly defined. Use gaussian, epanechnikov, rectangular, triangular, biweight or cosine"
    @assert lowercase(transformation) ∈ ["none", "log", "tan"] "Apply one of the following transformations: none, log, tan"

    # List/open alphas and summstat files
    targets = _dataset_files(output_folder, "alphas_")
    param_summaries = _dataset_files(output_folder, "summstat_")

    tol = _as_count(tol, "tol")

    # Change BLAS threads to 1 to proper parallel; restored even if a fit throws.
    nthreads_og = get_num_threads()
    set_num_threads(1)
    try

    # Per file: P, or P - 4 when that dataset's summaries have no ω block (see ABCreg)
    Ps = ThreadsX.map(f -> _count_cols(f) == P + S ? P : P - 4, param_summaries)

    @info "Running ABC"
    # Using mapi instead of map to limit tasks. Bash interaction not working as expected
    posteriors, posteriors_adjusted = unzip(ThreadsX.map((x, y, Pi) -> abc_loclinear(x,y, P=Pi, tol=tol, transformation=transformation, kernel=kernel), targets, param_summaries, Ps))

    @info "Filtering posterior distributions"
    # Remove if some file wasn't computed; a dataset without the ω block reports NaN ω
    posteriors = [_with_ω_columns(x, Pi) for (x, Pi) in zip(posteriors, Ps) if !iszero(x)]
    ok = .!iszero.(posteriors_adjusted)
    posteriors_adjusted = [_with_ω_columns(x, Pi) for (x, Pi) in zip(posteriors_adjusted[ok], Ps[ok])]
    _report_ω_missing(Ps[ok])

    # Remove summstat files
    if rm_summaries
        rm.(filter(x -> occursin("summstat", x) || occursin("alphas_", x), readdir(output_folder, join = true)))
    end

    c_names = _abc_param_names(size(first(posteriors_adjusted), 2))
    posteriors = map(x -> DataFrame(x,c_names),posteriors)
    posteriors_adjusted = map(x -> DataFrame(x,c_names),posteriors_adjusted)

    return OrderedDict("rejection"=>posteriors, "loclinear"=>posteriors_adjusted)
    finally
        set_num_threads(nthreads_og)
    end
end

function normalise(x::Union{SubArray,Float64},y::SubArray)
    if mad(y) != 0
        return x/mad(y)
    else
        return x
    end
end

function rejection(target::Vector{Float64},sumstat::Matrix{Float64},param::Matrix{Float64},tol::Float64,kernel::String)

    # Scale everything and euclidean distance
    sumstat_scaled = similar(sumstat)
    target_scaled  = similar(target)
    dist           = zeros(size(sumstat,1))

    for j in eachindex(target)
        sumstat_scaled[:,j] = normalise(view(sumstat,:,j),view(sumstat,:,j))
        target_scaled[j]    = normalise(target[j],view(sumstat,:,j))
        dist .+= (view(sumstat_scaled,:,j) .- target_scaled[j]).^2
    end
    dist = @. sqrt(dist)

    # Sort and get minimum distance to return values inside tolerance range
    n_accept = ceil(Int,length(dist)*tol)
    n_limit  = sort(dist)[n_accept]
    # Ensure get only n_limit, if more than one
    n_idx = findall(dist .<= n_limit)[1:n_accept]

    ## Weigthened distances
    wts = if(kernel == "epanechnikov")
        @.  1 - (dist[n_idx]/n_limit)^2
    elseif(kernel == "rectangular")
        @.  dist[n_idx]/n_limit
    elseif(kernel == "gaussian")
        @.  1/sqrt(2*pi)*exp(-0.5*(dist[n_idx]/(n_limit/2))^2)
    elseif(kernel == "triangular")
        @.  1 - abs(dist[n_idx]/n_limit)
    elseif(kernel == "biweight")
        @. (1 - (dist[n_idx]/n_limit)^2)^2
    elseif(kernel == "cosine")
        @. cos(pi/2*dist[n_idx]/n_limit)
    end


    return(target_scaled,sumstat[n_idx,:],sumstat_scaled[n_idx,:],param[n_idx,:],dist[n_idx],wts)
end

function regression(sumstat::Matrix{Float64},param::Matrix{Float64},weights::Vector{Float64})

    # Allocating results
    lm_residuals = zeros(size(sumstat,1),size(param,2))
    lm_coefficients = zeros(size(param,2),size(sumstat,2))

    # Loclinear regression
    for (i,v) in enumerate(eachcol(param))
        model = lm(sumstat,v,wts=weights)
        lm_coefficients[i,:] .= vec(coef(model))
        lm_residuals[:,i] .= vec(residuals(model))
    end

    return(lm_coefficients,lm_residuals)
end

function tangent_transfromation(x::Float64, minval::Float64, maxval::Float64)
    # Apply the tangent transformation of Hamilton et al. 2005 PNAS 7476
    return -log(1.0 / tan(((x - (minval - 1e-04)) / ((maxval + 1e-04) - (minval - 1e-04))) * (π / 2)))
end

function undo_tangent_transfromation(x::Float64, minval::Float64, maxval::Float64)
    # Undo the Hamilton tangent transformation
    return (minval - 1e-4) + (2 / π) * ((maxval + 1e-4) - (minval - 1e-4)) * atan(exp(x))
end

"""
    abc_loclinear(target_file, param_summaries_file; P, tol, transformation = "none",
                  kernel = "epanechnikov") -> (accepted, adjusted)

Rejection and local-linear regression ABC for one dataset, a port of `abc(...,
method = "loclinear")` from the R package abc.

`target_file` holds the observed summary statistics (one row) and
`param_summaries_file` the table, `P` parameter columns followed by the summary
statistics, as written by [`summary_statistics`](@ref).

1. Each statistic is scaled by its median absolute deviation over the table, and the
   `tol` rows closest to the target in Euclidean distance are accepted.
2. Accepted rows are weighted by `kernel` with bandwidth the largest accepted distance:
   `"epanechnikov"` (default), `"gaussian"`, `"rectangular"`, `"triangular"`,
   `"biweight"` or `"cosine"` (lower case).
3. The parameters are transformed: `"none"`, `"log"` (zeros replaced by the smallest
   non-zero value; throws an `ArgumentError` if any accepted value is negative, e.g. a
   negative `ωₙₐ`) or `"tan"` (the tangent transformation of Hamilton et al.
   2005, bounded by the table's range).
4. A weighted linear regression of each parameter on the statistics corrects the accepted
   values toward the target; a second regression of the log squared residuals rescales
   them (heteroscedastic correction).
5. Both sets are back-transformed.

# Keywords
 - `P::Int64`: number of parameter columns.
 - `tol::Int64`: number of accepted rows.
 - `transformation::String = "none"`, `kernel::String = "epanechnikov"`: as above.

# Returns
`(accepted, adjusted)`: two `tol x P` matrices, the rejection sample and the
regression-adjusted posterior.
"""
function abc_loclinear(target_file::String,param_summaries_file::String;P::Int64,tol::Int64,transformation::String="none",kernel::String= "epanechnikov")

    # Reading into matrix using Tables.matrix
    target = vec(CSV.read(target_file,matrix,ntasks=1,header=false))
    summaries = CSV.read(param_summaries_file,matrix,ntasks=1,header=false)
    # target = vec(readdlm(target_file,'\t',Float64,header=false))
    # summaries = readdlm(param_summaries_file,'\t',Float64,header=false))

    return abc_loclinear(target, summaries[:, 1:P], summaries[:, (P+1):end];
                         tol = tol, transformation = transformation, kernel = kernel)
end

# The same on matrices: `target` the observed summaries, `param` and `sumstat` the table's
# parameter and summary columns. The in-memory route (`abcmk(...; in_memory = true)`).
function abc_loclinear(target::AbstractVector, param::AbstractMatrix, sumstat::AbstractMatrix;
                       tol::Int64, transformation::String="none", kernel::String="epanechnikov")

    target  = Vector{Float64}(vec(target))
    param   = Matrix{Float64}(param)
    sumstat = Matrix{Float64}(sumstat)

    kernel, transformation = lowercase(kernel), lowercase(transformation)
    @assert lowercase(kernel) ∈ ["gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine"] "Kernel is incorrectly defined. Use gaussian, epanechnikov, rectangular, triangular, biweight or cosine"
    @assert lowercase(transformation) ∈ ["none", "log","tan"] "Apply one of the following transformations: none, log, tan"

    @assert length(target) == size(sumstat,2) "Number of summary statistics in target has to be the same as in sumstat."

    num_params = size(param,2)
    num_stats = size(sumstat,2)

    mins = minimum(param,dims=1)
    maxs = maximum(param,dims=1)

    ## Scale and euclidean distance
    ## #########################
    target_scaled, sumstat_accepted, sumstat_scaled, param_accepted, dist_accepted, wts = rejection(vec(target),sumstat,param,tol/size(sumstat,1),kernel)
    num_accepted = size(sumstat_accepted,1)

    ## Transform parameters
    ## ######################
    if transformation == "log" && any(<(0), param_accepted)
        throw(ArgumentError("log transformation needs non-negative parameters; use \"tan\" or \"none\""))
    end

    for (i,v) ∈ enumerate(eachcol(param_accepted))
        if transformation == "none"
            continue
        elseif transformation == "log"
            if(minimum(v) <= 0)
                v[findall(iszero, v)] .= minimum(filter(!iszero,v))
            end
            @. v .= log(v)
        elseif transformation == "tan"
            @. v = tangent_transfromation(v,mins[i],maxs[i])
        end
    end

    sumstat_intercept = hcat(ones(size(sumstat_scaled,1)),sumstat_scaled)

    # Linear regression
    lm_coefficients, lm_residuals = regression(sumstat_intercept,param_accepted,wts)


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
    lm_coefficients, lm_residuals = regression(sumstat_intercept,rsdl_log,wts)

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
            param_accepted[:,i] .= exp.(param_accepted[:,i])
            param_adjusted[:,i] .= exp.(param_adjusted[:,i])
        elseif(transformation == "tan")
            param_accepted[:,i] .= undo_tangent_transfromation.(param_accepted[:,i],mins[i],maxs[i])
            param_adjusted[:,i] .= undo_tangent_transfromation.(param_adjusted[:,i],mins[i],maxs[i])
        end
    end

    return param_accepted,param_adjusted
end

##################################
# ABCreg, in Julia and in memory #
##################################

# ABCreg's parameter transformations (transformations.cc): `-L` is log(x + 1e-4), `-T` the
# tangent transformation, both over the range of the whole table. Where the binary aborts
# (`-L` on a value <= -1e-4) the port gives NaN; each parameter is regressed on its own, so
# only that parameter's column of the posterior is NaN.
_abcreg_transform(x, t, lo, hi) =
    t == "log" ? (x + 1e-4 > 0 ? log(x + 1e-4) : NaN) : t == "tan" ? tangent_transfromation(x, lo, hi) : x
_abcreg_untransform(y, t, lo, hi) =
    t == "log" ? exp(y) - 1e-4 : t == "tan" ? undo_tangent_transfromation(y, lo, hi) : y

# gsl_stats_quantile_from_sorted_data
function _gsl_quantile(sorted::AbstractVector{Float64}, f::Float64)
    n = length(sorted)
    idx = f * (n - 1)
    lhs = floor(Int, idx)
    δ = idx - lhs
    lhs == n - 1 && return sorted[n]
    return (1 - δ) * sorted[lhs+1] + δ * sorted[lhs+2]
end

# The same on unsorted data: only the two order statistics it interpolates are needed.
function _gsl_quantile_unsorted(x::AbstractVector{Float64}, f::Float64)
    n = length(x)
    idx = f * (n - 1)
    lhs = floor(Int, idx)
    δ = idx - lhs
    v = copy(x)
    lhs == n - 1 && return maximum(v)
    a, b = partialsort!(v, (lhs+1):(lhs+2))
    return (1 - δ) * a + δ * b
end

# gsl_linalg_balance_columns: scale column j by the power of 2 that brings its L1 norm into
# [0.5, 1]. Returns the factors.
function _gsl_balance_columns!(A::AbstractMatrix{Float64})
    D = ones(size(A, 2))
    for j in axes(A, 2)
        s = sum(abs, view(A, :, j))
        f = 1.0
        (s == 0.0 || !isfinite(s)) && continue
        while s > 1.0
            s /= 2.0; f *= 2.0
        end
        while s < 0.5
            s *= 2.0; f /= 2.0
        end
        D[j] = f
        f != 1.0 && (view(A, :, j) ./= f)
    end
    return D
end

# gsl_multifit_wlinear_svd (GSL 2.x: balanced columns, singular values <= tol * s_max
# dropped), solved for every column of Y at once. Returns the (p x P) coefficients.
function _gsl_wlinear_svd(X::AbstractMatrix{Float64}, w::AbstractVector{Float64},
                          Y::AbstractMatrix{Float64}, tol::Float64)
    sw = sqrt.(max.(w, 0.0))
    A = X .* sw
    D = _gsl_balance_columns!(A)
    F = svd(A)
    inv_s = [s <= tol * F.S[1] ? 0.0 : 1 / s for s in F.S]
    return (F.V * (inv_s .* (F.U' * (Y .* sw)))) ./ D
end

"""
    abcreg_loclinear(target, param, sumstat; tol, transformation = "none") -> Matrix{Float64}

ABCreg ([Thornton 2009](https://doi.org/10.1186/1471-2156-10-35)) on one dataset, in
Julia and in memory: the algorithm of the binary that [`ABCreg`](@ref) runs, transcribed
from its source (`reg.cc`, `generate_posterior.cc`, `abc_algebra.cc`,
`transformations.cc`). Cite ABCreg when you use it.

`param` holds the table's parameter columns and `sumstat` its summary statistics, one row
per table row; `target` is the observed summary row.

1. Each statistic is centred and divided by its sample standard deviation over the table
   (the target with the same values).
2. The acceptance radius `q` is the `tol / nrow` quantile of the Euclidean distances
   (GSL's interpolated quantile), and rows with distance strictly below `q` are accepted.
   If `q` is 0, the exact matches are returned without regression, as ABCreg does.
3. Epanechnikov weights `1 - d^2 / q^2`.
4. Each parameter (transformed with the binary's `-L`, `log(x + 1e-4)`, or `-T`, the
   tangent transformation over the table's range, when `transformation` is `"log"` or
   `"tan"`) is regressed on the scaled statistics by weighted least squares through GSL's
   balanced SVD with tolerance `1 / n_accepted`.
5. Adjusted value = prediction at the target + the row's residual, back-transformed. No
   heteroscedastic correction (unlike [`abc_loclinear`](@ref)).

Under `"log"` a parameter at or below `-1e-4` (a negative `ωₙₐ`) makes the binary abort;
here that parameter's column is `NaN` and the other columns are unaffected.

Agrees with the binary to its printed precision (6 significant digits); the parity test is
`test/abc_in_memory.jl`.

# Returns
The adjusted posterior, one row per accepted table row in table order; zero rows if none
was accepted.
"""
function abcreg_loclinear(target::AbstractVector, param::AbstractMatrix, sumstat::AbstractMatrix;
                          tol::Integer, transformation::String = "none")
    transformation = lowercase(transformation)
    transformation ∈ ("none", "log", "tan") || throw(ArgumentError(
        "transformation must be \"none\", \"log\" or \"tan\", got \"$transformation\""))
    n, S = size(sumstat)
    length(target) == S || throw(DimensionMismatch(
        "target has $(length(target)) summaries, the table $S"))
    size(param, 1) == n || throw(DimensionMismatch("param and sumstat need the same rows"))
    1 <= tol <= n || throw(ArgumentError("tol must be in 1:$n, got $tol"))

    lo = vec(minimum(param, dims = 1))'
    hi = vec(maximum(param, dims = 1))'

    μ = vec(mean(sumstat, dims = 1))
    σ = vec(std(sumstat, dims = 1))
    z0 = (target .- μ) ./ σ
    # Distances without a scaled copy of the table: the same operations, in the same order,
    # one column at a time.
    dist = zeros(n)
    @inbounds for j in 1:S
        for i in 1:n
            dist[i] += ((sumstat[i, j] - μ[j]) / σ[j] - z0[j])^2
        end
    end
    dist .= sqrt.(dist)
    q = _gsl_quantile_unsorted(dist, tol / n)

    # Only the accepted rows are scaled and transformed (row by row, as before).
    if abs(q) <= eps(Float64)
        return _abcreg_untransform.(_abcreg_transform.(param[dist .<= q, :], transformation, lo, hi),
                                    transformation, lo, hi)
    end
    acc = findall(<(q), dist)
    isempty(acc) && return zeros(0, size(param, 2))

    w = @. 1.0 - dist[acc]^2 / q^2
    X = hcat(ones(length(acc)), (sumstat[acc, :] .- μ') ./ σ')
    Y = _abcreg_transform.(param[acc, :], transformation, lo, hi)
    C = _gsl_wlinear_svd(X, w, Y, 1 / length(acc))
    adjusted = (vcat(1.0, z0)' * C) .+ (Y .- X * C)
    return _abcreg_untransform.(adjusted, transformation, lo, hi)
end

"""
    _abc_in_memory(param, sfs, divergence; abc_method, tol, h5_file, rates, summstat_size, alpha,
                   B_bins, sampling, row_prior) -> Vector{DataFrame}

The ABC step of [`abcmk`](@ref) without files: the same table rows, Poisson realisations
and observed summaries as [`summary_statistics`](@ref) (the same random draws at the same
seed and thread count), handed per dataset to [`abcreg_loclinear`](@ref)
(`abc_method = "abcreg"`) or [`abc_loclinear`](@ref) (`"abc"`, its adjusted posterior).
One dataset per thread; only its posterior is kept. Datasets with an empty posterior are
dropped with a warning, as [`ABCreg`](@ref) does.
"""
function _abc_in_memory(param::parameters, sfs::Vector, divergence::Vector;
                        abc_method::String, tol::Int, h5_file::Union{Nothing,String} = nothing,
                        rates = nothing, summstat_size::Int,
                        alpha, B_bins, sampling::AbstractPoissonSampling,
                        row_prior::AbstractRowPrior)
    st = _summaries_setup(param, sfs, divergence; h5_file, rates, summstat_size, alpha, B_bins,
                          sampling, row_prior)
    S = length(param.dac)

    if abc_method == "abcreg"
        @info "Running ABCreg in memory (Julia port of ABCreg). Please cite https://doi.org/10.1186/1471-2156-10-35"
    else
        @info "Running ABC in memory"
    end

    # BLAS on one thread: the datasets run in parallel; restored even if a fit throws.
    # One buffer per running task, reused across datasets; the summaries are a view of it,
    # used before the buffer goes back to the pool.
    nthreads_og = get_num_threads()
    set_num_threads(1)
    pool = _BufferPool()
    posteriors = try
        ThreadsX.map(eachindex(st.sfs_p)) do i
            _with_buffer(pool, st) do buf
                out = _assemble_summaries!(buf, st.sfs_p[i], st.divergence_p[i],
                                           st.neut, st.sel, st.dsdn, st.gn, st.sh, st.gL, st.gH,
                                           st.B, st.summstat_size, st.samplers[i])
                # No ω block without site counts (Ln or Ls = 0): ABCreg / abc then read P = 8.
                P = size(out, 2) == 12 + S ? 12 : 8
                target = vec(st.α[i])
                pa, ss = view(out, :, 1:P), view(out, :, (P+1):size(out, 2))
                post = if abc_method == "abcreg"
                    abcreg_loclinear(target, pa, ss; tol)
                else
                    abc_loclinear(target, pa, ss; tol)[2]     # copies its inputs
                end
                (P, post)
            end
        end
    finally
        set_num_threads(nthreads_og)
    end

    ok = [!isempty(p) && !iszero(p) for (_, p) in posteriors]
    all(ok) || @warn "ABC gave no posterior for datasets $(findall(!, ok)); they are dropped"
    _report_ω_missing(first.(posteriors[ok]))
    # A dataset without the ω block reports NaN ω (see `_with_ω_columns`)
    full = [_with_ω_columns(p, P) for (P, p) in posteriors[ok]]
    c_names = _abc_param_names(size(first(full), 2))
    return [DataFrame(p, c_names) for p in full]
end

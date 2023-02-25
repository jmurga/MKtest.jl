
#################
# ABC functions #
#################

q_lower(x::Vector{Float64}) = quantile(x, [0.5])
q_upper(x::Vector{Float64}) = quantile(x, [0.95])
function write_files(x::Matrix{Float64}, name::String)
    CSV.write(name, Tables.table(x), delim = '\t', header = false)
end;

function write_files(x::DataFrame, name::String, newfile::Bool)
    CSV.write(name, x, delim = '\t', header = false, append = newfile)
end;

"""
	ABCreg(analysis_folder, S, P, tol, abcreg)

Performing ABC inference using ABCreg. Please, be sure the input folder contain the observed and summary statistics produced by MKtest.summary_statistics().

# Arguments
 - `analysis_folder::String` : folder containing the observed data and summary estatistics. It will be used to output the posterior distributions
 - `P::Int64` : number of parameters to infer.
 - `S::Int64` : number of summary stastitics to perform the inference.
 - `tol::Float64` : tolerance value. It defines the number of accepted values at ABC inference.
 - `abcreg::String` : path to ABCreg binary.
# Output
 - `posteriors:Vector{Matrix{Float64}}`: posteriors distributions estimated by ABCreg.
"""
function ABCreg(;
    analysis_folder::String,
    S::Int64,
    P::Int64 = 5,
    tol::Float64,
    rm_summaries::Bool = false,
)

    abcreg = abspath(joinpath(@__DIR__, "..", "scripts", "reg"))
    # List alphas and summstat files
    a_file = filter(x -> occursin("alphas", x), readdir(analysis_folder, join = true))
    sum_file = filter(x -> occursin("summstat", x), readdir(analysis_folder, join = true))

    # Creating output names
    out = analysis_folder .* "/out_" .* string.(1:size(a_file, 1))

    @info "Running ABCreg"
    # Using mapi instead map to limit tasks. Bash interaction not working as expected
    ThreadsX.mapi(
        (x, y, z) -> run(`$abcreg -d $x -p $y -P $P -S $S -t $tol -b $z`),
        a_file,
        sum_file,
        out,
        ntasks = Threads.nthreads(),
    )

    @info "Opening and filtering posteriors distributions"
    out = filter(x -> occursin("post", x), readdir(analysis_folder, join = true))
    out = filter(x -> !occursin(".1.", x), out)

    # Move from CSV.read to readdlm, bug when threads enable in node server.
    # open(x) = Array(CSV.read(x, DataFrame))
    # Control outlier inference. 2Nes non negative values
    posteriors = @. read_gz(out)
    # Remove is some file wasnt computed
    posteriors = posteriors[@. !iszero(posteriors)]

    # Remove summstat files
    if rm_summaries
        rm.(
            filter(
                x -> occursin("summstat", x) || occursin("alphas_", x),
                readdir(analysis_folder, join = true),
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

# read gz since CSV.read is bugged at node server
function read_gz(out_file::String)
    try
        x = readdlm(open(out_file), '\t', header = false)
        x = x[(x[:, 4].>0).&(x[:, 1].>0).&(x[:, 2].>0).&(x[:, 3].>0), :]
        return x
    catch
        @warn "$out_file is empty"
        return zeros(1, 5)
    end
end

"""
	summary_abc(posteriors,stat)
Posterior distributions statistics. The function estimates Min., 0.5% Perc., Median, Mean, Mode, "99.5% Perc., Max values. The mode is estimated following R package abc mode function. The argument `stat` define the output estimation.

# Arguments
 - `posterior::Matrix{Float64}` : posterior distribution.
 - `stat::String`: output estimation.
 - `plot_path{String,Nothing}`: path to save posterior plot.
# Output
 - `DataFrame`: posterior statistics.
 - `DataFrame`: chosen statistic inference.
"""
function summary_abc(posteriors::Vector{Matrix{Float64}}; stat::String = "Mean")
    @info "Computing statistics over posteriors distributions"

    p_min = vcat(map(x -> minimum(x, dims = 1), posteriors)...)
    p_max = vcat(map(x -> maximum(x, dims = 1), posteriors)...)
    p_mean = vcat(map(x -> mean(x, dims = 1), posteriors)...)
    p_median = vcat(map(x -> median(x, dims = 1), posteriors)...)
    p_mode = vcat(map(get_mode, posteriors)...)
    p_lower = vcat(map(x -> mapslices(q_lower, x, dims = 1), posteriors)...)
    p_upper = vcat(map(x -> mapslices(q_upper, x, dims = 1), posteriors)...)

    if length(posteriors) == 1
        stats = DataFrame(
            :Stats => repeat(
                ["Min.", "0.5% Perc.", "Median", "Mean", "Mode", "99.5% Perc.", "Max"],
                inner = length(posteriors),
            ),
        )

        tmp = DataFrame(
            vcat(p_min, p_lower, p_median, p_mean, p_mode, p_upper, p_max),
            [:α_weak, :α_strong, :α, :γ, :β],
        )

        df = hcat(stats, tmp)

        @info df

        return df, df[df.Stats.==uppercasefirst(stat), 2:end]
    else
        df = DataFrame[]
        for i::Int64 = 1:length(posteriors)
            stats = DataFrame(
                :Stats => [
                    "Min.",
                    "0.5% Perc.",
                    "Median",
                    "Mean",
                    "Mode",
                    "95% Perc.",
                    "Max",
                ],
            )

            tmp = DataFrame(
                hcat(
                    p_min[i, :],
                    p_lower[i, :],
                    p_median[i, :],
                    p_mean[i, :],
                    p_mode[i, :],
                    p_upper[i, :],
                    p_max[i, :],
                )',
                [:α_weak, :α_strong, :α, :γ, :β],
            )
            push!(df, hcat(stats, tmp))
        end

        df_vcat = vcat(df...)

        return df, df_vcat[df_vcat.Stats.==uppercasefirst(stat), 2:end]
    end
end

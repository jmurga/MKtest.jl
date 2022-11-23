function ABCreg(;
                analysis_folder::String,
                S::Int64,
                P::Int64 = 5,
                tol::Float64,
                abcreg::String,
                rm_summaries::Bool=false)

    # List alphas and summstat files
    a_file = filter(x -> occursin("alphas", x), readdir(analysis_folder, join = true))
    sum_file = filter(x -> occursin("summstat", x), readdir(analysis_folder, join = true))

    # Creating output names
    out = analysis_folder .* "/out_" .* string.(1:size(a_file, 1))

    @info "Running ABCreg"
    parallel_bin = CondaPkg.which("parallel")
    
    if isnothing(CondaPkg.which("parallel"))
        CondaPkg.add("parallel", channel = "conda-forge")
        parallel_bin = CondaPkg.which("parallel")
    end
    
    # Using GNU-parallel from CondaPkg
    n_jobs = Threads.nthreads();
    job_commands  = @. abcreg * " -d " * a_file * " -p " * sum_file* " -P " * string(P) * " -S " * string(S) * " -t " * string(tol) * " -b " * out;
    job_file = tempname(analysis_folder);

    CSV.write(job_file,Tables.table(job_commands),header=false);
    run(`$parallel_bin -j $n_jobs -u -a $job_file`)

    rm(job_file)

    @info "Opening and filtering posteriors distributions"
    out = filter(x -> occursin("post", x), readdir(analysis_folder, join = true))
    out = filter(x -> !occursin(".1.", x), out)

    # Control outlier inference. 2Nes non negative values
    # open(x) = Array(CSV.read(x, DataFrame))
    flt(x::Matrix{Float64}) = x[(x[:, 4] .> 0) .& (x[:, 1] .> 0) .& (x[:, 2] .> 0) .& (x[:, 3] .> 0), :]
    posteriors = map(x -> flt(readdlm(GZip.open(x))), out)
    # posteriors = flt.(open.(out))

    # Remove summstat files
    if rm_summaries
        rm.(filter(x -> occursin("summstat", x) || occursin("alphas_", x),
                   readdir(analysis_folder, join = true)))
    end

    return posteriors
end

function ABCreg(;
                analysis_folder::String,
                S::Int64,
                P::Int64 = 5,
                tol::Float64,
                abcreg::String,
                rm_summaries::Bool=false)

    # List alphas and summstat files
    a_file = filter(x -> occursin("alphas", x), readdir(analysis_folder, join = true))
    sum_file = filter(x -> occursin("summstat", x), readdir(analysis_folder, join = true))

    if(length(a_file)>100)
        a_file   = a_file[1:100];
        sum_file = sum_file[1:100];
    end

    # Creating output names
    out = analysis_folder .* "/out_" .* string.(1:size(a_file, 1))

    @info "Running ABCreg"
    # Using mapi to limit tasks. ThreadsX.map Not working properly with bash function. Try change to pipeline or something else.
    function r(a::String, s::String, o::String, abcreg::String = abcreg, P::Int64 = P,
               S::Int64 = S, tol::Float64 = tol)
        run(`$abcreg -d $a -p $s -P $P -S $S -t $tol -b $o`)
    end

    # Using mapi instead map. Bash interaction not working as expected
    # buffer_tasks = 1
    # if length(out) > 1
    #     buffer_tasks = Int(floor(length(out)/10))
    # end

    ThreadPools.tmap((x, y, z) -> r(x, y, z,abcreg,P,S,tol),
                  a_file,
                  sum_file,
                  out)

    @info "Opening and filtering posteriors distributions"
    out = filter(x -> occursin("post", x), readdir(analysis_folder, join = true))
    out = filter(x -> !occursin(".1.", x), out)

    # Control outlier inference. 2Nes non negative values
    open(x) = Array(CSV.read(x, DataFrame))
    flt(x) = x[(x[:, 4] .> 0) .& (x[:, 1] .> 0) .& (x[:, 2] .> 0) .& (x[:, 3] .> 0), :]
    posteriors = flt.(open.(out))

    # Remove summstat files
    if rm_summaries
        rm.(filter(x -> occursin("summstat", x) || occursin("alphas_", x),
                   readdir(analysis_folder, join = true)))
    end

    return posteriors
end

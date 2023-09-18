using Fire
"""
    julia --threads N abcmk_cli.jl rates 10000 661 --dac 1,2,4,5,10,20,50,100,200,400,500,661,925,1000 --gH 200,1000 --gL 1,10 --gam_flanking -1000,-500 --gam_dfe -1000,-100 --shape 0.184 --alpha 0.1,0.9 --iterations 1000000 --output rates.jld2

Function to solve analytical fixation rates and the expected SFS. The function will create N random models from prior values. Use the arguments to defined the input range for each parameter.

If rho and/or theta are not set, default values will be used (0.001).

To parallelize the estimations please be sure you start up Julia using --threads/-t option and set the number of cores.

The function returns a HDF5 file containing models solved and rates. The rates will be used to compute summary statistics required at ABC.

Please check the documentation to get more info about models parameters or detailed arguments description https://jmurga.github.io/MKtest.jl/dev/cli/ to check model
"""
@main function rates(N::Int64,
                     samples::Int64,
                     iterations::Int64;
                     alpha::String = "0.1,0.9",
                     gam_dfe::String = "-1000,-200",
                     gam_flanking::String = "-1000,-500",
                     gL::String = "5,10",
                     gH::String = "400,1000",
                     dac::String = "1,2,4,5,10,20,50,100,200,400,500,661,925,1000",
                     shape::Float64 = 0.184,
                     rho::Float64 = 0.001,
                     theta::Float64 = 0.001,
                     cutoff::String = "0.0,1.0",
                     output::String = "rates.jld2")

    @eval using MKtest

    alpha, cutoff = map(x -> parse.(Float64, x), split.([alpha, cutoff], ","))
    gam_dfe, gam_flanking, gL, gH, dac = map(x -> parse.(Int64, x),
                                             split.([gam_dfe, gam_flanking, gL, gH, dac],
                                                    ","))

    adap = MKtest.parameters(N = N, n = samples, dac = dac, cutoff = cutoff)

    MKtest.rates(adap,
                 gam_dfe = gam_dfe,
                 gam_flanking = gam_flanking,
                 gL = gL,
                 gH = gH,
                 alpha = alpha,
                 iterations = iterations,
                 output = output)
end

"""
    summary_stat

Estimate summary statistics using observed data and analytical rates.
"""
@main function summaries(N::Int64,
                         samples::Int64,
                         data::String;
                         genes::String = "<genes.txt>",
                         dac::String = "2,4,5,10,20,50,200,661,925",
                         cutoff::String = "0.0,1.0",
                         rates::String = "rates.jld2",
                         folder::String = "<folder>",
                         summsize::Int64 = 100000)
    @eval using MKtest

    cutoff, dac = map(x -> parse.(Float64, x), split.([cutoff, dac], ","))

    dac = Int64.(dac)

    adap = MKtest.parameters(N = N, n = samples, dac = dac, cutoff = cutoff)

    if genes == ""
        genes = nothing
    end

    @info "Parsing $data"
    alpha, sfs, divergence = MKtest.parse_sfs(adap, data = data, gene_list = genes)

    summstat = MKtest.summary_statistics(adap,
                                         sfs = sfs,
                                         divergence = divergence,
                                         h5_file = rates,
                                         analysis_folder = folder,
                                         summstat_size = summsize)
end

"""
    abc(nsumm,tol;folder,abcreg)

ABC inference and statistic. Please, be sure your analysis_folder contain the files produced by summaries. If you are using Docker or Singularity images you don't need to provide ABCreg path.

"""
@main function abc(nsumm::Int64,
                   tol::Float64;
                   folder::String = "<folder>",
                   abcreg::String = "reg")
    @eval using MKtest

    posteriors = MKtest.ABCreg(analysis_folder = folder,
                               P = 5,
                               S = nsumm,
                               tol = tol)

    out = MKtest.summary_abc(posteriors)

    @info "Posterior(s) distribution(s) was deposited at $folder"
end

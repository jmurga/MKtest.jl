using Fire   # Fire.@main: Julia >= 1.11 has its own Base.@main
using MKtest, CSV

# Command-line interface to ABC-MK. Two steps:
#
#   julia -t auto abcmk_cli.jl rates 1000 661 100000 --dac 2,4,5,10,20,50,200,661,925 --output rates_tgp.jld2
#   julia -t auto abcmk_cli.jl abcmk 1000 661 data/MKdata_tgp_may2023.txt rates_tgp.jld2 abc_out --dac 2,4,5,10,20,50,200,661,925
#
# `rates` builds the reference table once per population, sample size and cutoff;
# `abcmk` estimates alpha, alpha_W, alpha_S and the DFE of every dataset against it.
# The `dac` and `cutoff` of both calls must agree. See docs/src/abcmk.md.

_floats(x::String) = parse.(Float64, split(x, ","))
_ints(x::String) = parse.(Int64, split(x, ","))

"""
    rates N n iterations [--dac ...] [--alpha lo,hi] [--gam_dfe lo,hi] [--gL lo,hi]
          [--gH lo,hi] [--gam_flanking lo,hi] [--shape 0.184] [--cutoff 0.0,1.0]
          [--output rates.jld2]

Build the ABC-MK reference table with `MKtest.rates`: `iterations` models drawn from
the priors, for population size `N` and `n` diploid individuals. Start Julia with
`-t auto` to solve the models in parallel.
"""
Fire.@main function rates(N::Int64,
                     n::Int64,
                     iterations::Int64;
                     dac::String = "2,4,5,10,20,50,200,500,700",
                     alpha::String = "0.1,0.9",
                     gam_dfe::String = "-1000,-200",
                     gL::String = "1,10",
                     gH::String = "100,1000",
                     gam_flanking::String = "-1000,-500",
                     shape::Float64 = 0.184,
                     cutoff::String = "0.0,1.0",
                     output::String = "rates.jld2")
    param = MKtest.parameters(N = N, n = n, dac = _ints(dac), cutoff = _floats(cutoff),
                              shape = shape)
    MKtest.rates(param; iterations = iterations, output_file = output,
                 α = _floats(alpha), γ_deleterious = _ints(gam_dfe), γ_weak = _ints(gL),
                 γ_strong = _ints(gH), γ_flanking = _ints(gam_flanking))
    @info "Reference table written to $output"
end

"""
    abcmk N n data rates_file output_folder [--genes file] [--dac ...] [--cutoff 0.0,1.0]
          [--summsize 100000] [--tol 1000] [--method abcreg] [--stat median]

Run `MKtest.abcmk` on the datasets of `data` (one dataset, or one per line of `--genes`,
read with `MKtest.parse_sfs`) against the table `rates_file`. Datasets that fail
`MKtest.valid_data` are dropped. Writes `abcmk_inference.tsv` (one row per dataset) and
`abcmk_quantiles.tsv` to `output_folder`.
"""
Fire.@main function abcmk(N::Int64,
                     n::Int64,
                     data::String,
                     rates_file::String,
                     output_folder::String;
                     genes::String = "",
                     dac::String = "2,4,5,10,20,50,200,500,700",
                     cutoff::String = "0.0,1.0",
                     summsize::Int64 = 100000,
                     tol::Int64 = 1000,
                     method::String = "abcreg",
                     stat::String = "median")
    param = MKtest.parameters(N = N, n = n, dac = _ints(dac), cutoff = _floats(cutoff))
    sfs, div, _ = MKtest.parse_sfs(data; sample_size = n, cutoff = _floats(cutoff),
                                   gene_list = isempty(genes) ? nothing : genes)
    if sfs isa Matrix
        sfs, div = [sfs], [div]
    end
    ok = MKtest.valid_data(param, sfs, div)
    all(ok) || @warn "dropping $(count(!, ok)) of $(length(ok)) datasets (valid_data)"
    mkpath(output_folder)
    res = MKtest.abcmk(param, sfs[ok], div[ok]; h5_file = rates_file,
                       summstat_size = summsize, output_folder = output_folder,
                       tol = tol, abc_method = method, stat = stat)
    CSV.write(joinpath(output_folder, "abcmk_inference.tsv"), res[:inference]; delim = '\t')
    CSV.write(joinpath(output_folder, "abcmk_quantiles.tsv"), res[:quantiles]; delim = '\t')
    @info "Results written to $output_folder"
end

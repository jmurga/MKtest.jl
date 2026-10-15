using MKtest, CSV, DataFrames, DelimitedFiles
out = ARGS[1]
tab = CSV.read(joinpath(@__DIR__, "..", "..", "..", "test/fixtures/grapes_parity_tgp.csv"), DataFrame)
sfs, div, _ = MKtest.parse_sfs(joinpath(@__DIR__, "..", "..", "..", "data/MKdata_tgp_may2023.txt"); sample_size = 661)
div = Float64.(vec(div))
writedlm(joinpath(out, "div.tsv"), permutedims(div), '\t')
# the column names are GRAPES' (as in the recorded params.tsv)
rows = String["k\tsource\tnegGshape\tnegGmean\tpos_prop\tposGmean\tsyno_orientation_error"]
for k in (10, 20)
    s = MKtest.project(sfs, k)
    writedlm(joinpath(out, "sfs_k$(k).tsv"), s, '\t')
    j = MKtest.dfe_ml(s, div, MKtest.GammaExpo(); integrator = MKtest.RombergIntegrator())
    push!(rows, join([k, "julia", MKtest.param_values(j, :del_shape)[1], MKtest.param_values(j, :del_mean)[1], MKtest.param_values(j, :ben_prop)[1], MKtest.param_values(j, :ben_mean)[1], MKtest.param_values(j, :misorient)[1]], '\t'))
    g = tab[(tab.n_ind .== k) .& (tab.model .== "GammaExpo"), :][1, :]
    push!(rows, join([k, "grapes", g.p_negGshape, g.p_negGmean, g.p_pos_prop, g.p_posGmean, g.p_syno_orientation_error], '\t'))
end
write(joinpath(out, "params.tsv"), join(rows, "\n") * "\n")

using MKtest, CSV, DataFrames, DelimitedFiles
out = ARGS[1]
tab = CSV.read(expanduser("~/.julia/dev/MKtest_grapes/test/fixtures/grapes_parity_tgp.csv"), DataFrame)
sfs, div, _ = MKtest.parse_sfs(expanduser("~/.julia/dev/MKtest_grapes/data/MKdata_tgp_may2023.txt"); sample_size = 661)
div = Float64.(vec(div))
writedlm(joinpath(out, "div.tsv"), permutedims(div), '\t')
rows = String["k\tsource\tnegGshape\tnegGmean\tpos_prop\tposGmean\tsyno_orientation_error"]
for k in (10, 20)
    s = MKtest.project(sfs, k)
    writedlm(joinpath(out, "sfs_k$(k).tsv"), s, '\t')
    j = MKtest.dfe_ml(s, div, MKtest.GammaExpo(); integrator = MKtest.RombergIntegrator())
    push!(rows, join([k, "julia", j.negGshape[1], j.negGmean[1], j.pos_prop[1], j.posGmean[1], j.syno_orientation_error[1]], '\t'))
    g = tab[(tab.n_ind .== k) .& (tab.model .== "GammaExpo"), :][1, :]
    push!(rows, join([k, "grapes", g.p_negGshape, g.p_negGmean, g.p_pos_prop, g.p_posGmean, g.p_syno_orientation_error], '\t'))
end
write(joinpath(out, "params.tsv"), join(rows, "\n") * "\n")

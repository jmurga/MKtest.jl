# GRAPES binary parity fixture for the dfe_ml engine (Task 1, decision D10).
#
# Oracle: the GRAPES binary (conda package grapes-static, genomedk), at the version
# recorded in each output row. One job per (n, model), so the jobs run in parallel:
#
#   cd ~/.julia/dev/MKtest_grapes
#   for n in 10 20; do for m in GammaZero GammaExpo GammaGamma DisplGamma ScaledBeta; do
#       echo "$n $m"; done; done |
#     xargs -P 10 -n 2 julia --project=. test/fixtures/make_grapes_parity_fixture.jl
#   julia --project=. test/fixtures/make_grapes_parity_fixture.jl merge
#
# Input: data/MKdata_tgp_may2023.txt, parse_sfs(...; sample_size = 661) (1322
# chromosomes), projected by `MKtest.grapes(...; n)` itself to n individuals (2n
# chromosomes), exactly as `MKtest.project(sfs, n)`. The .dofe the binary read is
# re-created and kept under test/fixtures/dofe/ for the record.
#
# Output: test/fixtures/grapes_parity/<n>_<model>.csv (one row), merged into
# test/fixtures/grapes_parity_tgp.csv. The Julia side is NOT run here: the test
# (test/ml/grapes_parity.jl) reads this file, so CI needs no binary.

using MKtest
using CSV
using DataFrames
using Dates

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const TGP = joinpath(ROOT, "data", "MKdata_tgp_may2023.txt")
const JOBDIR = joinpath(@__DIR__, "grapes_parity")
const MERGED = joinpath(@__DIR__, "grapes_parity_tgp.csv")
const MODELS = ("GammaZero", "GammaExpo", "GammaGamma", "DisplGamma", "ScaledBeta")

src_sha() = strip(read(pipeline(`bash -c "cd '$ROOT' && find src ext -type f -name '*.jl' -print0 | sort -z | xargs -0 sha256sum | sha256sum | cut -d' ' -f1"`), String))

function grapes_version()
    meta = joinpath(MKtest.CondaPkg.envdir(), "conda-meta")
    hits = filter(f -> startswith(f, "grapes-static-") && endswith(f, ".json"), readdir(meta))
    isempty(hits) && return "unknown"
    return replace(first(hits), ".json" => "")
end

function write_dofe(path, sfs_n, div, n_ind)
    pn = permutedims(sfs_n[:, 2]); ps = permutedims(sfs_n[:, 3])
    dofe = MKtest.sfs_to_dofe(pn, ps, [div[1]], [div[2]], [div[3]], [div[4]], "1", 2 * n_ind)
    open(path, "w") do io
        print(io, "\n#unfolded\n")     # as MKtest.grapes() writes it
    end
    MKtest.write_files(dofe, path, true)
    return path
end

function job(n_ind::Int, model::String)
    model in MODELS || error("model must be one of $MODELS")
    sfs, div, _ = MKtest.parse_sfs(TGP; sample_size = 661)
    sfs_n = MKtest.project(sfs, n_ind)
    mkpath(joinpath(@__DIR__, "dofe"))
    write_dofe(joinpath(@__DIR__, "dofe", "tgp_n$(n_ind).dofe"), sfs_n, div, n_ind)

    t0 = time()
    df = MKtest.grapes(sfs, div, model; n = n_ind)
    wall = time() - t0
    nrow(df) == 1 || error("GRAPES returned $(nrow(df)) rows for n = $n_ind, $model")

    row = Dict{String,Any}(
        "n_ind" => n_ind, "n_chrom" => 2 * n_ind, "model" => model,
        "lnL" => df[1, "lnL"], "alpha" => df[1, "alpha"], "n_param_grapes" => df[1, "#param"],
        "theta" => df[1, "theta"],
        "grapes_version" => grapes_version(),
        # reconstructed from grapes()'s defaults (tmp paths elided)
        "cmd" => "grapes -in <dofe> -out <out> -model $model -nearly_neutral 5 -FWW_threshold 0.15 " *
                 "-nb_rand_start 0 -anc_to_rec_Ne_ratio 1.0",
        "wall_s" => wall, "src_sha" => src_sha(), "julia" => string(VERSION),
        "date" => string(now()))
    for c in names(df)
        startswith(c, model * ":") || continue
        row["p_" * split(c, ":"; limit = 2)[2]] = df[1, c]
    end
    mkpath(JOBDIR)
    out = joinpath(JOBDIR, "$(n_ind)_$(model).csv")
    CSV.write(out, DataFrame(row))
    println("wrote $out: lnL = $(row["lnL"]), alpha = $(row["alpha"]), $(round(wall; digits = 1)) s")
end

function merge_jobs()
    files = sort(filter(endswith(".csv"), readdir(JOBDIR; join = true)))
    isempty(files) && error("no job outputs in $JOBDIR")
    all_rows = vcat([CSV.read(f, DataFrame) for f in files]...; cols = :union)
    CSV.write(MERGED, all_rows)
    println("wrote $MERGED with $(nrow(all_rows)) rows")
    expected = [(n, m) for n in (10, 20) for m in MODELS]
    got = Set(zip(all_rows.n_ind, all_rows.model))
    missing_cells = [c for c in expected if !(c in got)]
    isempty(missing_cells) || println("MISSING cells (report them, do not raise n): $missing_cells")
end

if length(ARGS) == 1 && ARGS[1] == "merge"
    merge_jobs()
elseif length(ARGS) == 2
    job(parse(Int, ARGS[1]), ARGS[2])
else
    error("usage: make_grapes_parity_fixture.jl <n_individuals> <model> | merge")
end

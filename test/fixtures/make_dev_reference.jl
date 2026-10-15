# Refactor guard for the dfe_ml move: the dev engine at commit 9d752e40 on fixed
# inputs. A frozen record of the code the engine was moved FROM, not a correctness
# oracle (that is the GRAPES binary, make_grapes_parity_fixture.jl).
#
# Runs in the DEV environment, on a dev tree at 9d752e40:
#
#   laptop:      git -C ~/.julia/dev/MKtest.jl worktree add ~/.julia/dev/MKtest_dev9d75 9d752e40
#                rsync -a --delete --exclude .git ~/.julia/dev/MKtest_dev9d75/ $WS:.julia/dev/MKtest_dev9d75/
#   workstation: cd ~/.julia/dev/MKtest_dev9d75 && julia --project=. -e 'using Pkg; Pkg.instantiate()'
#                JULIA_NUM_THREADS=1 julia --project=. \
#                    ~/.julia/dev/MKtest_grapes/test/fixtures/make_dev_reference.jl \
#                    ~/.julia/dev/MKtest_grapes/test/fixtures/dev_reference
#
# Writes into the output folder: inputs_<name>.tsv (the exact SFS used, so the test
# does not depend on dev's projection), inputs_div.tsv, and dev_reference.csv with
#   kind = "objective": -loglik at each model start point (natural values given),
#   kind = "fit":       dfe_ml's loglik, alpha and parameters,
# for 4 models x {RombergIntegrator, ExactSojourn} x 2 inputs.

using MKtest
using DelimitedFiles

outdir = length(ARGS) >= 1 ? ARGS[1] : error("usage: make_dev_reference.jl <output folder>")
mkpath(outdir)

src_sha(root) = strip(read(pipeline(`bash -c "cd '$root' && find src ext -type f -name '*.jl' -print0 | sort -z | xargs -0 sha256sum | sha256sum | cut -d' ' -f1"`), String))
const DEV_ROOT = normpath(joinpath(pathof(MKtest), "..", ".."))

js = collect(1.0:39.0)
inputs = Dict{String,Tuple{Matrix{Float64},Vector{Float64}}}()
inputs["synthetic_n40"] = (hcat(js, round.(2000.0 ./ js), round.(3000.0 ./ js)), [900.0, 1200.0, 1.0e6, 3.0e5])
sfs, div, _ = MKtest.parse_sfs(joinpath(DEV_ROOT, "data", "MKdata_tgp_may2023.txt"); sample_size = 661)
inputs["tgp_n20"] = (MKtest.project(sfs, 20), Float64.(vec(div)))

models = [("GammaZero", MKtest.GammaZero()), ("GammaExpo", MKtest.GammaExpo()),
          ("DisplGamma", MKtest.DisplGamma()), ("WithLethal{GammaZero}", MKtest.WithLethal(MKtest.GammaZero()))]
integs = [("RombergIntegrator", MKtest.RombergIntegrator()), ("ExactSojourn", MKtest.ExactSojourn())]

fmt(x) = repr(Float64(x))
lines = ["input\tmodel\tintegrator\tkind\tpoint\tloglik\talpha\ttheta\tparams"]
for (iname, (s, d)) in sort(collect(inputs); by = first)
    writedlm(joinpath(outdir, "inputs_$(iname).tsv"), s, '\t')
    writedlm(joinpath(outdir, "inputs_div_$(iname).tsv"), permutedims(d), '\t')
    n = Int(round(maximum(s[:, 1]))) + 1
    for (mname, m) in models, (gname, integ) in integs
        pre = MKtest.get_precomputed(n, integ; expectation = MKtest.PRFSojourn())
        data = MKtest.build_data_input(s, reshape(d, 1, 4), false, false; neutral_shape = pre.neutral_shape)
        spec = MKtest.build_model_spec(m)
        tr = MKtest.build_transform(spec, Dict{String,Float64}(), true)
        obj = MKtest.make_objective(data, pre, spec, tr, false, false)
        for (k, sp) in enumerate(spec.start_points)
            nat = vcat(sp, 0.01)                      # + syno_orientation_error
            x = MKtest.to_opt_space(nat, tr.opt_log_scale)
            params = join(["$(nm)=$(fmt(v))" for (nm, v) in zip(tr.all_names, nat)], ";")
            push!(lines, join([iname, mname, gname, "objective", k, fmt(-obj(x)), "NaN", "NaN", params], '\t'))
        end
        t0 = time()
        r = MKtest.dfe_ml(s, d, m; integrator = integ)
        params = join(["$(nm)=$(fmt(r.params[nm]))" for nm in tr.all_names], ";")
        push!(lines, join([iname, mname, gname, "fit", 0, fmt(r.loglik), fmt(r.alpha), fmt(r.theta), params], '\t'))
        println("$iname $mname $gname: lnL $(r.loglik) alpha $(r.alpha) ($(round(time() - t0; digits = 1)) s)")
    end
end
open(joinpath(outdir, "dev_reference.csv"), "w") do io
    println(io, "# dev dfe_ml reference; dev src sha $(src_sha(DEV_ROOT)); julia $(VERSION); threads $(Threads.nthreads())")
    foreach(l -> println(io, l), lines)
end
println("wrote $(joinpath(outdir, "dev_reference.csv"))")

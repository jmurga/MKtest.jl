# MKtest test suite. `Pkg.test()` runs every group; `Pkg.test(test_args = ["abc"])` one
# of them. `scripts/test_parallel.sh` runs the groups as parallel single-thread processes
# (most of the suite's time is compilation, which every process then pays only for what it
# uses). The frozen fixtures reproduce at one thread.
using MKtest
using Test
using DataFrames: nrow, names, propertynames, Not
using JLD2
using Random: seed!, MersenneTwister
using DelimitedFiles: readdlm, writedlm
import Logging

quiet(f) = Logging.with_logger(f, Logging.NullLogger())
const BEST = JLD2.load(joinpath(@__DIR__, "fixtures", "abcmk_best_model.jld2"))
P500() = MKtest.parameters(N = 500, n = 50, dac = [2, 4, 5, 10, 20, 50])

const GROUPS = ["rates", "abc", "mk", "engine", "extension"]
const SELECTED = isempty(ARGS) ? GROUPS : ARGS
let unknown = setdiff(SELECTED, GROUPS)
    isempty(unknown) || error("unknown test group(s) $unknown; known: $GROUPS")
end

"rates" in SELECTED && include("abcmk_rates.jl")       # the forward model and its table
"abc" in SELECTED && include("abcmk_abc.jl")           # samplers, priors, ABC, in-memory route
"mk" in SELECTED && include("mk.jl")                   # parse_sfs, MK tests, bootstrap, small methods
"engine" in SELECTED && include("ml/engine.jl")        # dfe_ml engine (src/ml/)
if "extension" in SELECTED
    include("ml/extension.jl")
    include("ml/grapes_parity.jl")
end

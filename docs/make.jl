using MKtest
using Documenter
ENV["GKSwstype"] = "100"   # GR without a display, for the figures of the @example blocks

DocMeta.setdocmeta!(MKtest, :DocTestSetup, :(using MKtest); recursive=true)

makedocs(;
    modules=[MKtest],
    authors="Jesús Murga-Moreno, Lawrence Uricchio, David Enard",
    repo="https://github.com/jmurga/MKtest.jl/blob/{commit}{path}#{line}",
    sitename="MKtest.jl",
    checkdocs=:exports,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jmurga.github.io/MKtest.jl",
        assets=String[],
        size_threshold=400 * 2^10,         # the reference page lists the whole API
        size_threshold_warn=250 * 2^10,
    ),
        pages     = [
                "Home" => "index.md",
                "Input data" => "data.md",
                "MK tests" => "mk.md",
                "DFE and alpha by maximum likelihood" => "dfe_ml.md",
                "How to add a new DFE" => "new_dfe.md",
                "ABC-MK" => "abcmk.md",
                "Case and control sets" => "case_control.md",
                "Reference" => "reference.md"
        ]
)

deploydocs(;
    repo="github.com/jmurga/MKtest.jl",
    devbranch="main",
)

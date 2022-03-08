using MKtest
using Documenter

DocMeta.setdocmeta!(MKtest, :DocTestSetup, :(using MKtest); recursive=true)

makedocs(;
    modules=[MKtest],
    authors="Jesús Murga-Moreno, Lawrence Uricchio, David Enard",
    repo="https://github.com/jmurga/MKtest.jl/blob/{commit}{path}#{line}",
    sitename="MKtest.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jmurga.github.io/MKtest.jl",
        assets=String[],
    ),
        pages     = [
                "Home" => "index.md",
                "Package Overview" => [
                        "Analytical estimations" => "analytical.md",
                        "Processing data" => "data.md",
                        "MK approaches" => "mk.md"
                ],
                "Infering the rate and strength of adaptation" =>[
                        "Empirical estimation" => "empirical.md",
                        "Rates" => "rates.md",
                        "Input data" => "input.md",
                        "Summary statistics" => "summstat.md",
                        "ABC inference" => "abc.md",
                        "Multiple dataset" => "multiple.md",
                ],
                "Command-Line interface" => "cli.md",
                "Reference" => "reference.md"
        ]

)

deploydocs(;
    repo="github.com/jmurga/MKtest.jl",
    devbranch="main",
)

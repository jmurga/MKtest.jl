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
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jmurga/MKtest.jl",
    devbranch="main",
)

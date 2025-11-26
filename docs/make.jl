using TimeSeriesUtils
using Documenter

DocMeta.setdocmeta!(TimeSeriesUtils, :DocTestSetup, :(using TimeSeriesUtils); recursive=true)

makedocs(;
    modules=[TimeSeriesUtils],
    authors="Mattias Villani",
    sitename="TimeSeriesUtils.jl",
    format=Documenter.HTML(;
        canonical="https://mattiasvillani.github.io/TimeSeriesUtils.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mattiasvillani/TimeSeriesUtils.jl",
    devbranch="main",
)

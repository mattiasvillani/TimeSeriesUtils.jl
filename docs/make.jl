using Pkg

# Activate the docs environment in docs/Project.toml
Pkg.activate(@__DIR__)
Pkg.instantiate()  # optional but good to keep

@show Base.active_project()  # temporary debug, can be removed later

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
        "ARIMA utilities" => "ARIMAUtils.md",
        "ARTFIMA utilities" => "ARTFIMAUtils.md",
        "Spectral utilities" => "SpectralUtils.md",
        "General utilities" => "GeneralUtils.md"
    ],
)

deploydocs(;
    repo="github.com/mattiasvillani/TimeSeriesUtils.jl",
    devbranch="main",
)

using QuantumMeasurements
using Documenter

DocMeta.setdocmeta!(QuantumMeasurements, :DocTestSetup, :(using QuantumMeasurements); recursive=true)

makedocs(;
    modules=[QuantumMeasurements],
    authors="Mario Dagrada <mario.dagrada@pasqal.com>, " *
            "Sebastian Grijalva <sebastian.grijalva@pasqal.com>, " *
            "Niklas Heim <heim.niklas@pasqal.com>, " *
            "Louis Vignoli <louis.vignoli@pasqal.com>",
    repo="https://github.com/pasqal-io/QuantumMeasurements.jl/blob/{commit}{path}#{line}",
    sitename="QuantumMeasurements.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pasqal-io.github.io/QuantumMeasurements.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pasqal-io/QuantumMeasurements.jl",
    devbranch="main",
    push_preview=true
)

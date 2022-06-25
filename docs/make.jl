using Sniper
using Documenter

DocMeta.setdocmeta!(Sniper, :DocTestSetup, :(using Sniper); recursive=true)

makedocs(;
    modules=[Sniper],
    authors="Niklas Heim <heim.niklas@gmail.com>",
    repo="https://github.com/pasqal/Sniper.jl/blob/{commit}{path}#{line}",
    sitename="Sniper.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pasqal.github.io/Sniper.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pasqal/Sniper.jl",
    devbranch="main",
)
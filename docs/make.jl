using Sniper
using Documenter

DocMeta.setdocmeta!(Sniper, :DocTestSetup, :(using Sniper); recursive=true)

makedocs(;
    modules=[Sniper],
    authors=[
        "Mario Dagrada <mario.dagrada@pasqal.com>",
        "Sebastian Grijalva <sebastian.grijalva@pasqal.com>",
        "Niklas Heim <heim.niklas@gmail.com>",
        "Louis Vignoli <louis.vignoli@pasqal.com>",
    ],
    repo="https://github.com/pasqal-io/Sniper.jl/blob/{commit}{path}#{line}",
    sitename="Sniper.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pasqal-io.github.io/Sniper.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pasqal-io/Sniper.jl",
    devbranch="main",
    push_preview=true
)

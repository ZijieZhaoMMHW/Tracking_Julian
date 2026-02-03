using ExtremeTracking
using Documenter

DocMeta.setdocmeta!(ExtremeTracking, :DocTestSetup, :(using ExtremeTracking); recursive=true)

makedocs(;
    modules=[ExtremeTracking],
    authors="gaelforget <gforget@mit.edu> and contributors",
    sitename="ExtremeTracking.jl",
    format=Documenter.HTML(;
        canonical="https://gaelforget.github.io/ExtremeTracking.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/gaelforget/ExtremeTracking.jl",
    devbranch="main",
)

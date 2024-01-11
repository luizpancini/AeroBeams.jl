using AeroBeams
using Documenter

DocMeta.setdocmeta!(AeroBeams, :DocTestSetup, :(using AeroBeams); recursive=true)

push!(LOAD_PATH,"../src/")
makedocs(;
    modules=[AeroBeams],
    authors="Luiz Pancini <luizpancini@usp.br>",
    sitename="AeroBeams.jl",
    format=Documenter.HTML(;
        canonical="https://luizpancini.github.io/AeroBeams.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "References" => "references.md"
    ],
)

deploydocs(;
    repo="github.com/luizpancini/AeroBeams.jl.git",
    devbranch="main",
)

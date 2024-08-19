using AeroBeams, Documenter, DocumenterCitations

# DocMeta.setdocmeta!(AeroBeams, :DocTestSetup, :(using AeroBeams); recursive=true)

bib_filepath = joinpath(@__DIR__, "src/references.bib")
bib = CitationBibliography(bib_filepath, style=:authoryear)

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
        "Home" => "home.md",
        "API" => [
            "Public" => "public.md"
            "Private" => "private.md"
            ]
    ],
    plugins = [bib]
)

deploydocs(;
    repo="github.com/luizpancini/AeroBeams.jl.git",
    devbranch="main",
)

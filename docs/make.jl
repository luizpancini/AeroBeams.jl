using AeroBeams, Documenter, Literate

# Pre-install matplotlib
import Plots; Plots.pyplot()

DocMeta.setdocmeta!(AeroBeams, :DocTestSetup, :(using AeroBeams); recursive=true)

# Examples to be included in the documentation
const included = ["archUnderDeadPressure.jl","archUnderFollowerPressure.jl"]

# Literate .md files output
for ex in included
    inputPath = pkgdir(AeroBeams)*"/test/examples/"*ex
    outputPath = "src/"
    exName = splitext(ex)[1]
    Literate.markdown(inputPath, outputPath; name=exName, mdstrings=true)
end

# Make documentation
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
        "Examples" => ["Arch under dead pressure" => "archUnderDeadPressure.md",
        "Arch under follower pressure" => "archUnderFollowerPressure.md"],
        "Public API" => "publicAPI.md"
    ],
)

# CI
deploydocs(;
    repo="github.com/luizpancini/AeroBeams.jl.git",
    devbranch="main",
)
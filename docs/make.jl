using AeroBeams, Documenter, Literate

# Set plot backend for documentation
import Plots; Plots.pyplot()

DocMeta.setdocmeta!(AeroBeams, :DocTestSetup, :(using AeroBeams); recursive=true)

# Examples to be included in the documentation
global included = ["archUnderFollowerPressure.jl","BWBflutter.jl","conventionalHALECheckedPitchManeuver.jl","conventionalHALEmodel.jl","initialDispAndVelBeam.jl","twoStoryFrame.jl"]

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
        "Examples" => [
        "Arch under follower pressure" => "archUnderFollowerPressure.md",
        "Creating a HALE aircraft model" => "conventionalHALEmodel.md",
        "Motion of a simply supported beam under initial conditions" => "initialDispAndVelBeam.md",
        "Flutter of a Blended-Wing-Body" => "BWBflutter.md",
        "Pitch maneuver of a HALE aircraft" => "conventionalHALECheckedPitchManeuver.md",
        "Two-story frame" => "twoStoryFrame.md"
        ],
        "Public API" => "publicAPI.md"
    ],
)

# CI
deploydocs(;
    repo="github.com/luizpancini/AeroBeams.jl.git",
    devbranch="main",
)
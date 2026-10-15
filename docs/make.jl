using AeroBeams, Documenter, Literate

# Headless GR backend, so that plots in the examples are rendered without opening windows
ENV["GKSwstype"] = "100"

# Examples to be included in the documentation (paths relative to test/examples)
included = [
    "aerodynamic/OMCgustTests.jl",
    "dynamicAeroelastic/cHALECheckedPitchManeuver.jl",
    "dynamicAeroelastic/cHALECheckedRollManeuver.jl",
    "dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating.jl",
    "dynamicStructural/flyingScissors.jl",
    "dynamicStructural/initialDispAndVelBeam.jl",
    "eigenAeroelastic/BWBflutter.jl",
    "eigenAeroelastic/HealyBaselineFFWTfreeFlutterAoARangeURange.jl",
    "eigenAeroelastic/PazyWingFlutterPitchRange.jl",
    "eigenAeroelastic/typicalSectionFlutterAndDivergence.jl",
    "eigenStructural/sweptTipRotor.jl",
    "eigenStructural/twoStoryFrame.jl",
    "staticStructural/archUnderFollowerPressure.jl",
    "staticStructural/compositeCantileverMD.jl",
    "staticStructural/PazyWingTorsionTest.jl",
    "steadyAeroelastic/cHALEmodel.jl",
    "steadyAeroelastic/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast.jl",
    "steadyAeroelastic/PazyWingPitchRange.jl",
    "trimAeroelastic/heliosTrim.jl",
]

# The "## Plots" section at the end of each example is for local use only: remove it from the docs
hide_local_plots(content) = replace(content, r"^## Plots\b.*\z"ms => "")

# Generate .md files with Literate
outputPath = joinpath(@__DIR__, "src", "literate")
for ex in included
    inputPath = joinpath(pkgdir(AeroBeams), "test", "examples", ex)
    exName = splitext(basename(ex))[1]
    Literate.markdown(inputPath, outputPath; name=exName, mdstrings=true, preprocess=hide_local_plots)
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
        size_threshold_ignore = ["publicAPI.md"]
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Examples" => [
        "Static structural analysis of an arch under follower pressure" => "literate/archUnderFollowerPressure.md",
        "Static structural analysis of composite beams" => "literate/compositeCantileverMD.md",
        "Static structural analysis of the Pazy wing" => "literate/PazyWingTorsionTest.md",
        "Modal analysis of a rotating beam" => "literate/sweptTipRotor.md",
        "Modal analysis of a two-story frame" => "literate/twoStoryFrame.md",
        "Motion of a simply supported beam under initial conditions" => "literate/initialDispAndVelBeam.md",
        "Dynamic analysis of an articulated beam" => "literate/flyingScissors.md",
        "Flutter and divergence of a typical aeroelastic section" => "literate/typicalSectionFlutterAndDivergence.md",
        "Gust response of an airfoil section" => "literate/OMCgustTests.md",
        "Steady aeroelastic analysis of the Pazy Wing" => "literate/PazyWingPitchRange.md",
        "Steady aeroelastic analysis of a wing with flared folding wingtip (FFWT)" => "literate/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast.md",
        "Flutter analysis of the Pazy wing" => "literate/PazyWingFlutterPitchRange.md",
        "Flutter analysis of a wing with flared folding wingtip (FFWT)" => "literate/HealyBaselineFFWTfreeFlutterAoARangeURange.md",
        "Gust response of a wing with flared folding wingtip (FFWT)" => "literate/HealyBaselineFFWTOMCGustFloating.md",
        "Creating a HALE aircraft model" => "literate/cHALEmodel.md",
        "Trimming a flying-wing HALE" => "literate/heliosTrim.md",
        "Flutter of a Blended-Wing-Body" => "literate/BWBflutter.md",
        "Pitch maneuver of a HALE aircraft" => "literate/cHALECheckedPitchManeuver.md",
        "Roll maneuver of a HALE aircraft" => "literate/cHALECheckedRollManeuver.md"
        ],
        "Public API" => "publicAPI.md"
    ],
)

# CI
deploydocs(;
    repo="github.com/luizpancini/AeroBeams.jl.git",
    devbranch="main",
)
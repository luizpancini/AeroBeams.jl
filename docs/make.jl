using AeroBeams, Documenter, Literate

# Examples to be included in the documentation
global included = ["archUnderFollowerPressure.jl","BWBflutter.jl","compositeCantileverMD.jl","conventionalHALECheckedPitchManeuver.jl","conventionalHALEmodel.jl","flyingScissors.jl","HealyBaselineFFWTfreeFlutterAoARangeURange.jl","HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast.jl","heliosTrim.jl","initialDispAndVelBeam.jl","OMCgustTests.jl","PazyWingPitchRange.jl","PazyWingFlutterPitchRange.jl","PazyWingTorsionTest.jl","sweptTipRotor.jl","twoStoryFrame.jl","typicalSectionFlutterAndDivergence.jl"]

# Generate .md files with Literate
for ex in included
    inputPath = pkgdir(AeroBeams)*"/test/examples/"*ex
    outputPath = "src/literate/"
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
        size_threshold_ignore = ["publicAPI.md"]
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Examples" => [
        "Creating a HALE aircraft model" => "literate/conventionalHALEmodel.md",
        "Gust response of an airfoil section" => "literate/OMCgustTests.md",
        "Flutter and divergence of a typical aeroelastic section" => "literate/typicalSectionFlutterAndDivergence.md",
        "Static structural analysis of an arch under follower pressure" => "literate/archUnderFollowerPressure.md",
        "Static structural analysis of composite beams" => "literate/compositeCantileverMD.md",
        "Static structural analysis of the Pazy wing" => "literate/PazyWingTorsionTest.md",
        "Steady aeroelastic analysis of the Pazy Wing" => "literate/PazyWingPitchRange.md",
        "Flutter analysis of a wing with flared folding wingtip (FFWT)" => "literate/HealyBaselineFFWTfreeFlutterAoARangeURange.md",
        "Steady aeroelastic analysis of a wing with flared folding wingtip (FFWT)" => "literate/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast.md",
        "Flutter analysis of the Pazy wing" => "literate/PazyWingFlutterPitchRange.md",
        "Modal analysis of a rotating beam" => "literate/sweptTipRotor.md",
        "Modal analysis of a two-story frame" => "literate/twoStoryFrame.md",
        "Motion of a simply supported beam under initial conditions" => "literate/initialDispAndVelBeam.md",
        "Dynamic analysis of an articulated beam" => "literate/flyingScissors.md",
        "Trimming a flying-wing HALE" => "literate/heliosTrim.md",
        "Flutter of a Blended-Wing-Body" => "literate/BWBflutter.md",
        "Pitch maneuver of a HALE aircraft" => "literate/conventionalHALECheckedPitchManeuver.md"
        ],
        "Public API" => "publicAPI.md"
    ],
)

# CI
deploydocs(;
    repo="github.com/luizpancini/AeroBeams.jl.git",
    devbranch="main",
)
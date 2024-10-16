# New reference data for eigen aeroelastic tests

# Flutter analysis of the Blended-Wing-Body flying wing
include("../examples/BWBflutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/BWBflutter"))
writedlm("test/newTestDataGenerators/BWBflutter/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/BWBflutter/damps.txt", damps)

# Flutter analysis of the conventional HALE aircraft in free flight with structural stiffness as the varying parameter
include("../examples/conventionalHALELambdaRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALELambdaRange"))
writedlm("test/newTestDataGenerators/conventionalHALELambdaRange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/conventionalHALELambdaRange/damps.txt", damps)

# Flutter analysis of the conventional HALE aircraft in free flight with airspeed and structural stiffness as the varying parameters
include("../examples/conventionalHALELURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALELURange"))
for (i,λ) in enumerate(λRange)
    writedlm(string("test/newTestDataGenerators/conventionalHALELURange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/conventionalHALELURange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Flutter analysis of the conventional HALE aircraft in free flight with airspeed as the varying parameter
include("../examples/conventionalHALEURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALEURange"))
writedlm("test/newTestDataGenerators/conventionalHALEURange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/conventionalHALEURange/damps.txt", damps)

# Flutter analysis of the Goland wing
include("../examples/GolandWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/GolandWingFlutter"))
writedlm("test/newTestDataGenerators/GolandWingFlutter/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/GolandWingFlutter/damps.txt", damps)

# Flutter analysis of the Helios flying-wing in free flight with payload and structural stiffness as the varying parameters
include("../examples/heliosFlutterPLambdaRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosFlutterPLambdaRange"))
for (i,λ) in enumerate(λRange)
    writedlm(string("test/newTestDataGenerators/heliosFlutterPLambdaRange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/heliosFlutterPLambdaRange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Flutter analysis of the Helios flying-wing in free flight with payload as the varying parameter
include("../examples/heliosFlutterPRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosFlutterPRange"))
writedlm("test/newTestDataGenerators/heliosFlutterPRange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/heliosFlutterPRange/damps.txt", damps)

# Flutter analysis of the Helios flying-wing in free flight with airspeed as the varying parameter
include("../examples/heliosFlutterURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosFlutterURange"))
writedlm("test/newTestDataGenerators/heliosFlutterURange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/heliosFlutterURange/damps.txt", damps)

# Flutter analysis of the wing of the Helios flying-wing
include("../examples/heliosWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosWingFlutter"))
writedlm("test/newTestDataGenerators/heliosWingFlutter/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/heliosWingFlutter/damps.txt", damps)

# Eigen-analysis of the Pazy wing with flared folding wing tip (FFWT)
include("../examples/PazyFFWTeigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTeigen"))
writedlm("test/newTestDataGenerators/PazyFFWTeigen/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/PazyFFWTeigen/damps.txt", damps)

# Flutter analysis of the Pazy wing
include("../examples/PazyWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutter"))
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOnsetSpeedsOfMode.txt", flutterOnsetSpeedsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOnsetFreqsOfMode.txt", flutterOnsetFreqsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOnsetDispOfMode.txt", flutterOnsetDispOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOffsetSpeedsOfMode.txt", flutterOffsetSpeedsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOffsetFreqsOfMode.txt", flutterOffsetFreqsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOffsetDispOfMode.txt", flutterOffsetDispOfMode)

# Flutter and divergence analysis of the Pazy wing
include("../examples/PazyWingFlutterAndDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutterAndDivergence"))
writedlm("test/newTestDataGenerators/PazyWingFlutterAndDivergence/flutterOnsetSpeedsOfMode.txt", flutterOnsetSpeedsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutterAndDivergence/flutterOnsetFreqsOfMode.txt", flutterOnsetFreqsOfMode)

# Flutter analysis of the Pazy wing with varying root pitch angle
include("../examples/PazyWingFlutterPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutterPitchRange"))
for (i,θ) in enumerate(θRange)
    writedlm(string("test/newTestDataGenerators/PazyWingFlutterPitchRange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/PazyWingFlutterPitchRange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Flutter analysis of the Pazy wing with varying tip mass positions
include("../examples/PazyWingFlutterTipMassRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutterTipMassRange"))
for i=1:3
    writedlm(string("test/newTestDataGenerators/PazyWingFlutterTipMassRange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/PazyWingFlutterTipMassRange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Flutter analysis of the sixteen-meter-wing
include("../examples/SMWFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutter"))
writedlm("test/newTestDataGenerators/SMWFlutter/x1_def.txt",x1_def[end]/L)
writedlm("test/newTestDataGenerators/SMWFlutter/x3_def.txt",x3_def[end]/L)
writedlm("test/newTestDataGenerators/SMWFlutter/alpha.txt",α_of_x1[end]*180/pi)
for mode in 1:nModes
    writedlm(string("test/newTestDataGenerators/SMWFlutter/freqsMode",mode,".txt"), modeFrequencies[mode])
    writedlm(string("test/newTestDataGenerators/SMWFlutter/dampsMode",mode,".txt"), modeDampings[mode])
end

# Flutter boundary analysis of the sixteen-meter-wing as a function of the pitch angle
include("../examples/SMWFlutterPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterPitchRange"))
for i in 1:nModes
    writedlm(string("test/newTestDataGenerators/SMWFlutterPitchRange/freqsMode",i,".txt"),modeFrequencies[i])
    writedlm(string("test/newTestDataGenerators/SMWFlutterPitchRange/dampsMode",i,".txt"),modeDampings[i])
end

# Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with tip load as the varying parameter
include("../examples/SMWFlutterPrecurvatureRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterPrecurvatureRange"))
for i in eachindex(kRange)
    writedlm(string("test/newTestDataGenerators/SMWFlutterPrecurvatureRange/flutterSpeedk",i,".txt"),flutterSpeed[i,:])
end

# Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with root angle as the varying parameter
include("../examples/SMWFlutterPrecurvatureRange2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterPrecurvatureRange2"))
for (ki,k) in enumerate(kRange)
    for (i,θ) in enumerate(θRange)
        writedlm(string("test/newTestDataGenerators/SMWFlutterPrecurvatureRange2/freqs_k",ki,"th",i,".txt"),hcat(freqs[ki,i,:]...)')
        writedlm(string("test/newTestDataGenerators/SMWFlutterPrecurvatureRange2/damps_k",ki,"th",i,".txt"),hcat(damps[ki,i,:]...)')
    end
end

# Flutter boundary analysis of the sixteen-meter-wing as a function of the bending-torsion coupling factor
include("../examples/SMWFlutterStructuralCouplingRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterStructuralCouplingRange"))
writedlm("test/newTestDataGenerators/SMWFlutterStructuralCouplingRange/flutterSpeed.txt", flutterSpeed)

# Flutter boundary analysis of the sixteen-meter-wing as a function of the tip displacement
include("../examples/SMWFlutterTipDispRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterTipDispRange"))
writedlm("test/newTestDataGenerators/SMWFlutterTipDispRange/flutterSpeed.txt", flutterSpeed)
writedlm("test/newTestDataGenerators/SMWFlutterTipDispRange/flutterFreq.txt", flutterFreq)

# Linear flutter analysis of the sixteen-meter-wing
include("../examples/SMWLinearFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWLinearFlutter"))
writedlm("test/newTestDataGenerators/SMWLinearFlutter/flutterSpeed.txt", flutterSpeed)
writedlm("test/newTestDataGenerators/SMWLinearFlutter/flutterFreq.txt", flutterFreq)

# Flutter and divergence analysis of a typical section
include("../examples/typicalSectionFlutterAndDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/typicalSectionFlutterAndDivergence"))
writedlm("test/newTestDataGenerators/typicalSectionFlutterAndDivergence/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/typicalSectionFlutterAndDivergence/damps.txt", damps)
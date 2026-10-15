using Random

# Flutter analysis of the Blended-Wing-Body flying wing
include("../examples/eigenAeroelastic/BWBflutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/BWBflutter"))
writedlm("test/newTestDataGenerators/BWBflutter/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/BWBflutter/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the matched wing of the Blended-Wing-Body flying wing
include("../examples/eigenAeroelastic/BWBmatchedWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/BWBmatchedWingFlutter"))
writedlm("test/newTestDataGenerators/BWBmatchedWingFlutter/freqsOsc.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/BWBmatchedWingFlutter/dampsOsc.txt", test_data_matrix(dampsOsc))
writedlm("test/newTestDataGenerators/BWBmatchedWingFlutter/trimAoA.txt", test_data_matrix(trimAoA))
writedlm("test/newTestDataGenerators/BWBmatchedWingFlutter/trimThrust.txt", test_data_matrix(trimThrust))
writedlm("test/newTestDataGenerators/BWBmatchedWingFlutter/trimDelta.txt", test_data_matrix(trimδ))
writedlm("test/newTestDataGenerators/BWBmatchedWingFlutter/flutterSpeeds.txt", test_data_matrix(flutterSpeeds))

# Flutter analysis of the conventional HALE aircraft in free flight with structural stiffness as the varying parameter
include("../examples/eigenAeroelastic/cHALELambdaRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cHALELambdaRange"))
writedlm("test/newTestDataGenerators/cHALELambdaRange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/cHALELambdaRange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the conventional HALE aircraft in free flight with airspeed and structural stiffness as the varying parameters
include("../examples/eigenAeroelastic/cHALELURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cHALELURange"))
writedlm(string("test/newTestDataGenerators/cHALELURange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/cHALELURange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the conventional HALE aircraft in free flight with airspeed as the varying parameter
include("../examples/eigenAeroelastic/cHALEURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cHALEURange"))
writedlm("test/newTestDataGenerators/cHALEURange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/cHALEURange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the Goland wing
include("../examples/eigenAeroelastic/GolandWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/GolandWingFlutter"))
writedlm("test/newTestDataGenerators/GolandWingFlutter/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/GolandWingFlutter/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the baseline Healy free FFWT wing with tip loss, root pitch and airspeed as the varying parameters
include("../examples/eigenAeroelastic/HealyBaselineFFWTfreeFlutterAoARangeURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTfreeFlutterAoARangeURange"))
writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTfreeFlutterAoARangeURange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTfreeFlutterAoARangeURange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the baseline Healy free FFWT wing with flare angle and airspeed as the varying parameters
include("../examples/eigenAeroelastic/HealyBaselineFFWTfreeFlutterFlareRangeURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTfreeFlutterFlareRangeURange"))
writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTfreeFlutterFlareRangeURange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTfreeFlutterFlareRangeURange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the baseline Healy locked FFWT wing with airspeed as the varying parameters
include("../examples/eigenAeroelastic/HealyBaselineFFWTlockedFlutterURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTlockedFlutterURange"))
writedlm("test/newTestDataGenerators/HealyBaselineFFWTlockedFlutterURange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/HealyBaselineFFWTlockedFlutterURange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the Healy LCO FFWT wing with root pitch angle as the varying parameter
include("../examples/eigenAeroelastic/HealyLCOFFWTfreeFlutterAoARange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyLCOFFWTfreeFlutterAoARange"))
writedlm(string("test/newTestDataGenerators/HealyLCOFFWTfreeFlutterAoARange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/HealyLCOFFWTfreeFlutterAoARange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the Healy sideslip FFWT wing with sideslip angle as the varying parameter
include("../examples/eigenAeroelastic/HealySideslipFFWTflutterSideslipRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealySideslipFFWTflutterSideslipRange"))
writedlm("test/newTestDataGenerators/HealySideslipFFWTflutterSideslipRange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/HealySideslipFFWTflutterSideslipRange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the Healy sideslip FFWT wing with airspeed as the varying parameter
include("../examples/eigenAeroelastic/HealySideslipFFWTflutterURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealySideslipFFWTflutterURange"))
writedlm("test/newTestDataGenerators/HealySideslipFFWTflutterURange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/HealySideslipFFWTflutterURange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the Helios flying-wing in free flight with payload and structural stiffness as the varying parameters
Random.seed!(1234)
include("../examples/eigenAeroelastic/heliosFlutterPLambdaRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosFlutterPLambdaRange"))
writedlm(string("test/newTestDataGenerators/heliosFlutterPLambdaRange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/heliosFlutterPLambdaRange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the Helios flying-wing in free flight with payload as the varying parameter
Random.seed!(1234)
include("../examples/eigenAeroelastic/heliosFlutterPRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosFlutterPRange"))
writedlm("test/newTestDataGenerators/heliosFlutterPRange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/heliosFlutterPRange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the Helios flying-wing in free flight with airspeed as the varying parameter
include("../examples/eigenAeroelastic/heliosFlutterURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosFlutterURange"))
writedlm("test/newTestDataGenerators/heliosFlutterURange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/heliosFlutterURange/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the wing of the Helios flying-wing
include("../examples/eigenAeroelastic/heliosWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosWingFlutter"))
writedlm("test/newTestDataGenerators/heliosWingFlutter/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/heliosWingFlutter/damps.txt", test_data_matrix(dampsOsc))

# Flutter analysis of the Pazy FFWT wing with airspeed as the varying parameter
include("../examples/eigenAeroelastic/PazyFFWTflutterURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTflutterURange"))
writedlm("test/newTestDataGenerators/PazyFFWTflutterURange/freqs.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/PazyFFWTflutterURange/damps.txt", test_data_matrix(dampsOsc))
writedlm("test/newTestDataGenerators/PazyFFWTflutterURange/phiHinge.txt", test_data_matrix(ϕHinge))

# Flutter analysis of the Pazy wing
include("../examples/eigenAeroelastic/PazyWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutter"))
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOnsetSpeeds.txt", test_data_matrix(flutterOnsetSpeeds))
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOnsetFreqs.txt", test_data_matrix(flutterOnsetFreqs))
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOnsetDisp.txt", test_data_matrix(flutterOnsetDisp))
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOffsetSpeeds.txt", test_data_matrix(flutterOffsetSpeeds))
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOffsetFreqs.txt", test_data_matrix(flutterOffsetFreqs))
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOffsetDisp.txt", test_data_matrix(flutterOffsetDisp))

# Flutter and divergence analysis of the Pazy wing
include("../examples/eigenAeroelastic/PazyWingFlutterAndDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutterAndDivergence"))
writedlm("test/newTestDataGenerators/PazyWingFlutterAndDivergence/flutterOnsetSpeeds.txt", test_data_matrix(flutterOnsetSpeeds))
writedlm("test/newTestDataGenerators/PazyWingFlutterAndDivergence/flutterOnsetFreqs.txt", test_data_matrix(flutterOnsetFreqs))

# Flutter analysis of the Pazy wing with varying root pitch angle
include("../examples/eigenAeroelastic/PazyWingFlutterPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutterPitchRange"))
writedlm(string("test/newTestDataGenerators/PazyWingFlutterPitchRange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/PazyWingFlutterPitchRange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the Pazy wing with varying tip mass positions
include("../examples/eigenAeroelastic/PazyWingFlutterTipMassRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutterTipMassRange"))
writedlm(string("test/newTestDataGenerators/PazyWingFlutterTipMassRange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/PazyWingFlutterTipMassRange/damps.txt"), test_data_matrix(dampsOsc))

# Frequency analysis of the Pazy wing
include("../examples/eigenAeroelastic/PazyWingFreqsEvolution.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFreqsEvolution"))
writedlm("test/newTestDataGenerators/PazyWingFreqsEvolution/freqs.txt", test_data_matrix(freqsOsc))

# Flutter analysis of the sixteen-meter-wing
include("../examples/eigenAeroelastic/SMWFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutter"))
writedlm("test/newTestDataGenerators/SMWFlutter/x1_def.txt", test_data_matrix(x1_def[end]/L))
writedlm("test/newTestDataGenerators/SMWFlutter/x3_def.txt", test_data_matrix(x3_def[end]/L))
writedlm("test/newTestDataGenerators/SMWFlutter/alpha.txt", test_data_matrix(α_of_x1[end]))
writedlm(string("test/newTestDataGenerators/SMWFlutter/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/SMWFlutter/damps.txt"), test_data_matrix(dampsOsc))

# Flutter boundary analysis of the sixteen-meter-wing as a function of the pitch angle
include("../examples/eigenAeroelastic/SMWFlutterPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterPitchRange"))
writedlm(string("test/newTestDataGenerators/SMWFlutterPitchRange/freqs.txt"),test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/SMWFlutterPitchRange/damps.txt"),test_data_matrix(dampsOsc))

# Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with tip load as the varying parameter
include("../examples/eigenAeroelastic/SMWFlutterPrecurvatureRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterPrecurvatureRange"))
writedlm(string("test/newTestDataGenerators/SMWFlutterPrecurvatureRange/flutterSpeed.txt"), test_data_matrix(flutterSpeed))

# Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with root angle as the varying parameter
include("../examples/eigenAeroelastic/SMWFlutterPrecurvatureRange2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterPrecurvatureRange2"))
writedlm(string("test/newTestDataGenerators/SMWFlutterPrecurvatureRange2/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/SMWFlutterPrecurvatureRange2/damps.txt"), test_data_matrix(dampsOsc))

# Flutter boundary analysis of the sixteen-meter-wing as a function of the bending-torsion coupling factor
include("../examples/eigenAeroelastic/SMWFlutterStructuralCouplingRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterStructuralCouplingRange"))
writedlm("test/newTestDataGenerators/SMWFlutterStructuralCouplingRange/flutterSpeed.txt", test_data_matrix(flutterSpeed))

# Flutter boundary analysis of the sixteen-meter-wing as a function of the tip displacement
include("../examples/eigenAeroelastic/SMWFlutterTipDispRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterTipDispRange"))
writedlm("test/newTestDataGenerators/SMWFlutterTipDispRange/flutterSpeed.txt", test_data_matrix(flutterSpeed))
writedlm("test/newTestDataGenerators/SMWFlutterTipDispRange/flutterFreq.txt", test_data_matrix(flutterFreq))

# Linear flutter analysis of the sixteen-meter-wing
include("../examples/eigenAeroelastic/SMWLinearFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWLinearFlutter"))
writedlm("test/newTestDataGenerators/SMWLinearFlutter/flutterSpeed.txt", test_data_matrix(flutterSpeed))
writedlm("test/newTestDataGenerators/SMWLinearFlutter/flutterFreq.txt", test_data_matrix(flutterFreq))

# Torsional divergence analysis of a straight wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/straightWingTorsionalDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/straightWingTorsionalDivergence"))
writedlm("test/newTestDataGenerators/straightWingTorsionalDivergence/damps.txt", test_data_matrix(dampsNonOsc))

# Coupled bending-torsion divergence analysis of a straight wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/straightWingBendingTorsionalDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/straightWingBendingTorsionalDivergence"))
writedlm("test/newTestDataGenerators/straightWingBendingTorsionalDivergence/dampsNonOsc.txt", test_data_matrix(dampsNonOsc))
writedlm("test/newTestDataGenerators/straightWingBendingTorsionalDivergence/UD.txt", test_data_matrix(UD))

# Bending divergence analysis of a swept wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/sweptWingBendingDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptWingBendingDivergence"))
writedlm("test/newTestDataGenerators/sweptWingBendingDivergence/dampsNonOsc.txt", test_data_matrix(dampsNonOsc))
writedlm("test/newTestDataGenerators/sweptWingBendingDivergence/UD.txt", test_data_matrix(UD))

# Coupled bending-torsion divergence analysis of a swept wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/sweptWingBendingTorsionalDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptWingBendingTorsionalDivergence"))
writedlm("test/newTestDataGenerators/sweptWingBendingTorsionalDivergence/dampsNonOsc.txt", test_data_matrix(dampsNonOsc))
writedlm("test/newTestDataGenerators/sweptWingBendingTorsionalDivergence/UD.txt", test_data_matrix(UD))

# Torsional divergence analysis of a swept wing
Random.seed!(1234)
include("../examples/eigenAeroelastic/sweptWingTorsionalDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptWingTorsionalDivergence"))
writedlm("test/newTestDataGenerators/sweptWingTorsionalDivergence/dampsNonOsc.txt", test_data_matrix(dampsNonOsc))
writedlm("test/newTestDataGenerators/sweptWingTorsionalDivergence/UD.txt", test_data_matrix(UD))

# Flutter analysis of the swept Pazy wing
include("../examples/eigenAeroelastic/sweptPazyFlutterPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptPazyFlutterPitchRange"))
writedlm(string("test/newTestDataGenerators/sweptPazyFlutterPitchRange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/sweptPazyFlutterPitchRange/damps.txt"), test_data_matrix(dampsOsc))

# Flutter analysis of the undeformed swept Pazy wing
include("../examples/eigenAeroelastic/sweptPazyUndeformedFlutterSweepRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptPazyUndeformedFlutterSweepRange"))
writedlm(string("test/newTestDataGenerators/sweptPazyUndeformedFlutterSweepRange/freqs.txt"), test_data_matrix(freqsOsc))
writedlm(string("test/newTestDataGenerators/sweptPazyUndeformedFlutterSweepRange/damps.txt"), test_data_matrix(dampsOsc))

# Bending divergence analyses of swept wings
include("../examples/eigenAeroelastic/sweptWingBendingDivergenceSweepRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptWingBendingDivergenceSweepRange"))
writedlm(string("test/newTestDataGenerators/sweptWingBendingDivergenceSweepRange/UD.txt"), test_data_matrix(UD))

# Coupled bending-torsion divergence analyses of swept wings
include("../examples/eigenAeroelastic/sweptWingBendingTorsionalDivergenceSweepRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptWingBendingTorsionalDivergenceSweepRange"))
writedlm(string("test/newTestDataGenerators/sweptWingBendingTorsionalDivergenceSweepRange/UD.txt"), test_data_matrix(UD))

# Torsional divergence analyses of swept wings
include("../examples/eigenAeroelastic/sweptWingTorsionalDivergenceSweepRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptWingTorsionalDivergenceSweepRange"))
writedlm(string("test/newTestDataGenerators/sweptWingTorsionalDivergenceSweepRange/UD.txt"), test_data_matrix(UD))

# Eigen aeroelastic analysis of the Tang&Dowell wing at varying airspeed
include("../examples/eigenAeroelastic/TDWingAirspeedRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/TDWingAirspeedRange"))
writedlm("test/newTestDataGenerators/TDWingAirspeedRange/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/TDWingAirspeedRange/tip_twist.txt", test_data_matrix(tip_twist))

# Pitch divergence analysis of a typical section
Random.seed!(1234)
include("../examples/eigenAeroelastic/typicalSectionPitchDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/typicalSectionPitchDivergence"))
writedlm("test/newTestDataGenerators/typicalSectionPitchDivergence/damps.txt", test_data_matrix(dampsNonOsc))

# Flutter and divergence analysis of a typical section
include("../examples/eigenAeroelastic/typicalSectionFlutterAndDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/typicalSectionFlutterAndDivergence"))
writedlm("test/newTestDataGenerators/typicalSectionFlutterAndDivergence/freqsOsc.txt", test_data_matrix(freqsOsc))
writedlm("test/newTestDataGenerators/typicalSectionFlutterAndDivergence/dampsOsc.txt", test_data_matrix(dampsOsc))
writedlm("test/newTestDataGenerators/typicalSectionFlutterAndDivergence/freqsNonOsc.txt", test_data_matrix(freqsNonOsc))
writedlm("test/newTestDataGenerators/typicalSectionFlutterAndDivergence/dampsNonOsc.txt", test_data_matrix(dampsNonOsc))
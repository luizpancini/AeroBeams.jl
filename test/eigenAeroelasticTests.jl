# Seed the random initial vector of the Krylov eigensolver for reproducible modes
using Random

!is_ci && @testset "Flutter analysis of the Blended-Wing-Body flying wing" begin
    include_example("examples/eigenAeroelastic/BWBflutter.jl")
    @test matches_test_data(freqsOsc, "BWBflutter", "freqs.txt")
    @test matches_test_data(dampsOsc, "BWBflutter", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the matched wing of the Blended-Wing-Body flying wing" begin
    include_example("examples/eigenAeroelastic/BWBmatchedWingFlutter.jl")
    @test matches_test_data(freqsOsc, "BWBmatchedWingFlutter", "freqsOsc.txt")
    @test matches_test_data(dampsOsc, "BWBmatchedWingFlutter", "dampsOsc.txt")
    @test matches_test_data(trimAoA, "BWBmatchedWingFlutter", "trimAoA.txt")
    @test matches_test_data(trimThrust, "BWBmatchedWingFlutter", "trimThrust.txt")
    @test matches_test_data(trimδ, "BWBmatchedWingFlutter", "trimDelta.txt")
    @test matches_test_data(flutterSpeeds, "BWBmatchedWingFlutter", "flutterSpeeds.txt")
end

!is_ci && @testset "Flutter analysis of the conventional HALE aircraft in free flight with structural stiffness as the varying parameter" begin
    include_example("examples/eigenAeroelastic/cHALELambdaRange.jl")
    @test matches_test_data(freqsOsc, "cHALELambdaRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "cHALELambdaRange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the conventional HALE aircraft in free flight with airspeed and structural stiffness as the varying parameters" begin
    include_example("examples/eigenAeroelastic/cHALELURange.jl")
    @test matches_test_data(freqsOsc, "cHALELURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "cHALELURange", "damps.txt")
end

@testset "Flutter analysis of the conventional HALE aircraft in free flight with airspeed as the varying parameter" begin
    include_example("examples/eigenAeroelastic/cHALEURange.jl")
    @test matches_test_data(freqsOsc, "cHALEURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "cHALEURange", "damps.txt")
end

@testset "Flutter analysis of the Goland wing" begin
    include_example("examples/eigenAeroelastic/GolandWingFlutter.jl")
    @test matches_test_data(freqsOsc, "GolandWingFlutter", "freqs.txt")
    @test matches_test_data(dampsOsc, "GolandWingFlutter", "damps.txt")
end

@testset "Flutter analysis of the baseline Healy free FFWT wing with tip loss, root pitch and airspeed as the varying parameters" begin
    include_example("examples/eigenAeroelastic/HealyBaselineFFWTfreeFlutterAoARangeURange.jl")
    @test matches_test_data(freqsOsc, "HealyBaselineFFWTfreeFlutterAoARangeURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "HealyBaselineFFWTfreeFlutterAoARangeURange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the baseline Healy free FFWT wing with flare angle and airspeed as the varying parameters" begin
    include_example("examples/eigenAeroelastic/HealyBaselineFFWTfreeFlutterFlareRangeURange.jl")
    @test matches_test_data(freqsOsc, "HealyBaselineFFWTfreeFlutterFlareRangeURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "HealyBaselineFFWTfreeFlutterFlareRangeURange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the baseline Healy locked FFWT wing with airspeed as the varying parameters" begin
    include_example("examples/eigenAeroelastic/HealyBaselineFFWTlockedFlutterURange.jl")
    @test matches_test_data(freqsOsc, "HealyBaselineFFWTlockedFlutterURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "HealyBaselineFFWTlockedFlutterURange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Healy LCO FFWT wing with root pitch angle as the varying parameter" begin
    include_example("examples/eigenAeroelastic/HealyLCOFFWTfreeFlutterAoARange.jl")
    @test matches_test_data(freqsOsc, "HealyLCOFFWTfreeFlutterAoARange", "freqs.txt")
    @test matches_test_data(dampsOsc, "HealyLCOFFWTfreeFlutterAoARange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Healy FFWT wing with sideslip angle as the varying parameter" begin
    include_example("examples/eigenAeroelastic/HealySideslipFFWTflutterSideslipRange.jl")
    @test matches_test_data(freqsOsc, "HealySideslipFFWTflutterSideslipRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "HealySideslipFFWTflutterSideslipRange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Healy sideslip FFWT wing with airspeed as the varying parameter" begin
    include_example("examples/eigenAeroelastic/HealySideslipFFWTflutterURange.jl")
    @test matches_test_data(freqsOsc, "HealySideslipFFWTflutterURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "HealySideslipFFWTflutterURange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Helios flying-wing in free flight with payload and structural stiffness as the varying parameters" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/heliosFlutterPLambdaRange.jl")
    @test matches_test_data(freqsOsc, "heliosFlutterPLambdaRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "heliosFlutterPLambdaRange", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Helios flying-wing in free flight with payload as the varying parameter" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/heliosFlutterPRange.jl")
    @test matches_test_data(freqsOsc, "heliosFlutterPRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "heliosFlutterPRange", "damps.txt")
end

@testset "Flutter analysis of the Helios flying-wing in free flight with airspeed as the varying parameter" begin
    include_example("examples/eigenAeroelastic/heliosFlutterURange.jl")
    @test matches_test_data(freqsOsc, "heliosFlutterURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "heliosFlutterURange", "damps.txt")
end

@testset "Flutter analysis of the wing of the Helios flying-wing" begin
    include_example("examples/eigenAeroelastic/heliosWingFlutter.jl")
    @test matches_test_data(freqsOsc, "heliosWingFlutter", "freqs.txt")
    @test matches_test_data(dampsOsc, "heliosWingFlutter", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the Pazy FFWT wing with airspeed as the varying parameter" begin
    include_example("examples/eigenAeroelastic/PazyFFWTflutterURange.jl")
    @test matches_test_data(freqsOsc, "PazyFFWTflutterURange", "freqs.txt")
    @test matches_test_data(dampsOsc, "PazyFFWTflutterURange", "damps.txt")
    @test matches_test_data(ϕHinge, "PazyFFWTflutterURange", "phiHinge.txt")
end

@testset "Flutter analysis of the Pazy wing" begin
    include_example("examples/eigenAeroelastic/PazyWingFlutter.jl")
    @test matches_test_data(flutterOnsetSpeeds, "PazyWingFlutter", "flutterOnsetSpeeds.txt")
    @test matches_test_data(flutterOnsetFreqs, "PazyWingFlutter", "flutterOnsetFreqs.txt")
    @test matches_test_data(flutterOnsetDisp, "PazyWingFlutter", "flutterOnsetDisp.txt")
    @test matches_test_data(flutterOffsetSpeeds, "PazyWingFlutter", "flutterOffsetSpeeds.txt")
    @test matches_test_data(flutterOffsetFreqs, "PazyWingFlutter", "flutterOffsetFreqs.txt")
    @test matches_test_data(flutterOffsetDisp, "PazyWingFlutter", "flutterOffsetDisp.txt")
end

!is_ci && @testset "Flutter and divergence analysis of the Pazy wing" begin
    include_example("examples/eigenAeroelastic/PazyWingFlutterAndDivergence.jl")
    @test matches_test_data(flutterOnsetSpeeds, "PazyWingFlutterAndDivergence", "flutterOnsetSpeeds.txt")
    @test matches_test_data(flutterOnsetFreqs, "PazyWingFlutterAndDivergence", "flutterOnsetFreqs.txt")
end

!is_ci && @testset "Flutter analysis of the Pazy wing with varying root pitch angle" begin
    include_example("examples/eigenAeroelastic/PazyWingFlutterPitchRange.jl")
    @test matches_test_data(freqsOsc, "PazyWingFlutterPitchRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "PazyWingFlutterPitchRange", "damps.txt")
end

@testset "Flutter analysis of the Pazy wing with varying tip mass positions" begin
    include_example("examples/eigenAeroelastic/PazyWingFlutterTipMassRange.jl")
    @test matches_test_data(freqsOsc, "PazyWingFlutterTipMassRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "PazyWingFlutterTipMassRange", "damps.txt")
end

@testset "Frequency analysis of the Pazy wing" begin
    include_example("examples/eigenAeroelastic/PazyWingFreqsEvolution.jl")
    @test matches_test_data(freqsOsc, "PazyWingFreqsEvolution", "freqs.txt")
end

!is_ci && @testset "Flutter analysis of the sixteen-meter-wing" begin
    include_example("examples/eigenAeroelastic/SMWFlutter.jl")
    @test matches_test_data(x1_def[end]/L, "SMWFlutter", "x1_def.txt")
    @test matches_test_data(x3_def[end]/L, "SMWFlutter", "x3_def.txt")
    @test matches_test_data(α_of_x1[end], "SMWFlutter", "alpha.txt")
    @test matches_test_data(freqsOsc, "SMWFlutter", "freqs.txt")
    @test matches_test_data(dampsOsc, "SMWFlutter", "damps.txt")
end

!is_ci && @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the pitch angle" begin
    include_example("examples/eigenAeroelastic/SMWFlutterPitchRange.jl")
    @test matches_test_data(freqsOsc, "SMWFlutterPitchRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "SMWFlutterPitchRange", "damps.txt")
end
    
!is_ci && @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with tip load as the varying parameter" begin
    include_example("examples/eigenAeroelastic/SMWFlutterPrecurvatureRange.jl")
    @test matches_test_data(flutterSpeed, "SMWFlutterPrecurvatureRange", "flutterSpeed.txt")
end

!is_ci && @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with root angle as the varying parameter" begin
    include_example("examples/eigenAeroelastic/SMWFlutterPrecurvatureRange2.jl")
    @test matches_test_data(freqsOsc, "SMWFlutterPrecurvatureRange2", "freqs.txt")
    @test matches_test_data(dampsOsc, "SMWFlutterPrecurvatureRange2", "damps.txt")
end

!is_ci && @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the bending-torsion coupling factor" begin
    include_example("examples/eigenAeroelastic/SMWFlutterStructuralCouplingRange.jl")
    @test matches_test_data(flutterSpeed, "SMWFlutterStructuralCouplingRange", "flutterSpeed.txt")
end

!is_ci && @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the tip displacement" begin
    include_example("examples/eigenAeroelastic/SMWFlutterTipDispRange.jl")
    @test matches_test_data(flutterSpeed, "SMWFlutterTipDispRange", "flutterSpeed.txt")
    @test matches_test_data(flutterFreq, "SMWFlutterTipDispRange", "flutterFreq.txt")
end

!is_ci && @testset "Linear flutter analysis of the sixteen-meter-wing" begin
    include_example("examples/eigenAeroelastic/SMWLinearFlutter.jl")
    @test matches_test_data(flutterSpeed, "SMWLinearFlutter", "flutterSpeed.txt")
    @test matches_test_data(flutterFreq, "SMWLinearFlutter", "flutterFreq.txt")
end

!is_ci && @testset "Coupled bending-torsion divergence analysis of a straight wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/straightWingBendingTorsionalDivergence.jl")
    @test matches_test_data(dampsNonOsc, "straightWingBendingTorsionalDivergence", "dampsNonOsc.txt")
    @test matches_test_data(UD, "straightWingBendingTorsionalDivergence", "UD.txt")
end

!is_ci && @testset "Torsional divergence analysis of a straight wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/straightWingTorsionalDivergence.jl")
    @test matches_test_data(dampsNonOsc, "straightWingTorsionalDivergence", "damps.txt")
end

!is_ci && @testset "Flutter analysis of the swept Pazy wing" begin
    include_example("examples/eigenAeroelastic/sweptPazyFlutterPitchRange.jl")
    @test matches_test_data(freqsOsc, "sweptPazyFlutterPitchRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "sweptPazyFlutterPitchRange", "damps.txt")
end

@testset "Flutter analysis of the undeformed swept Pazy wing" begin
    include_example("examples/eigenAeroelastic/sweptPazyUndeformedFlutterSweepRange.jl")
    @test matches_test_data(freqsOsc, "sweptPazyUndeformedFlutterSweepRange", "freqs.txt")
    @test matches_test_data(dampsOsc, "sweptPazyUndeformedFlutterSweepRange", "damps.txt")
end

@testset "Bending divergence analyses of swept wing" begin
    include_example("examples/eigenAeroelastic/sweptWingBendingDivergenceSweepRange.jl")
    @test matches_test_data(UD, "sweptWingBendingDivergenceSweepRange", "UD.txt")
end

@testset "Coupled bending-torsion divergence analyses of swept wings" begin
    include_example("examples/eigenAeroelastic/sweptWingBendingTorsionalDivergenceSweepRange.jl")
    @test matches_test_data(UD, "sweptWingBendingTorsionalDivergenceSweepRange", "UD.txt")
end

@testset "Torsional divergence analyses of swept wings" begin
    include_example("examples/eigenAeroelastic/sweptWingTorsionalDivergenceSweepRange.jl")
    @test matches_test_data(UD, "sweptWingTorsionalDivergenceSweepRange", "UD.txt")
end

@testset "Bending divergence analysis of a swept wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/sweptWingBendingDivergence.jl")
    @test matches_test_data(dampsNonOsc, "sweptWingBendingDivergence", "dampsNonOsc.txt")
    @test matches_test_data(UD, "sweptWingBendingDivergence", "UD.txt")
end

@testset "Coupled bending-torsion divergence analysis of a swept wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/sweptWingBendingTorsionalDivergence.jl")
    @test matches_test_data(dampsNonOsc, "sweptWingBendingTorsionalDivergence", "dampsNonOsc.txt")
    @test matches_test_data(UD, "sweptWingBendingTorsionalDivergence", "UD.txt")
end

@testset "Torsional divergence analysis of a swept wing" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/sweptWingTorsionalDivergence.jl")
    @test matches_test_data(dampsNonOsc, "sweptWingTorsionalDivergence", "dampsNonOsc.txt")
    @test matches_test_data(UD, "sweptWingTorsionalDivergence", "UD.txt")
end

@testset "Steady aeroelastic analysis of the Tang&Dowell wing at varying airspeed" begin
    include_example("examples/eigenAeroelastic/TDWingAirspeedRange.jl")
    @test matches_test_data(tip_u3, "TDWingAirspeedRange", "tip_u3.txt")
    @test matches_test_data(tip_twist, "TDWingAirspeedRange", "tip_twist.txt")
end

@testset "Flutter and divergence analysis of a typical section" begin
    include_example("examples/eigenAeroelastic/typicalSectionFlutterAndDivergence.jl")
    @test matches_test_data(freqsOsc, "typicalSectionFlutterAndDivergence", "freqsOsc.txt")
    @test matches_test_data(dampsOsc, "typicalSectionFlutterAndDivergence", "dampsOsc.txt")
    @test matches_test_data(freqsNonOsc, "typicalSectionFlutterAndDivergence", "freqsNonOsc.txt")
    @test matches_test_data(dampsNonOsc, "typicalSectionFlutterAndDivergence", "dampsNonOsc.txt")
end

@testset "Torsional divergence analysis of a typical section" begin
    Random.seed!(1234)
    include_example("examples/eigenAeroelastic/typicalSectionPitchDivergence.jl")
    @test matches_test_data(dampsNonOsc, "typicalSectionPitchDivergence", "damps.txt")
end
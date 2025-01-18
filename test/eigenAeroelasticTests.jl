# Eigen aeroelastic problems

@testset "Flutter analysis of the conventional HALE aircraft in free flight with structural stiffness as the varying parameter" begin
    include("examples/conventionalHALELambdaRange.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALELambdaRange", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALELambdaRange", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Flutter analysis of the conventional HALE aircraft in free flight with airspeed as the varying parameter" begin
    include("examples/conventionalHALEURange.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEURange", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEURange", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Flutter analysis of the Goland wing" begin
    include("examples/GolandWingFlutter.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "GolandWingFlutter", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "GolandWingFlutter", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Flutter analysis of the baseline Healy free FFWT wing with tip loss, root pitch and airspeed as the varying parameters" begin
    include("examples/HealyBaselineFFWTfreeFlutterAoARangeURange.jl")
    # Self-comparison
    for (i,τ) in enumerate(τRange)
        for (j,θ) in enumerate(θRange)
            freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyBaselineFFWTfreeFlutterAoARangeURange", string("freqs",i,j,".txt")))
            damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyBaselineFFWTfreeFlutterAoARangeURange", string("damps",i,j,".txt")))
            @test isapprox(hcat(freqs[i,j,:]...)', freqs_, atol=SELFatol, nans=true)
            @test isapprox(hcat(damps[i,j,:]...)', damps_, atol=SELFatol, nans=true)
        end
    end
end

@testset "Flutter analysis of the baseline Healy free FFWT wing with flare angle and airspeed as the varying parameters" begin
    include("examples/HealyBaselineFFWTfreeFlutterFlareRangeURange.jl")
    # Self-comparison
    for (i,Λ) in enumerate(ΛRange)
        freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyBaselineFFWTfreeFlutterFlareRangeURange", string("freqs",i,".txt")))
        damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyBaselineFFWTfreeFlutterFlareRangeURange", string("damps",i,".txt")))
        @test hcat(freqs[i,:]...)' ≈ freqs_ atol=SELFatol
        @test hcat(damps[i,:]...)' ≈ damps_ atol=SELFatol
    end
end

@testset "Flutter analysis of the baseline Healy locked FFWT wing with airspeed as the varying parameters" begin
    include("examples/HealyBaselineFFWTlockedFlutterURange.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyBaselineFFWTlockedFlutterURange", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyBaselineFFWTlockedFlutterURange", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Flutter analysis of the Healy FFWT wing with sideslip angle as the varying parameter" begin
    include("examples/HealyFFWTflutterSideslipRange.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyFFWTflutterSideslipRange", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyFFWTflutterSideslipRange", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Flutter analysis of the Healy FFWT wing with airspeed as the varying parameter" begin
    include("examples/HealyFFWTflutterURange.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyFFWTflutterURange", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyFFWTflutterURange", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Flutter analysis of the Helios flying-wing in free flight with airspeed as the varying parameter" begin
    include("examples/heliosFlutterURange.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosFlutterURange", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosFlutterURange", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=1e-2
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Flutter analysis of the wing of the Helios flying-wing" begin
    include("examples/heliosWingFlutter.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosWingFlutter", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosWingFlutter", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Flutter analysis of the Pazy FFWT wing with airspeed as the varying parameter" begin
    include("examples/PazyFFWTflutterURange.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTflutterURange", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTflutterURange", "damps.txt"))
    ϕHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTflutterURange", "phiHinge.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
    @test ϕHinge ≈ ϕHinge_ atol=SELFatol
end

@testset "Flutter analysis of the Pazy wing" begin
    include("examples/PazyWingFlutter.jl")
    # Self-comparison
    flutterOnsetSpeedsOfMode_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingFlutter", "flutterOnsetSpeedsOfMode.txt"))
    flutterOnsetFreqsOfMode_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingFlutter", "flutterOnsetFreqsOfMode.txt"))
    flutterOnsetDispOfMode_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingFlutter", "flutterOnsetDispOfMode.txt"))
    flutterOffsetSpeedsOfMode_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingFlutter", "flutterOffsetSpeedsOfMode.txt"))
    flutterOffsetFreqsOfMode_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingFlutter", "flutterOffsetFreqsOfMode.txt"))
    flutterOffsetDispOfMode_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingFlutter", "flutterOffsetDispOfMode.txt"))
    @test hcat(filter(x->!isempty(x),flutterOnsetSpeedsOfMode)...)' ≈ flutterOnsetSpeedsOfMode_ atol=SELFatol
    @test hcat(filter(x->!isempty(x),flutterOnsetFreqsOfMode)...)' ≈ flutterOnsetFreqsOfMode_ atol=SELFatol
    @test hcat(filter(x->!isempty(x),flutterOnsetDispOfMode)...)' ≈ flutterOnsetDispOfMode_ atol=SELFatol
    @test hcat(filter(x->!isempty(x),flutterOffsetSpeedsOfMode)...)' ≈ flutterOffsetSpeedsOfMode_ atol=SELFatol
    @test hcat(filter(x->!isempty(x),flutterOffsetFreqsOfMode)...)' ≈ flutterOffsetFreqsOfMode_ atol=SELFatol
    @test hcat(filter(x->!isempty(x),flutterOffsetDispOfMode)...)' ≈ flutterOffsetDispOfMode_ atol=SELFatol
end

@testset "Flutter analysis of the Pazy wing with varying tip mass positions" begin
    include("examples/PazyWingFlutterTipMassRange.jl")
    # Self-comparison
    for i=1:3
        freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingFlutterTipMassRange", string("freqs",i,".txt")))
        damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingFlutterTipMassRange", string("damps",i,".txt")))
        @test hcat(freqs[i,:]...)' ≈ freqs_ atol=SELFatol
        @test hcat(damps[i,:]...)' ≈ damps_ atol=SELFatol
    end
end

@testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the bending-torsion coupling factor" begin
    include("examples/SMWFlutterStructuralCouplingRange.jl")
    # Reference comparison (flutter at zero displacement)
    @test flutterSpeed[1,1] ≈ flutterSpeedVsTipLoadΨm02[2,1] rtol=2e-2
    @test flutterSpeed[2,1] ≈ flutterSpeedVsTipLoadΨ0[2,1] rtol=2e-2
    @test flutterSpeed[3,1] ≈ flutterSpeedVsTipLoadΨp02[2,1] rtol=2e-2
    # Self-comparison
    flutterSpeed_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutterStructuralCouplingRange", "flutterSpeed.txt"))
    @test flutterSpeed ≈ flutterSpeed_ atol=SELFatol
end

@testset "Linear flutter analysis of the sixteen-meter-wing" begin
    include("examples/SMWLinearFlutter.jl")
    # Reference comparison
    @test flutterSpeed ≈ flutterSpeedRef rtol=2e-2
    @test flutterFreq ≈ flutterFreqRef rtol=2e-2
    # Self-comparison
    flutterSpeed_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWLinearFlutter", "flutterSpeed.txt"))[1]
    flutterFreq_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWLinearFlutter", "flutterFreq.txt"))[1]
    @test flutterSpeed ≈ flutterSpeed_ atol=SELFatol
    @test flutterFreq ≈ flutterFreq_ atol=SELFatol
end

@testset "Flutter and divergence analysis of a typical section" begin
    include("examples/typicalSectionFlutterAndDivergence.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "typicalSectionFlutterAndDivergence", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "typicalSectionFlutterAndDivergence", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

# Run tests that currently fail in CI, if applicable
if runCIfailingTests && (get(ENV, "CI", false) || get(ENV, "GITHUB_ACTIONS", false))

    @testset "Flutter analysis of the Blended-Wing-Body flying wing" begin
        include("examples/BWBflutter.jl")
        # Self-comparison
        freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "BWBflutter", "freqs.txt"))
        damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "BWBflutter", "damps.txt"))
        @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
        @test hcat(damps...)' ≈ damps_ atol=SELFatol
    end

    @testset "Flutter analysis of the conventional HALE aircraft in free flight with airspeed and structural stiffness as the varying parameters" begin
        include("examples/conventionalHALELURange.jl")
        # Self-comparison
        for (i,λ) in enumerate(λRange)
            freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALELURange", string("freqs",i,".txt")))
            damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALELURange", string("damps",i,".txt")))
            @test hcat(freqs[i,:]...)' ≈ freqs_ atol=SELFatol
            @test hcat(damps[i,:]...)' ≈ damps_ atol=SELFatol
        end
    end

    @testset "Flutter analysis of the Helios flying-wing in free flight with payload and structural stiffness as the varying parameters" begin
        include("examples/heliosFlutterPLambdaRange.jl")
        # Self-comparison
        for (i,λ) in enumerate(λRange)
            freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosFlutterPLambdaRange", string("freqs",i,".txt")))
            damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosFlutterPLambdaRange", string("damps",i,".txt")))
            @test hcat(freqs[i,:]...)' ≈ freqs_ atol=SELFatol
            @test hcat(damps[i,:]...)' ≈ damps_ atol=SELFatol
        end
    end
    
    @testset "Flutter analysis of the Helios flying-wing in free flight with payload as the varying parameter" begin
        include("examples/heliosFlutterPRange.jl")
        # Self-comparison
        freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosFlutterPRange", "freqs.txt"))
        damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosFlutterPRange", "damps.txt"))
        @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
        @test hcat(damps...)' ≈ damps_ atol=SELFatol
    end

    @testset "Flutter and divergence analysis of the Pazy wing" begin
        include("examples/PazyWingFlutterAndDivergence.jl")
        # Self-comparison
        flutterOnsetSpeedsOfMode_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingFlutterAndDivergence", "flutterOnsetSpeedsOfMode.txt"))
        flutterOnsetFreqsOfMode_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingFlutterAndDivergence", "flutterOnsetFreqsOfMode.txt"))
        @test hcat(filter(x->!isempty(x),flutterOnsetSpeedsOfMode)...)' ≈ flutterOnsetSpeedsOfMode_ atol=SELFatol
        @test hcat(filter(x->!isempty(x),flutterOnsetFreqsOfMode)...)' ≈ flutterOnsetFreqsOfMode_ atol=SELFatol
    end
    
    @testset "Flutter analysis of the Pazy wing with varying root pitch angle" begin
        include("examples/PazyWingFlutterPitchRange.jl")
        # Self-comparison
        for (i,θ) in enumerate(θRange)
            freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingFlutterPitchRange", string("freqs",i,".txt")))
            damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingFlutterPitchRange", string("damps",i,".txt")))
            @test hcat(freqs[i,:]...)' ≈ freqs_ atol=SELFatol
            @test hcat(damps[i,:]...)' ≈ damps_ atol=SELFatol
        end
    end

    @testset "Flutter analysis of the sixteen-meter-wing" begin
        include("examples/SMWFlutter.jl")
         # Self-comparison
        x1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutter", "x1_def.txt"))
        x3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutter", "x3_def.txt"))
        α_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutter", "alpha.txt"))
        @test x1_def[end]/L ≈ x1_ atol=SELFatol
        @test x3_def[end]/L ≈ x3_ atol=SELFatol
        @test α_of_x1[end]*180/pi ≈ α_ atol=SELFatol
        for mode in 1:nModes
            freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutter", string("freqsMode",mode,".txt")))
            damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutter", string("dampsMode",mode,".txt")))
            @test modeFrequencies[mode] ≈ freqs_ atol=SELFatol
            @test modeDampings[mode] ≈ damps_ atol=SELFatol
        end
    end
    
    @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the pitch angle" begin
        include("examples/SMWFlutterPitchRange.jl")
        # Reference comparison
        @test freqs[21,151][2] ≈ freqVsSpeedRootAoA2[2,end] rtol = 0.1
        @test freqs[31,151][2] ≈ freqVsSpeedRootAoA3[2,end] rtol = 0.1
        @test freqs[51,151][2] ≈ freqVsSpeedRootAoA5[2,end] rtol = 0.1
        @test damps[21,151][2] ≈ dampVsSpeedRootAoA2[2,end] atol = 0.1
        @test damps[31,151][2] ≈ dampVsSpeedRootAoA3[2,end] atol = 0.1
        @test damps[51,151][2] ≈ dampVsSpeedRootAoA5[2,end] atol = 0.1
        # Self-comparison
        for i in 1:nModes
            freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutterPitchRange", string("freqsMode",i,".txt")))
            damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutterPitchRange", string("dampsMode",i,".txt")))
            @test modeFrequencies[i] ≈ freqs_ atol=SELFatol
            @test modeDampings[i] ≈ damps_ atol=SELFatol
        end
    end
    
    @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with tip load as the varying parameter" begin
        include("examples/SMWFlutterPrecurvatureRange.jl")
        # Self-comparison
        for i in eachindex(kRange)
            flutterSpeed_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutterPrecurvatureRange", string("flutterSpeedk",i,".txt")))
            @test flutterSpeed[i,:] ≈ flutterSpeed_ atol=SELFatol
        end
    end

    @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with root angle as the varying parameter" begin
        include("examples/SMWFlutterPrecurvatureRange2.jl")
        # Self-comparison
        for (ki,k) in enumerate(kRange)
            for (i,θ) in enumerate(θRange)
                freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutterPrecurvatureRange2", string("freqs_k",ki,"th",i,".txt")))
                damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutterPrecurvatureRange2", string("damps_k",ki,"th",i,".txt")))
                @test hcat(freqs[ki,i,:]...)' ≈ freqs_ atol=SELFatol
                @test hcat(damps[ki,i,:]...)' ≈ damps_ atol=SELFatol
            end
        end
    end

    @testset "Flutter boundary analysis of the sixteen-meter-wing as a function of the tip displacement" begin
        include("examples/SMWFlutterTipDispRange.jl")
        # Reference comparison (flutter at zero displacement)
        iZeroDisp = div(length(F3Range),2)+1
        @test flutterSpeed[iZeroDisp] ≈ flutterSpeedRef[2,1] rtol=2e-2
        @test flutterFreq[iZeroDisp] ≈ flutterFreqRef[2,1] rtol=3e-2
        # Self-comparison
        flutterSpeed_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutterTipDispRange", "flutterSpeed.txt"))
        flutterFreq_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWFlutterTipDispRange", "flutterFreq.txt"))
        @test flutterSpeed ≈ flutterSpeed_ atol=SELFatol
        @test flutterFreq ≈ flutterFreq_ atol=SELFatol
    end

end
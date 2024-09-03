# Steady aeroelastic problems

@testset "Steady analysis of the Pazy wing with flared folding wing tip (FFWT)" begin
    include("examples/PazyFFWTsteady.jl")
    # Self-comparison
    u3_of_x1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteady", "u3_of_x1.txt"))
    p2_of_x1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteady", "p2_of_x1.txt"))
    M2_of_x1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteady", "M2_of_x1.txt"))
    @test u3_of_x1 ≈ u3_of_x1_ atol=SELFatol
    @test p2_of_x1 ≈ p2_of_x1_ atol=SELFatol
    @test M2_of_x1 ≈ M2_of_x1_ atol=SELFatol
end

@testset "Steady analysis of the Pazy wing with varying root pitch angle" begin
    include("examples/PazyWingPitchRange.jl")
    # Self-comparison
    tip_AoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingPitchRange", "tip_AoA.txt"))
    tip_OOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingPitchRange", "tip_OOP.txt"))
    tip_IP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingPitchRange", "tip_IP.txt"))
    tip_twist_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingPitchRange", "tip_twist.txt"))
    # @test filter(!isnan,tip_AoA) ≈ filter(!isnan,tip_AoA_) atol=SELFatol
    @test tip_OOP ≈ tip_OOP_ atol=SELFatol
    @test tip_IP ≈ tip_IP_ atol=SELFatol
    @test tip_twist ≈ tip_twist_ atol=SELFatol
end

@testset "Steady aeroelastic analysis of the sixteen-meter-wing" begin
    include("examples/SMWSteady.jl")
    # Self-comparison
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWSteady", "tip_u3.txt"))
    tip_twist_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWSteady", "tip_twist.txt"))
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
    @test tip_twist ≈ tip_twist_ atol=SELFatol
end

@testset "Steady aeroelastic analysis of the Tang&Dowell wing at varying airspeed" begin
    include("examples/TDWingAirspeedRange.jl")
    # Reference comparison (chordwise bending and torsion aeroelastic frequencies at last airspeed)
    @test freqs[1,end][[2,4]] ≈ freqs_ref[[3,4],end] rtol=5e-2
    # Self-comparison
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "TDWingAirspeedRange", "tip_u3.txt"))
    tip_twist_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "TDWingAirspeedRange", "tip_twist.txt"))
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
    @test tip_twist ≈ tip_twist_ atol=SELFatol
end
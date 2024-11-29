# Steady aeroelastic problems

@testset "Steady aeroelastic analysis of the clamped conventional HALE" begin
    include("examples/conventionalHALEclampedSteady.jl")
    # Self-comparison
    x1_def_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEclampedSteady", "x1_def.txt"))
    x3_def_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEclampedSteady", "x3_def.txt"))
    @test x1_def ≈ x1_def_ atol=SELFatol
    @test x3_def ≈ x3_def_ atol=SELFatol
end

@testset "Steady analysis of the Pazy wing with flared folding wing tip (FFWT)" begin
    include("examples/PazyFFWTsteady.jl")
    # Self-comparison
    u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteady", "u1.txt"))
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteady", "u3.txt"))
    p1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteady", "p1.txt"))
    p2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteady", "p2.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteady", "M2.txt"))
    α_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteady", "alpha.txt"))
    cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteady", "cn.txt"))
    hingeBalanceM_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteady", "hingeBalanceM.txt"))[1]
    @test u1 ≈ u1_ atol=SELFatol
    @test u3 ≈ u3_ atol=SELFatol
    @test p1 ≈ p1_ atol=SELFatol
    @test p2 ≈ p2_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
    @test α ≈ α_ atol=SELFatol
    @test cn ≈ cn_ atol=SELFatol
    @test hingeBalanceM ≈ hingeBalanceM_ atol=SELFatol
end

@testset "Steady analysis of the Pazy wing with flared folding wing tip (FFWT) and varying airspeed" begin
    include("examples/PazyFFWTsteadyURange.jl")
    # Self-comparison
    u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyURange", "u1.txt"))
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyURange", "u3.txt"))
    p1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyURange", "p1.txt"))
    p2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyURange", "p2.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyURange", "M2.txt"))
    α_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyURange", "alpha.txt"))
    cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyURange", "cn.txt"))
    @test vcat(u1'...) ≈ u1_ atol=SELFatol
    @test vcat(u3'...) ≈ u3_ atol=SELFatol
    @test vcat(p1'...) ≈ p1_ atol=SELFatol
    @test vcat(p2'...) ≈ p2_ atol=SELFatol
    @test vcat(M2'...) ≈ M2_ atol=SELFatol
    @test vcat(α'...) ≈ α_ atol=SELFatol
    @test vcat(cn'...) ≈ cn_ atol=SELFatol
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
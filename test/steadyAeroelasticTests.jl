# Steady aeroelastic problems

@testset "Steady aeroelastic analysis of the clamped conventional HALE" begin
    include("examples/conventionalHALEclampedSteady.jl")
    # Self-comparison
    x1_def_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEclampedSteady", "x1_def.txt"))
    x3_def_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEclampedSteady", "x3_def.txt"))
    @test x1_def ≈ x1_def_ atol=SELFatol
    @test x3_def ≈ x3_def_ atol=SELFatol
end

@testset "Steady analysis of the baseline Healy free FFWT wing with varying root pitch and airspeed" begin
    include("examples/HealyBaselineFFWTsteadyAoARangeURangeCoast.jl")
    # Self-comparison
    ϕHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyBaselineFFWTsteadyAoARangeURangeCoast", "phiHinge.txt"))
    u3Hinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyBaselineFFWTsteadyAoARangeURangeCoast", "u3Hinge.txt"))
    M2root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyBaselineFFWTsteadyAoARangeURangeCoast", "M2root.txt"))
    @test ϕHinge ≈ ϕHinge_ atol=SELFatol
    @test u3Hinge ≈ u3Hinge_ atol=SELFatol
    @test M2root ≈ M2root_ atol=SELFatol
end

@testset "Steady analysis of the Healy FFWT wing with varying flare angle and root pitch angle" begin
    include("examples/HealyFFWTsteadyFlareRangeAoARangeCoast.jl")
    # Self-comparison
    ϕHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyFFWTsteadyFlareRangeAoARangeCoast", "phiHinge.txt"))
    @test ϕHinge ≈ ϕHinge_ atol=SELFatol
end

@testset "Steady analysis of the Healy FFWT wing with varying flare angle, root pitch angle and airspeed" begin
    include("examples/HealyFFWTsteadyFlareRangeURangeAoARangeCoast.jl")
    # Self-comparison
    ϕHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyFFWTsteadyFlareRangeURangeAoARangeCoast", "phiHinge.txt"))
    @test hcat(ϕHinge...)' ≈ ϕHinge_ atol=SELFatol
end

@testset "Steady analysis of the Healy FFWT wing with varying wingtip twist, root pitch angle and sideslip angle" begin
    include("examples/HealyFFWTsteadyTwistRangeAoARangeSideslipRangeCoast.jl")
    # Self-comparison
    ϕHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyFFWTsteadyTwistRangeAoARangeSideslipRangeCoast", "phiHinge.txt"))
    @test isapprox(hcat(ϕHinge...)', ϕHinge_, atol=SELFatol, nans=true)
end

@testset "Steady analysis of the Pazy wing with a coasting FFWT" begin
    include("examples/PazyFFWTsteadyCoast.jl")
    # Self-comparison
    ϕHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyCoast", "phiHinge.txt"))[1]
    @test ϕHinge ≈ ϕHinge_ atol=SELFatol
end

@testset "Steady analysis of the Pazy wing with a FFWT at a fixed fold angle" begin
    include("examples/PazyFFWTsteadyFixedFold.jl")
    # Self-comparison
    ϕHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyFixedFold", "phiHinge.txt"))[1]
    hingeBalanceM_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyFixedFold", "hingeBalanceM.txt"))[1]
    @test ϕHinge ≈ ϕHinge_ atol=SELFatol
    @test hingeBalanceM ≈ hingeBalanceM_ atol=SELFatol
end

@testset "Steady analysis of the Pazy wing with a coasting FFWT, at varying airspeed and root pitch angle" begin
    include("examples/PazyFFWTsteadyURangeAoARangeCoast.jl")
    # Self-comparison
    ϕHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyURangeAoARangeCoast", "phiHinge.txt"))
    @test isapprox(ϕHinge, ϕHinge_, atol=SELFatol, nans=true)
end

@testset "Steady analysis of the Pazy wing with a coasting FFWT, at varying airspeed" begin
    include("examples/PazyFFWTsteadyURangeCoast.jl")
    # Self-comparison
    ϕHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyURangeCoast", "phiHinge.txt"))
    @test ϕHinge ≈ ϕHinge_ atol=SELFatol
end

@testset "Steady analysis of the Pazy wing with a FFWT at a fixed fold angle, at varying airspeed" begin
    include("examples/PazyFFWTsteadyURangeFixedFold.jl")
    # Self-comparison
    ϕHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyURangeFixedFold", "phiHinge.txt"))
    hingeMoment_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTsteadyURangeFixedFold", "hingeMoment.txt"))
    @test ϕHinge ≈ ϕHinge_ atol=SELFatol
    @test hingeMoment ≈ hingeMoment_ atol=SELFatol
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
@testset "Steady aeroelastic analysis of the clamped conventional HALE" begin
    include_example("examples/steadyAeroelastic/cHALEclampedSteady.jl")
    @test matches_test_data(x1_def, "cHALEclampedSteady", "x1_def.txt")
    @test matches_test_data(x3_def, "cHALEclampedSteady", "x3_def.txt")
end

@testset "Steady analysis of the baseline Healy free FFWT wing with varying root pitch and airspeed" begin
    include_example("examples/steadyAeroelastic/HealyBaselineFFWTsteadyAoARangeURangeCoast.jl")
    @test matches_test_data(ϕHinge, "HealyBaselineFFWTsteadyAoARangeURangeCoast", "phiHinge.txt")
    @test matches_test_data(u3Hinge, "HealyBaselineFFWTsteadyAoARangeURangeCoast", "u3Hinge.txt")
    @test matches_test_data(M2root, "HealyBaselineFFWTsteadyAoARangeURangeCoast", "M2root.txt")
end

!is_ci && @testset "Steady analysis of the Healy sideslip FFWT wing with varying flare angle and root pitch angle" begin
    include_example("examples/steadyAeroelastic/HealySideslipFFWTsteadyFlareRangeAoARangeCoast.jl")
    @test matches_test_data(ϕHinge, "HealySideslipFFWTsteadyFlareRangeAoARangeCoast", "phiHinge.txt")
end

!is_ci && @testset "Steady analysis of the Healy sideslip FFWT wing with varying flare angle, root pitch angle and airspeed" begin
    include_example("examples/steadyAeroelastic/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast.jl")
    @test matches_test_data(ϕHinge, "HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast", "phiHinge.txt")
end

@testset "Steady analysis of the Healy sideslip FFWT wing with varying wingtip twist, root pitch angle and sideslip angle" begin
    include_example("examples/steadyAeroelastic/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast.jl")
    @test matches_test_data(ϕHinge, "HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast", "phiHinge.txt")
end

!is_ci && @testset "Steady analysis of the Pazy wing with a coasting FFWT" begin
    include_example("examples/steadyAeroelastic/PazyFFWTsteadyCoast.jl")
    @test matches_test_data(ϕHinge, "PazyFFWTsteadyCoast", "phiHinge.txt")
end

!is_ci && @testset "Steady analysis of the Pazy wing with a FFWT at a fixed fold angle" begin
    include_example("examples/steadyAeroelastic/PazyFFWTsteadyFixedFold.jl")
    @test matches_test_data(ϕHinge, "PazyFFWTsteadyFixedFold", "phiHinge.txt")
    @test matches_test_data(hingeMoment, "PazyFFWTsteadyFixedFold", "hingeMoment.txt")
end

!is_ci && @testset "Steady analysis of the Pazy wing with a coasting FFWT, at varying airspeed and root pitch angle" begin
    include_example("examples/steadyAeroelastic/PazyFFWTsteadyURangeAoARangeCoast.jl")
    @test matches_test_data(ϕHinge, "PazyFFWTsteadyURangeAoARangeCoast", "phiHinge.txt")
end

!is_ci && @testset "Steady analysis of the Pazy wing with a coasting FFWT, at varying airspeed" begin
    include_example("examples/steadyAeroelastic/PazyFFWTsteadyURangeCoast.jl")
    @test matches_test_data(ϕHinge, "PazyFFWTsteadyURangeCoast", "phiHinge.txt")
end

!is_ci && @testset "Steady analysis of the Pazy wing with a FFWT at a fixed fold angle, at varying airspeed" begin
    include_example("examples/steadyAeroelastic/PazyFFWTsteadyURangeFixedFold.jl")
    @test matches_test_data(ϕHinge, "PazyFFWTsteadyURangeFixedFold", "phiHinge.txt")
    @test matches_test_data(hingeMoment, "PazyFFWTsteadyURangeFixedFold", "hingeMoment.txt")
end

@testset "Steady analysis of the Pazy wing with varying root pitch angle" begin
    include_example("examples/steadyAeroelastic/PazyWingPitchRange.jl")
    @test matches_test_data(tip_AoA, "PazyWingPitchRange", "tip_AoA.txt")
    @test matches_test_data(tip_OOP, "PazyWingPitchRange", "tip_OOP.txt")
    @test matches_test_data(tip_IP, "PazyWingPitchRange", "tip_IP.txt")
    @test matches_test_data(tip_twist, "PazyWingPitchRange", "tip_twist.txt")
end

@testset "Steady aeroelastic analysis of the sixteen-meter-wing" begin
    include_example("examples/steadyAeroelastic/SMWSteady.jl")
    @test matches_test_data(tip_u3, "SMWSteady", "tip_u3.txt")
    @test matches_test_data(tip_twist, "SMWSteady", "tip_twist.txt")
end

!is_ci && @testset "Steady analysis of the Healy LCO FFWT wing with varying airspeed and root pitch angle" begin
    include_example("examples/steadyAeroelastic/HealyLCOFFWTsteadyURangeAoARangeCoast.jl")
    @test matches_test_data(ϕHinge, "HealyLCOFFWTsteadyURangeAoARangeCoast", "phiHinge.txt")
end

@testset "Creating the conventional HALE model" begin
    include_example("examples/steadyAeroelastic/cHALEmodel.jl")
    @test typeof(cHALE).name.name == :Model
end

!is_ci && @testset "Steady analysis of the swept Pazy wing with varying root pitch angle" begin
    include_example("examples/steadyAeroelastic/sweptPazySteadyPitchRange.jl")
    @test matches_test_data(tipOOP, "sweptPazySteadyPitchRange", "tipOOP.txt")
    @test matches_test_data(tipTwist, "sweptPazySteadyPitchRange", "tipTwist.txt")
    @test matches_test_data(tipAoA, "sweptPazySteadyPitchRange", "tipAoA.txt")
end

!is_ci && @testset "Steady analysis of the swept Pazy wing with varying sweep and root pitch angles" begin
    include_example("examples/steadyAeroelastic/sweptPazySteadySweepRangePitchRange.jl")
    @test matches_test_data(tipOOP, "sweptPazySteadySweepRangePitchRange", "tipOOP.txt")
    @test matches_test_data(tipTwist, "sweptPazySteadySweepRangePitchRange", "tipTwist.txt")
end

@testset "Steady lift distribution of swept wings" begin
    include_example("examples/steadyAeroelastic/sweptWingsLiftDistribution.jl")
    @test matches_test_data(cn, "sweptWingsLiftDistribution", "cn.txt")
end
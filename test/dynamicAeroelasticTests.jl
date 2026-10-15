@testset "Dynamic analysis of the Blended-Wing-Body vehicle undergoing a checked pitch maneuver" begin
    include_example("examples/dynamicAeroelastic/BWBcheckedPitchManeuver.jl")
    @test matches_test_data(rootAoA, "BWBcheckedPitchManeuver", "rootAoA.txt")
    @test matches_test_data(Δu3, "BWBcheckedPitchManeuver", "Deltau3.txt")
end

!is_ci && @testset "Dynamic analysis of the conventional HALE aircraft undergoing a checked pitch maneuver" begin
    include_example("examples/dynamicAeroelastic/cHALECheckedPitchManeuver.jl")
    @test matches_test_data(wingAoA, "cHALECheckedPitchManeuver", "wingAoA.txt")
    @test matches_test_data(Δu3, "cHALECheckedPitchManeuver", "Deltau3.txt")
end

@testset "Dynamic analysis of the conventional HALE aircraft undergoing a coordinated turn maneuver" begin
    include_example("examples/dynamicAeroelastic/cHALECheckedRollManeuver.jl")
    @test matches_test_data(wingAoA, "cHALECheckedRollManeuver", "wingAoA.txt")
    @test matches_test_data(Δu3, "cHALECheckedRollManeuver", "Deltau3.txt")
end

!is_ci && @testset "Dynamic analysis of the baseline Healy FFWT wing under a series of one-minus-cosine gusts" begin
    include_example("examples/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating.jl")
    @test matches_test_data(M2root, "HealyBaselineFFWTOMCGustFloating", "M2root.txt")
    @test matches_test_data(ϕ, "HealyBaselineFFWTOMCGustFloating", "fold.txt")
end

@testset "Dynamic analysis of the baseline Healy FFWT wing for a specific one-minus-cosine gust" begin
    include_example("examples/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating2.jl")
    @test matches_test_data(M2root, "HealyBaselineFFWTOMCGustFloating2", "M2root.txt")
    @test matches_test_data(ϕ, "HealyBaselineFFWTOMCGustFloating2", "fold.txt")
end

@testset "Dynamic analysis of the Helios flying-wing undergoing a checked pitch maneuver" begin
    include_example("examples/dynamicAeroelastic/heliosCheckedPitchManeuver.jl")
    @test matches_test_data(rootAoA, "heliosCheckedPitchManeuver", "rootAoA.txt")
    @test matches_test_data(Δu3, "heliosCheckedPitchManeuver", "Deltau3.txt")
end

!is_ci && @testset "Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over time" begin
    include_example("examples/dynamicAeroelastic/PazyWingContinuous1DGust.jl")
    @test matches_test_data(tipAoA, "PazyWingContinuous1DGust", "tipAoA.txt")
    @test matches_test_data(tipOOP, "PazyWingContinuous1DGust", "tipOOP.txt")
    @test matches_test_data(tqSpan_cn, "PazyWingContinuous1DGust", "tqSpan_cn.txt")
    @test matches_test_data(tqSpan_cm, "PazyWingContinuous1DGust", "tqSpan_cm.txt")
    @test matches_test_data(tqSpan_ct, "PazyWingContinuous1DGust", "tqSpan_ct.txt")
end

!is_ci && @testset "Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over space" begin
    include_example("examples/dynamicAeroelastic/PazyWingContinuous1DSpaceGust.jl")
    @test matches_test_data(tipAoA, "PazyWingContinuous1DSpaceGust", "tipAoA.txt")
    @test matches_test_data(tipOOP, "PazyWingContinuous1DSpaceGust", "tipOOP.txt")
    @test matches_test_data(tqSpan_cn, "PazyWingContinuous1DSpaceGust", "tqSpan_cn.txt")
    @test matches_test_data(tqSpan_cm, "PazyWingContinuous1DSpaceGust", "tqSpan_cm.txt")
    @test matches_test_data(tqSpan_ct, "PazyWingContinuous1DSpaceGust", "tqSpan_ct.txt")
end

!is_ci && @testset "Dynamic analysis of the Pazy wing encountering a continuous, 2-dimensional gust" begin
    include_example("examples/dynamicAeroelastic/PazyWingContinuous2DSpaceGust.jl")
    @test matches_test_data(tipAoA, "PazyWingContinuous2DSpaceGust", "tipAoA.txt")
    @test matches_test_data(tipOOP, "PazyWingContinuous2DSpaceGust", "tipOOP.txt")
    @test matches_test_data(tqSpan_cn, "PazyWingContinuous2DSpaceGust", "tqSpan_cn.txt")
    @test matches_test_data(tqSpan_cm, "PazyWingContinuous2DSpaceGust", "tqSpan_cm.txt")
    @test matches_test_data(tqSpan_ct, "PazyWingContinuous2DSpaceGust", "tqSpan_ct.txt")
end

!is_ci && @testset "Dynamic analysis of the Pazy wing encountering a DARPA gust" begin
    include_example("examples/dynamicAeroelastic/PazyWingDARPAGust.jl")
    @test matches_test_data(tipAoA, "PazyWingDARPAGust", "tipAoA.txt")
    @test matches_test_data(tipOOP, "PazyWingDARPAGust", "tipOOP.txt")
    @test matches_test_data(tqSpan_cn, "PazyWingDARPAGust", "tqSpan_cn.txt")
    @test matches_test_data(tqSpan_cm, "PazyWingDARPAGust", "tqSpan_cm.txt")
    @test matches_test_data(tqSpan_ct, "PazyWingDARPAGust", "tqSpan_ct.txt")
end

!is_ci && @testset "Dynamic analysis of the Pazy wing encountering a one-minus-cosine gust" begin
    include_example("examples/dynamicAeroelastic/PazyWingOMCGust.jl")
    @test matches_test_data(tipAoA, "PazyWingOMCGust", "tipAoA.txt")
    @test matches_test_data(tipOOP, "PazyWingOMCGust", "tipOOP.txt")
    @test matches_test_data(tqSpan_cn, "PazyWingOMCGust", "tqSpan_cn.txt")
    @test matches_test_data(tqSpan_cm, "PazyWingOMCGust", "tqSpan_cm.txt")
    @test matches_test_data(tqSpan_ct, "PazyWingOMCGust", "tqSpan_ct.txt")
end

# Reduce CI time
!is_ci && @testset "Dynamic analysis of the Pazy wing with a tip impulse force" begin
    include_example("examples/dynamicAeroelastic/PazyWingTipImpulse.jl")
    @test matches_test_data(tipAoA, "PazyWingTipImpulse", "tipAoA.txt")
    @test matches_test_data(tipOOP, "PazyWingTipImpulse", "tipOOP.txt")
    @test matches_test_data(tqSpan_cn, "PazyWingTipImpulse", "tqSpan_cn.txt")
    @test matches_test_data(tqSpan_cm, "PazyWingTipImpulse", "tqSpan_cm.txt")
    @test matches_test_data(tqSpan_ct, "PazyWingTipImpulse", "tqSpan_ct.txt")
end

!is_ci && @testset "Dynamic analysis of the baseline Healy locked FFWT wing under a series of one-minus-cosine gusts" begin
    include_example("examples/dynamicAeroelastic/HealyBaselineFFWTlockedOMCGustFloating.jl")
    @test matches_test_data(M2root, "HealyBaselineFFWTlockedOMCGustFloating", "M2root.txt")
end

!is_ci && @testset "Dynamic analysis of the baseline Healy free and locked FFWT wings under one-minus-cosine gusts" begin
    include_example("examples/dynamicAeroelastic/HealyBaselineFFWTlockedVfreeOMCGustFloating.jl")
    @test matches_test_data(ΔM2peaks, "HealyBaselineFFWTlockedVfreeOMCGustFloating", "dM2peaks.txt")
end

!is_ci && @testset "Dynamic analysis of the Healy sideslip FFWT wing under a one-minus-cosine gust" begin
    include_example("examples/dynamicAeroelastic/HealySideslipFFWTOMCGustFloating.jl")
    @test matches_test_data(tipOOP, "HealySideslipFFWTOMCGustFloating", "tipOOP.txt")
    @test matches_test_data(M2_root, "HealySideslipFFWTOMCGustFloating", "M2_root.txt")
    @test matches_test_data(root_αₑ, "HealySideslipFFWTOMCGustFloating", "root_alphae.txt")
    @test matches_test_data(root_cn, "HealySideslipFFWTOMCGustFloating", "root_cn.txt")
end

!is_ci && @testset "Dynamic analysis of the Pazy FFWT wing under a one-minus-cosine gust" begin
    include_example("examples/dynamicAeroelastic/PazyFFWTOMCGustFloating.jl")
    @test matches_test_data(ϕ, "PazyFFWTOMCGustFloating", "fold.txt")
    @test matches_test_data(tipOOP, "PazyFFWTOMCGustFloating", "tipOOP.txt")
    @test matches_test_data(M2_root, "PazyFFWTOMCGustFloating", "M2_root.txt")
    @test matches_test_data(tqSpan_cn, "PazyFFWTOMCGustFloating", "tqSpan_cn.txt")
end

!is_ci && @testset "Dynamic analysis of a rotating wing using the dynamic stall model" begin
    include_example("examples/dynamicAeroelastic/rotorDStest.jl")
    @test matches_test_data(cl, "rotorDStest", "cl.txt")
    @test matches_test_data(cm, "rotorDStest", "cm.txt")
    @test matches_test_data(cdrag, "rotorDStest", "cdrag.txt")
end

@testset "Dynamic analysis of a pitching wing with dynamic stall models" begin
    include_example("examples/dynamicAeroelastic/wingDStest.jl")
    @test matches_test_data(cl, "wingDStest", "cl.txt")
    @test matches_test_data(cm, "wingDStest", "cm.txt")
    @test matches_test_data(cdrag, "wingDStest", "cdrag.txt")
end
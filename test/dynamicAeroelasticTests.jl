# Dynamic aeroelastic problems

@testset "Dynamic analysis of the Blended-Wing-Body vehicle undergoing a checked pitch maneuver" begin
    include("examples/BWBcheckedPitchManeuver.jl")
    # Self-comparison
    rootAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "BWBcheckedPitchManeuver", "rootAoA.txt"))
    Δu3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "BWBcheckedPitchManeuver", "Deltau3.txt"))
    @test rootAoA ≈ rootAoA_ atol=SELFatol
    @test Δu3 ≈ Δu3_ atol=SELFatol
end

## Reduce CI time
# @testset "Dynamic analysis of the conventional HALE aircraft undergoing a checked pitch maneuver" begin
#     include("examples/conventionalHALECheckedPitchManeuver.jl")
#     # Self-comparison
#     wingAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALECheckedPitchManeuver", "wingAoA.txt"))
#     Δu3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALECheckedPitchManeuver", "Deltau3.txt"))
#     @test wingAoA ≈ wingAoA_ atol=SELFatol
#     @test Δu3 ≈ Δu3_ atol=SELFatol
# end

@testset "Dynamic analysis of the conventional HALE aircraft undergoing a coordinated turn maneuver" begin
    include("examples/conventionalHALECheckedRollManeuver.jl")
    # Self-comparison
    wingAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALECheckedRollManeuver", "wingAoA.txt"))
    Δu3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALECheckedRollManeuver", "Deltau3.txt"))
    @test wingAoA ≈ wingAoA_ atol=SELFatol
    @test Δu3 ≈ Δu3_ atol=SELFatol
end

@testset "Dynamic analysis of the baseline Healy FFWT wing under a series of one-minus-cosine gusts" begin
    include("examples/HealyBaselineFFWTOMCGustFloating.jl")
    # Self-comparison
    for (i,ω) in enumerate(ωRange)
        M2root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyBaselineFFWTOMCGustFloating", string("M2root_omega",ω,".txt")))
        ϕ_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyBaselineFFWTOMCGustFloating", string("fold_omega",ω,".txt")))
        @test M2root[i] ≈ M2root_ atol=SELFatol
        @test ϕ[i] ≈ ϕ_ atol=SELFatol
    end
end

@testset "Dynamic analysis of the Helios flying-wing undergoing a checked pitch maneuver" begin
    include("examples/heliosCheckedPitchManeuver.jl")
    # Self-comparison
    rootAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosCheckedPitchManeuver", "rootAoA.txt"))
    Δu3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosCheckedPitchManeuver", "Deltau3.txt"))
    @test rootAoA ≈ rootAoA_ atol=SELFatol
    @test Δu3 ≈ Δu3_ atol=SELFatol
end

## Reduce CI time
# @testset "Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over time" begin
#     include("examples/PazyWingContinuous1DGust.jl")
#     # Self-comparison
#     tipAoA_ = vec(readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DGust", "tipAoA.txt")))
#     tipOOP_ = vec(readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DGust", "tipOOP.txt")))
#     tqSpan_cn_ = vec(readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DGust", "tqSpan_cn.txt")))
#     tqSpan_cm_ = vec(readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DGust", "tqSpan_cm.txt")))
#     tqSpan_ct_ = vec(readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DGust", "tqSpan_ct.txt")))
#     @test tipAoA[end] ≈ tipAoA_[end] atol=5e-3
#     @test tipOOP[end] ≈ tipOOP_[end] atol=5e-3
#     @test tqSpan_cn[end] ≈ tqSpan_cn_[end] atol=5e-3
#     @test tqSpan_cm[end] ≈ tqSpan_cm_[end] atol=5e-3
#     @test tqSpan_ct[end] ≈ tqSpan_ct_[end] atol=5e-3
# end

## Reduce CI time
# @testset "Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over space" begin
#     include("examples/PazyWingContinuous1DSpaceGust.jl")
#     # Self-comparison
#     tipAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DSpaceGust", "tipAoA.txt"))
#     tipOOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DSpaceGust", "tipOOP.txt"))
#     tqSpan_cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DSpaceGust", "tqSpan_cn.txt"))
#     tqSpan_cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DSpaceGust", "tqSpan_cm.txt"))
#     tqSpan_ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DSpaceGust", "tqSpan_ct.txt"))
#     @test tipAoA[end] ≈ tipAoA_[end] atol=5e-3
#     @test tipOOP[end] ≈ tipOOP_[end] atol=5e-3
#     @test tqSpan_cn[end] ≈ tqSpan_cn_[end] atol=5e-3
#     @test tqSpan_cm[end] ≈ tqSpan_cm_[end] atol=5e-3
#     @test tqSpan_ct[end] ≈ tqSpan_ct_[end] atol=5e-3
# end

## Reduce CI time
# @testset "Dynamic analysis of the Pazy wing encountering a continuous, 2-dimensional gust" begin
#     include("examples/PazyWingContinuous2DSpaceGust.jl")
#     # Self-comparison
#     tipAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous2DSpaceGust", "tipAoA.txt"))
#     tipOOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous2DSpaceGust", "tipOOP.txt"))
#     tqSpan_cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous2DSpaceGust", "tqSpan_cn.txt"))
#     tqSpan_cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous2DSpaceGust", "tqSpan_cm.txt"))
#     tqSpan_ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous2DSpaceGust", "tqSpan_ct.txt"))
#     @test tipAoA[end] ≈ tipAoA_[end] atol=5e-3
#     @test tipOOP[end] ≈ tipOOP_[end] atol=5e-3
#     @test tqSpan_cn[end] ≈ tqSpan_cn_[end] atol=5e-3
#     @test tqSpan_cm[end] ≈ tqSpan_cm_[end] atol=5e-3
#     @test tqSpan_ct[end] ≈ tqSpan_ct_[end] atol=5e-3
# end

## Reduce CI time
# @testset "Dynamic analysis of the Pazy wing encountering a DARPA gust" begin
#     include("examples/PazyWingDARPAGust.jl")
#     # Self-comparison
#     tipAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingDARPAGust", "tipAoA.txt"))
#     tipOOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingDARPAGust", "tipOOP.txt"))
#     tqSpan_cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingDARPAGust", "tqSpan_cn.txt"))
#     tqSpan_cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingDARPAGust", "tqSpan_cm.txt"))
#     tqSpan_ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingDARPAGust", "tqSpan_ct.txt"))
#     @test tipAoA[end] ≈ tipAoA_[end] atol=5e-3
#     @test tipOOP[end] ≈ tipOOP_[end] atol=5e-3
#     @test tqSpan_cn[end] ≈ tqSpan_cn_[end] atol=5e-3
#     @test tqSpan_cm[end] ≈ tqSpan_cm_[end] atol=5e-3
#     @test tqSpan_ct[end] ≈ tqSpan_ct_[end] atol=5e-3
# end

## Reduce CI time
# @testset "Dynamic analysis of the Pazy wing encountering a one-minus-cosine gust" begin
#     include("examples/PazyWingOMCGust.jl")
#     # Self-comparison
#     tipAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingOMCGust", "tipAoA.txt"))
#     tipOOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingOMCGust", "tipOOP.txt"))
#     tqSpan_cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingOMCGust", "tqSpan_cn.txt"))
#     tqSpan_cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingOMCGust", "tqSpan_cm.txt"))
#     tqSpan_ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingOMCGust", "tqSpan_ct.txt"))
#     @test tipAoA[end] ≈ tipAoA_[end] atol=5e-3
#     @test tipOOP[end] ≈ tipOOP_[end] atol=5e-3
#     @test tqSpan_cn[end] ≈ tqSpan_cn_[end] atol=5e-3
#     @test tqSpan_cm[end] ≈ tqSpan_cm_[end] atol=5e-3
#     @test tqSpan_ct[end] ≈ tqSpan_ct_[end] atol=5e-3
# end

## Reduce CI time
# @testset "Dynamic analysis of the Pazy wing with a tip impulse force" begin
#     include("examples/PazyWingTipImpulse.jl")
#     # Self-comparison
#     tipAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTipImpulse", "tipAoA.txt"))
#     tipOOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTipImpulse", "tipOOP.txt"))
#     tqSpan_cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTipImpulse", "tqSpan_cn.txt"))
#     tqSpan_cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTipImpulse", "tqSpan_cm.txt"))
#     tqSpan_ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTipImpulse", "tqSpan_ct.txt"))
#     @test tipAoA[end] ≈ tipAoA_[end] atol=5e-3
#     @test tipOOP[end] ≈ tipOOP_[end] atol=5e-3
#     @test tqSpan_cn[end] ≈ tqSpan_cn_[end] atol=5e-3
#     @test tqSpan_cm[end] ≈ tqSpan_cm_[end] atol=5e-3
#     @test tqSpan_ct[end] ≈ tqSpan_ct_[end] atol=5e-3
# end
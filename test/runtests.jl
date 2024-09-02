using Test, BenchmarkTools, DelimitedFiles

# Default absolute tolerance for self-comparison
SELFatol = 1e-4

# @testset "Static analysis of an arch under a dead pressure load" begin
#     include("examples/archUnderDeadPressure.jl")
#     # Self-comparison
#     mid_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "archUnderDeadPressure", "mid_u3.txt"))
#     @test mid_u3 ≈ mid_u3_ atol=SELFatol
# end

# @testset "Static analysis of an arch under a follower pressure load" begin
#     include("examples/archUnderFollowerPressure.jl")
#     # Self-comparison
#     mid_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "archUnderFollowerPressure", "mid_u3.txt"))
#     @test mid_u3 ≈ mid_u3_ atol=SELFatol
# end

# @testset "Dynamic analysis of the axial vibration of a beam under a traction force applied suddenly" begin
#     include("examples/axialTractionCantilever.jl")
#     # Reference comparison
#     @test u1_08[end]*1e3 ≈ u1_08_ref[end] rtol=1e-2
#     @test u1_10[end]*1e3 ≈ u1_10_ref[end] rtol=2e-2
#     # Self-comparison
#     u1_08_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "axialTractionCantilever", "u1_08.txt"))
#     u1_10_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "axialTractionCantilever", "u1_10.txt"))
#     @test u1_08 ≈ u1_08_ atol=SELFatol
#     @test u1_10 ≈ u1_10_ atol=SELFatol
# end

# @testset "Modal analysis of the axial vibration of a beam under clamped-clamped boundary conditions" begin
#     include("examples/beamAxialVibrationCC.jl")
#     # Analytical comparison
#     @test freqs ≈ freqsAnalytical rtol=5e-3
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamAxialVibrationCC", "freqs.txt"))
#     u1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamAxialVibrationCC", "u1_modeShapes.txt"))
#     @test freqs ≈ freqs_ atol=SELFatol
#     # @test hcat(u1_modeShapes...)' ≈ u1_modeShapes_ atol=SELFatol
# end

# @testset "Modal analysis of the axial vibration of a beam under clamped-free boundary conditions" begin
#     include("examples/beamAxialVibrationCF.jl")
#     # Analytical comparison
#     @test freqs ≈ freqsAnalytical rtol=5e-3
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamAxialVibrationCF", "freqs.txt"))
#     u1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamAxialVibrationCF", "u1_modeShapes.txt"))
#     @test freqs ≈ freqs_ atol=SELFatol
#     # @test hcat(u1_modeShapes...)' ≈ u1_modeShapes_ atol=SELFatol
# end

# @testset "Modal analysis of the axial vibration of a beam under free-free boundary conditions" begin
#     include("examples/beamAxialVibrationFF.jl")
#     # Analytical comparison
#     @test freqs ≈ freqsAnalytical rtol=5e-3
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamAxialVibrationFF", "freqs.txt"))
#     u1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamAxialVibrationFF", "u1_modeShapes.txt"))
#     @test freqs ≈ freqs_ atol=SELFatol
#     # @test hcat(u1_modeShapes...)' ≈ u1_modeShapes_ atol=SELFatol
# end

# @testset "Modal analysis of the bending vibration of a beam under clamped-clamped boundary conditions" begin
#     include("examples/beamBendingVibrationCC.jl")
#     # Analytical comparison
#     @test freqs ≈ freqsAnalytical rtol=1e-2
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCC", "freqs.txt"))
#     u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCC", "u3_modeShapes.txt"))
#     @test freqs ≈ freqs_ atol=SELFatol
#     # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
# end

# @testset "Modal analysis of the bending vibration of a beam under clamped-free boundary conditions" begin
#     include("examples/beamBendingVibrationCF.jl")
#     # Analytical comparison
#     @test freqs ≈ freqsAnalytical rtol=1e-2
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCF", "freqs.txt"))
#     u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCF", "u3_modeShapes.txt"))
#     @test freqs ≈ freqs_ atol=SELFatol
#     # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
# end

# @testset "Modal analysis of the bending vibration of a beam under clamped-pinned boundary conditions" begin
#     include("examples/beamBendingVibrationCP.jl")
#     # Analytical comparison
#     @test freqs ≈ freqsAnalytical rtol=1e-2
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCP", "freqs.txt"))
#     u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCP", "u3_modeShapes.txt"))
#     @test freqs ≈ freqs_ atol=SELFatol
#     # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
# end

# @testset "Modal analysis of the bending vibration of a beam under clamped-sliding boundary conditions" begin
#     include("examples/beamBendingVibrationCS.jl")
#     # Analytical comparison
#     @test freqs ≈ freqsAnalytical rtol=1e-2
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCS", "freqs.txt"))
#     u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCS", "u3_modeShapes.txt"))
#     @test freqs ≈ freqs_ atol=SELFatol
#     # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
# end

# @testset "Modal analysis of the bending vibration of a beam under free-free boundary conditions" begin
#     include("examples/beamBendingVibrationFF.jl")
#     # Analytical comparison
#     @test freqs ≈ freqsAnalytical rtol=1e-2
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationFF", "freqs.txt"))
#     u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationFF", "u3_modeShapes.txt"))
#     @test freqs ≈ freqs_ atol=SELFatol
#     # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
# end

# @testset "Modal analysis of the bending vibration of a beam under pinned-pinned boundary conditions" begin
#     include("examples/beamBendingVibrationPP.jl")
#     # Analytical comparison
#     @test freqs ≈ freqsAnalytical rtol=1e-2
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationPP", "freqs.txt"))
#     u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationPP", "u3_modeShapes.txt"))
#     @test freqs ≈ freqs_ atol=SELFatol
#     # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
# end

# @testset "Modal analysis of the torsional vibration of a beam under clamped-clamped boundary conditions" begin
#     include("examples/beamTorsionalVibrationCC.jl")
#     # Analytical comparison
#     @test freqs ≈ freqsAnalytical rtol=5e-3
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamTorsionalVibrationCC", "freqs.txt"))
#     p1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamTorsionalVibrationCC", "p1_modeShapes.txt"))
#     @test freqs ≈ freqs_ atol=SELFatol
#     # @test hcat(p1_modeShapes...)' ≈ p1_modeShapes_ atol=SELFatol
# end

# @testset "Modal analysis of the torsional vibration of a beam under clamped-free boundary conditions" begin
#     include("examples/beamTorsionalVibrationCF.jl")
#     # Analytical comparison
#     @test freqs ≈ freqsAnalytical rtol=5e-3
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamTorsionalVibrationCF", "freqs.txt"))
#     p1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamTorsionalVibrationCF", "p1_modeShapes.txt"))
#     @test freqs ≈ freqs_ atol=SELFatol
#     # @test hcat(p1_modeShapes...)' ≈ p1_modeShapes_ atol=SELFatol
# end

# @testset "Modal analysis of the torsional vibration of a beam under free-free boundary conditions" begin
#     include("examples/beamTorsionalVibrationFF.jl")
#     # Analytical comparison
#     @test freqs ≈ freqsAnalytical rtol=5e-3
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamTorsionalVibrationFF", "freqs.txt"))
#     p1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamTorsionalVibrationFF", "p1_modeShapes.txt"))
#     @test freqs ≈ freqs_ atol=SELFatol
#     # @test hcat(p1_modeShapes...)' ≈ p1_modeShapes_ atol=SELFatol
# end

# @testset "Dynamic analysis of the free response of a beam clamped at both ends and subjected to an initial displacement profile" begin
#     include("examples/biclampedBeam.jl")
#     # Analytical comparison
#     @test u3_mid ≈ u3_mid_analytic rtol=1e-2
#     @test V3_mid ≈ V3_mid_analytic rtol=5e-2
#     @test θ2_quarter ≈ θ2_quarter_analytic rtol=2e-2
#     @test Ω2_quarter ≈ Ω2_quarter_analytic rtol=10e-2
#     # Self-comparison
#     u3_mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedBeam", "u3_mid.txt"))
#     V3_mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedBeam", "V3_mid.txt"))
#     Vdot3_mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedBeam", "Vdot3_mid.txt"))
#     θ2_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedBeam", "theta2_quarter.txt"))
#     Ω2_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedBeam", "Omega2_quarter.txt"))
#     Ωdot2_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedBeam", "Omegadot2_quarter.txt"))
#     @test u3_mid ≈ u3_mid_ atol=SELFatol
#     @test V3_mid ≈ V3_mid_ atol=SELFatol
#     @test Vdot3_mid ≈ Vdot3_mid_ atol=SELFatol
#     @test θ2_quarter ≈ θ2_quarter_ atol=SELFatol
#     @test Ω2_quarter ≈ Ω2_quarter_ atol=SELFatol
#     @test Ωdot2_quarter ≈ Ωdot2_quarter_ atol=SELFatol
# end

# @testset "Flutter analysis of the Blended-Wing-Body flying wing" begin
#     include("examples/BWBflutter.jl")
#     # Self-comparison
#     freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "BWBflutter", "freqs.txt"))
#     damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "BWBflutter", "damps.txt"))
#     @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
#     @test hcat(damps...)' ≈ damps_ atol=SELFatol
# end

@testset "Trim analysis of the Blended-Wing-Body flying wing in free flight" begin
    include("examples/BWBtrim.jl")
    # Self-comparison
    trimAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "BWBtrim", "trimAoA.txt"))
    trimThrust_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "BWBtrim", "trimThrust.txt"))
    trimδ_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "BWBtrim", "trimDelta.txt"))
    @test trimAoA ≈ trimAoA_ atol=SELFatol
    @test trimThrust ≈ trimThrust_ atol=SELFatol
    @test trimδ ≈ trimδ_ atol=SELFatol
end

@testset "Static analysis of a cantilever beam bending under self weight" begin
    include("examples/cantileverUnderSelfWeight.jl")
    # Self-comparison
    u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverUnderSelfWeight", "u1.txt"))
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverUnderSelfWeight", "u3.txt"))
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverUnderSelfWeight", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverUnderSelfWeight", "M2.txt"))
    @test u1 ≈ u1_ atol=SELFatol
    @test u3 ≈ u3_ atol=SELFatol
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end

@testset "Modal analysis of a cantilever beam with a tip axial inertia" begin
    include("examples/cantileverWithTipAxialMassEigen.jl")
    # Analytical comparison
    @test freqsNorm ≈ freqsNormAnalytical rtol=2e-2
    # Self-comparison
    freqsNorm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipAxialMassEigen", "freqsNorm.txt"))
    u1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipAxialMassEigen", "u1_modeShapes.txt"))
    @test freqsNorm ≈ freqsNorm_ atol=SELFatol
    @test hcat(u1_modeShapes...)' ≈ u1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of a cantilever beam with a tip axial spring" begin
    include("examples/cantileverWithTipAxialSpringEigen.jl")
    # Analytical comparison
    @test freqsNorm ≈ freqsNormAnalytical rtol=5e-3
    # Self-comparison
    freqsNorm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipAxialSpringEigen", "freqsNorm.txt"))
    u1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipAxialSpringEigen", "u1_modeShapes.txt"))
    @test freqsNorm ≈ freqsNorm_ atol=SELFatol
    @test hcat(u1_modeShapes...)' ≈ u1_modeShapes_ atol=SELFatol
end

@testset "Static analysis of a cantilever beam with a tip spring in bending" begin
    include("examples/cantileverWithTipSpring.jl")
    # Analytical comparison
    @test u3[end] ≈ 2.5e-3 rtol=1e-3
    # Self-comparison
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipSpring", "u3.txt"))
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipSpring", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipSpring", "M2.txt"))
    @test u3 ≈ u3_ atol=SELFatol
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end

@testset "Modal analysis of a cantilever beam with a tip spring in bending" begin
    include("examples/cantileverWithTipSpringEigen.jl")
    # Analytical comparison
    @test freqsNorm ≈ freqsNormRef rtol=1e-2
    # Self-comparison
    freqsNorm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipSpringEigen", "freqsNorm.txt"))
    u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipSpringEigen", "u3_modeShapes.txt"))
    @test freqsNorm ≈ freqsNorm_ atol=SELFatol
    @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of a cantilever beam with a tip torsional inertia" begin
    include("examples/cantileverWithTipTorsionalInertiaEigen.jl")
    # Analytical comparison
    @test freqsNorm ≈ freqsNormAnalytical rtol=2e-2
    # Self-comparison
    freqsNorm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipTorsionalInertiaEigen", "freqsNorm.txt"))
    p1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipTorsionalInertiaEigen", "p1_modeShapes.txt"))
    @test freqsNorm ≈ freqsNorm_ atol=SELFatol
    @test hcat(p1_modeShapes...)' ≈ p1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of a cantilever beam with a tip torsional spring" begin
    include("examples/cantileverWithTipTorsionalSpringEigen.jl")
    # Analytical comparison
    @test freqsNorm ≈ freqsNormAnalytical rtol=1e-2
    # Self-comparison
    freqsNorm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipTorsionalSpringEigen", "freqsNorm.txt"))
    p1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipTorsionalSpringEigen", "p1_modeShapes.txt"))
    @test freqsNorm ≈ freqsNorm_ atol=SELFatol
    @test hcat(p1_modeShapes...)' ≈ p1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis a beam clamped at one end, simply-supported at the other and with a tip inertia" begin
    include("examples/clampedSSBeamWIthTipInertia.jl")
    # Reference (analytical) comparison
    @test freqs[1] ≈ freqRef rtol=2e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedSSBeamWIthTipInertia", "freqs.txt"))
    u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedSSBeamWIthTipInertia", "u3_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
end

@testset "Dynamic analysis of a composite cantilever beam under a tip sinusoidal load" begin
    include("examples/compositeCantilever.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantilever", "u1_tip.txt"))
    u2_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantilever", "u2_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantilever", "u3_tip.txt"))
    p1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantilever", "p1_tip.txt"))
    p2_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantilever", "p2_tip.txt"))
    p3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantilever", "p3_tip.txt"))
    F1_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantilever", "F1_root.txt"))
    F2_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantilever", "F2_root.txt"))
    F3_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantilever", "F3_root.txt"))
    M1_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantilever", "M1_root.txt"))
    M2_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantilever", "M2_root.txt"))
    M3_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantilever", "M3_root.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol
    @test u2_tip ≈ u2_tip_ atol=SELFatol
    @test u3_tip ≈ u3_tip_ atol=SELFatol
    @test p1_tip ≈ p1_tip_ atol=SELFatol
    @test p2_tip ≈ p2_tip_ atol=SELFatol
    @test p3_tip ≈ p3_tip_ atol=SELFatol
    @test F1_root ≈ F1_root_ atol=SELFatol
    @test F2_root ≈ F2_root_ atol=SELFatol
    @test F3_root ≈ F3_root_ atol=SELFatol
    @test M1_root ≈ M1_root_ atol=SELFatol
    @test M2_root ≈ M2_root_ atol=SELFatol
    @test M3_root ≈ M3_root_ atol=SELFatol
end

@testset "Static analysis of composite laminates subjected to tip loads" begin
    include("examples/compositeCantileverMD.jl")
    # Self-comparison
    for i in 1:3
        u1_500mm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantileverMD", string("u1_500mm_b",i,".txt")))
        u2_500mm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantileverMD", string("u2_500mm_b",i,".txt")))
        u3_500mm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "compositeCantileverMD", string("u3_500mm_b",i,".txt")))
        @test hcat(u1_500mm[i,:]...)' ≈ u1_500mm_ atol=SELFatol
        @test hcat(u2_500mm[i,:]...)' ≈ u2_500mm_ atol=SELFatol
        @test hcat(u3_500mm[i,:]...)' ≈ u3_500mm_ atol=SELFatol
    end
end

@testset "Dynamic analysis the conventional HALE aircraft undergoing a checked pitch maneuver" begin
    include("examples/conventionalHALECheckedPitchManeuver.jl")
    # Self-comparison
    rootAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALECheckedPitchManeuver", "rootAoA.txt"))
    Δu3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALECheckedPitchManeuver", "Deltau3.txt"))
    @test rootAoA ≈ rootAoA_ atol=SELFatol
    @test Δu3 ≈ Δu3_ atol=SELFatol
end

@testset "Trim analysis the conventional HALE aircraft in free flight (considering aerodynamics from stabilizers and thrust)" begin
    include("examples/conventionalHALEfullTrim.jl")
    # Self-comparison
    trimAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEfullTrim", "trimAoA.txt"))
    trimThrust_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEfullTrim", "trimThrust.txt"))
    trimδ_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEfullTrim", "trimDelta.txt"))
    @test trimAoA ≈ trimAoA_ atol=SELFatol
    @test trimThrust ≈ trimThrust_ atol=SELFatol
    @test trimδ ≈ trimδ_ atol=SELFatol
end

@testset "Flutter analysis the conventional HALE aircraft in free flight with structural stiffness as the varying parameter" begin
    include("examples/conventionalHALELambdaRange.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALELambdaRange", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALELambdaRange", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Flutter analysis the conventional HALE aircraft in free flight with airspeed and structural stiffness as the varying parameters" begin
    include("examples/conventionalHALELURange.jl")
    # Self-comparison
    for (i,λ) in enumerate(λRange)
        freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALELURange", string("freqs",i,".txt")))
        damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALELURange", string("damps",i,".txt")))
        @test hcat(freqs[i,:]...)' ≈ freqs_ atol=SELFatol
        @test hcat(damps[i,:]...)' ≈ damps_ atol=SELFatol
    end
end

@testset "Trim analysis the conventional HALE aircraft in free flight at rigid and flexible configurations (neglecting aerodynamics from stabilizers and thrust)" begin
    include("examples/conventionalHALEtrim.jl")
    # Reference comparison
    @test trimAoA[1,1] ≈ trimAoAERef[2,1] rtol=2e-2
    @test trimAoA[2,1] ≈ trimAoARRef[2,1] rtol=2e-2
    # Self-comparison
    trimAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEtrim", "trimAoA.txt"))
    @test trimAoA ≈ trimAoA_ atol=SELFatol
end

@testset "Flutter analysis the conventional HALE aircraft in free flight with airspeed as the varying parameter" begin
    include("examples/conventionalHALEURange.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEURange", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "conventionalHALEURange", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Static analysis of a curved cantilever subjected to a tip dead force" begin
    include("examples/curvedCantileverDeadLoad.jl")
    # Reference comparison
    @test tip_u1[end] ≈ u1Ref rtol=1e-2
    @test tip_u2[end] ≈ u2Ref rtol=2e-2
    @test tip_u3[end] ≈ u3Ref rtol=1e-2
    # Self-comparison
    tip_u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverDeadLoad", "tip_u1.txt"))[1]
    tip_u2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverDeadLoad", "tip_u2.txt"))[1]
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverDeadLoad", "tip_u3.txt"))[1]
    @test tip_u1 ≈ tip_u1_ atol=SELFatol
    @test tip_u2 ≈ tip_u2_ atol=SELFatol
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
end

@testset "Dynamic analysis of a curved cantilever subjected to a tip follower force" begin
    include("examples/curvedCantileverDynamicFollower.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverDynamicFollower", "u1_tip.txt"))
    u2_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverDynamicFollower", "u2_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverDynamicFollower", "u3_tip.txt"))
    F1_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverDynamicFollower", "F1_root.txt"))
    F2_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverDynamicFollower", "F2_root.txt"))
    F3_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverDynamicFollower", "F3_root.txt"))
    M1_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverDynamicFollower", "M1_root.txt"))
    M2_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverDynamicFollower", "M2_root.txt"))
    M3_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverDynamicFollower", "M3_root.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol
    @test u2_tip ≈ u2_tip_ atol=SELFatol
    @test u3_tip ≈ u3_tip_ atol=SELFatol
    @test F1_root ≈ F1_root_ atol=SELFatol
    @test F2_root ≈ F2_root_ atol=SELFatol
    @test F3_root ≈ F3_root_ atol=SELFatol
    @test M1_root ≈ M1_root_ atol=SELFatol
    @test M2_root ≈ M2_root_ atol=SELFatol
    @test M3_root ≈ M3_root_ atol=SELFatol
end

@testset "Static analysis of a curved cantilever subjected to a tip follower force" begin
    include("examples/curvedCantileverStaticFollower.jl")
    # Reference comparison
    @test -tip_u1[end] ≈ u1_ref[2,end] rtol=5e-2
    @test tip_u2[end] ≈ u2_ref[2,end] rtol=5e-2
    @test tip_u3[end] ≈ u3_ref[2,end] atol=6
    # Self-comparison
    tip_u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverStaticFollower", "tip_u1.txt"))
    tip_u2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverStaticFollower", "tip_u2.txt"))
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "curvedCantileverStaticFollower", "tip_u3.txt"))
    @test tip_u1 ≈ tip_u1_ atol=SELFatol
    @test tip_u2 ≈ tip_u2_ atol=SELFatol
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
end

@testset "Static analysis of a cantilever with distributed follower force" begin
    include("examples/distributedLoadCantilever.jl")
    # Reference comparison
    @test -tip_u1[end]/L ≈ u1Ref[1,end] rtol=0.1
    @test tip_u3[end]/L ≈ u3Ref[1,end] atol=0.1
    @test -tip_angle[end]/π ≈ θRef[1,end] rtol=0.1
    # Self-comparison
    tip_u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "distributedLoadCantilever", "tip_u1.txt"))
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "distributedLoadCantilever", "tip_u3.txt"))
    tip_angle_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "distributedLoadCantilever", "tip_angle.txt"))
    @test tip_u1 ≈ tip_u1_ atol=SELFatol
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
    @test tip_angle ≈ tip_angle_ atol=SELFatol
end

@testset "Dynamic analysis of a double pendulum released from rest" begin
    include("examples/doublePendulum.jl")
    # Self-comparison
    u1_hinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "doublePendulum", "u1_hinge.txt"))
    u3_hinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "doublePendulum", "u3_hinge.txt"))
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "doublePendulum", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "doublePendulum", "u3_tip.txt"))
    @test u1_hinge ≈ u1_hinge_ atol=SELFatol
    @test u3_hinge ≈ u3_hinge_ atol=SELFatol
    @test u1_tip ≈ u1_tip_ atol=SELFatol
    @test u3_tip ≈ u3_tip_ atol=SELFatol
end

@testset "Dynamic analysis of a harmonically pitching airfoil, using the dynamic stall model" begin
    include("examples/DSModelTest.jl")
    # Self-comparison
    cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "DSModelTest", "cn.txt"))
    cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "DSModelTest", "cm.txt"))
    ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "DSModelTest", "ct.txt"))
    @test cn ≈ cn_ atol=SELFatol
    @test cm ≈ cm_ atol=SELFatol
    @test ct ≈ ct_ atol=SELFatol
end

@testset "Dynamic analysis of a right-angled frame subjected to an out-of-plane force" begin
    include("examples/elbowFrame.jl")
    # Analytical comparison
    @test u3_elbow[end] ≈ u3ElbowRef[2,end] atol=0.5
    @test u3_tip[end] ≈ u3TipRef[2,end] atol=1
    # Self-comparison
    u3_elbow_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "elbowFrame", "u3_elbow.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "elbowFrame", "u3_tip.txt"))
    @test u3_elbow ≈ u3_elbow_ atol=SELFatol
    @test u3_tip ≈ u3_tip_ atol=SELFatol
end

@testset "Dynamic analysis of an airfoils with harmonic flap deflection profile" begin
    include("examples/flapOscillation.jl")
    # Self-comparison
    cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flapOscillation", "cn.txt"))
    cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flapOscillation", "cm.txt"))
    @test cn ≈ cn_ atol=SELFatol
    @test cm ≈ cm_ atol=SELFatol
end

@testset "Dynamic analysis of two airfoils with linked harmonic flap deflection profiles" begin
    include("examples/flapOscillationLinked.jl")
    # Self-comparison
    cnMaster_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flapOscillationLinked", "cnMaster.txt"))
    cmMaster_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flapOscillationLinked", "cmMaster.txt"))
    cnSlave_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flapOscillationLinked", "cnSlave.txt"))
    cmSlave_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flapOscillationLinked", "cmSlave.txt"))
    @test cnMaster ≈ cnMaster_ atol=SELFatol
    @test cmMaster ≈ cmMaster_ atol=SELFatol
    @test cnSlave ≈ cnSlave_ atol=SELFatol
    @test cmSlave ≈ cmSlave_ atol=SELFatol
end

@testset "Dynamic analysis of a flexible beam subjected to loads yielding two-dimensional motion" begin
    include("examples/flyingFlexibleBeam2D.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingFlexibleBeam2D", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingFlexibleBeam2D", "u3_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol
    @test u3_tip ≈ u3_tip_ atol=SELFatol
end

@testset "Dynamic analysis of a hinged beam in free flight" begin
    include("examples/flyingScissors.jl")
    # Self-comparison
    u1_tipA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingScissors", "u1_tipA.txt"))
    u3_tipA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingScissors", "u3_tipA.txt"))
    u1_tipB_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingScissors", "u1_tipB.txt"))
    u3_tipB_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingScissors", "u3_tipB.txt"))
    u1_hinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingScissors", "u1_hinge.txt"))
    u3_hinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingScissors", "u3_hinge.txt"))
    @test u1_tipA ≈ u1_tipA_ atol=SELFatol
    @test u3_tipA ≈ u3_tipA_ atol=SELFatol
    @test u1_tipB ≈ u1_tipB_ atol=SELFatol
    @test u3_tipB ≈ u3_tipB_ atol=SELFatol
    @test u1_hinge ≈ u1_hinge_ atol=SELFatol
    @test u3_hinge ≈ u3_hinge_ atol=SELFatol
end

@testset "Dynamic analysis of a very flexible beam subjected to loads yielding two-dimensional motion" begin
    include("examples/flyingSpaghetti2D.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingSpaghetti2D", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingSpaghetti2D", "u3_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol
    @test u3_tip ≈ u3_tip_ atol=SELFatol
end

@testset "Dynamic analysis of a very flexible beam subjected to loads yielding tri-dimensional motion" begin
    include("examples/flyingSpaghetti3D.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingSpaghetti3D", "u1_tip.txt"))
    u2_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingSpaghetti3D", "u2_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingSpaghetti3D", "u3_tip.txt"))
    θ_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingSpaghetti3D", "theta_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol
    @test u2_tip ≈ u2_tip_ atol=SELFatol
    @test u3_tip ≈ u3_tip_ atol=SELFatol
    @test θ_tip ≈ θ_tip_ atol=SELFatol
end

@testset "Trim analysis (reaction loads check) of a beam loaded at the middle" begin
    include("examples/freeBeamTrim.jl")
    # Analytical comparison
    @test trimForce ≈ trimForceAnalytical rtol=1e-4
    # Self-comparison
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "freeBeamTrim", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "freeBeamTrim", "M2.txt"))
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end

@testset "Dynamic analysis the Helios flying-wing undergoing a checked pitch maneuver" begin
    include("examples/heliosCheckedPitchManeuver.jl")
    # Self-comparison
    rootAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosCheckedPitchManeuver", "rootAoA.txt"))
    Δu3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosCheckedPitchManeuver", "Deltau3.txt"))
    @test rootAoA ≈ rootAoA_ atol=SELFatol
    @test Δu3 ≈ Δu3_ atol=SELFatol
end

@testset "Flutter analysis the Helios flying-wing in free flight with payload and structural stiffness as the varying parameters" begin
    include("examples/heliosFlutterPLambdaRange.jl")
    # Self-comparison
    for (i,λ) in enumerate(λRange)
        freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosFlutterPLambdaRange", string("freqs",i,".txt")))
        damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosFlutterPLambdaRange", string("damps",i,".txt")))
        @test hcat(freqs[i,:]...)' ≈ freqs_ atol=SELFatol
        @test hcat(damps[i,:]...)' ≈ damps_ atol=SELFatol
    end
end

@testset "Flutter analysis the Helios flying-wing in free flight with payload as the varying parameter" begin
    include("examples/heliosFlutterPRange.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosFlutterPRange", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosFlutterPRange", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Flutter analysis the Helios flying-wing in free flight with airspeed as the varying parameter" begin
    include("examples/heliosFlutterURange.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosFlutterURange", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosFlutterURange", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Trim analysis of the Helios flying-wing" begin
    include("examples/heliosTrim.jl")
    # Self-comparison
    trimAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosTrim", "trimAoA.txt"))
    trimThrust_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosTrim", "trimThrust.txt"))
    trimδ_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosTrim", "trimDelta.txt"))
    @test trimAoA ≈ trimAoA_ atol=SELFatol
    @test trimThrust ≈ trimThrust_ atol=SELFatol
    @test trimδ ≈ trimδ_ atol=SELFatol
end

@testset "Flutter analysis of the wing of the Helios flying-wing" begin
    include("examples/heliosWingFlutter.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosWingFlutter", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "heliosWingFlutter", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

@testset "Static analysis of a hinged beam subjected to a distributed load" begin
    include("examples/hingedBeam.jl")
    # Self-comparison
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedBeam", "u3.txt"))
    p2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedBeam", "p2.txt"))
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedBeam", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedBeam", "M2.txt"))
    @test u3 ≈ u3_ atol=SELFatol
    @test p2 ≈ p2_ atol=SELFatol
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end

@testset "Static analysis of a hinged beam subjected to a distributed load and a rotational spring" begin
    include("examples/hingedSpringedBeam.jl")
    # Self-comparison
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedSpringedBeam", "u3.txt"))
    p2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedSpringedBeam", "p2.txt"))
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedSpringedBeam", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedSpringedBeam", "M2.txt"))
    @test u3 ≈ u3_ atol=SELFatol
    @test p2 ≈ p2_ atol=SELFatol
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end

@testset "Static analysis of a T-frame hinged at the connection subjected to a distributed load" begin
    include("examples/hingedTFrame.jl")
    # Self-comparison
    u1_beam1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedTFrame", "u1_beam1.txt"))
    u1_beam2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedTFrame", "u1_beam2.txt"))
    u3_beam1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedTFrame", "u3_beam1.txt"))
    u3_beam2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedTFrame", "u3_beam2.txt"))
    F3_beam1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedTFrame", "F3_beam1.txt"))
    M2_beam1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "hingedTFrame", "M2_beam1.txt"))
    @test u1_beam1 ≈ u1_beam1_ atol=SELFatol
    @test u1_beam2 ≈ u1_beam2_ atol=SELFatol
    @test u3_beam1 ≈ u3_beam1 atol=SELFatol
    @test u3_beam2 ≈ u3_beam2_ atol=SELFatol
    @test F3_beam1 ≈ F3_beam1_ atol=SELFatol
    @test M2_beam1 ≈ M2_beam1_ atol=SELFatol
end

@testset "Dynamic analysis of the free response of a beam subjected to initial displacement and velocity profiles" begin
    include("examples/initialDispAndVelBeam.jl")
    # Analytical comparison
    @test u3_quarter ≈ u3_quarter_analytic rtol=1e-2
    @test V3_quarter ≈ V3_quarter_analytic rtol=1e-2
    @test Vdot3_quarter ≈ Vdot3_quarter_analytic rtol=1e-2
    @test θ2_root ≈ θ2_root_analytic rtol=1e-2
    @test Ω2_mid ≈ Ω2_mid_analytic rtol=1e-2
    # @test Ωdot2_mid ≈ Ωdot2_mid_analytic atol=20
    # Self-comparison
    u3_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialDispAndVelBeam", "u3_quarter.txt"))
    V3_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialDispAndVelBeam", "V3_quarter.txt"))
    Vdot3_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialDispAndVelBeam", "Vdot3_quarter.txt"))
    θ2_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialDispAndVelBeam", "theta2_root.txt"))
    Ω2_mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialDispAndVelBeam", "Omega2_mid.txt"))
    Ωdot2_mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialDispAndVelBeam", "Omegadot2_mid.txt"))
    @test u3_quarter ≈ u3_quarter_ atol=SELFatol
    @test V3_quarter ≈ V3_quarter_ atol=SELFatol
    @test Vdot3_quarter ≈ Vdot3_quarter_ atol=SELFatol
    @test θ2_root ≈ θ2_root_ atol=SELFatol
    @test Ω2_mid ≈ Ω2_mid_ atol=SELFatol
    @test Ωdot2_mid ≈ Ωdot2_mid_ atol=SELFatol
end

@testset "Dynamic analysis of the free response of a beam subjected to an initial displacement profile" begin
    include("examples/initialDisplacementBeam.jl")
    # Analytical comparison
    @test u3_quarter ≈ u3_quarter_analytic rtol=1e-2
    @test V3_quarter ≈ V3_quarter_analytic rtol=1e-2
    @test Vdot3_quarter ≈ Vdot3_quarter_analytic rtol=1e-2
    @test θ2_root ≈ θ2_root_analytic rtol=1e-2
    @test Ω2_mid ≈ Ω2_mid_analytic rtol=1e-2
    @test Ωdot2_mid ≈ Ωdot2_mid_analytic rtol=1e-2
    # Self-comparison
    u3_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialDisplacementBeam", "u3_quarter.txt"))
    V3_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialDisplacementBeam", "V3_quarter.txt"))
    Vdot3_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialDisplacementBeam", "Vdot3_quarter.txt"))
    θ2_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialDisplacementBeam", "theta2_root.txt"))
    Ω2_mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialDisplacementBeam", "Omega2_mid.txt"))
    Ωdot2_mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialDisplacementBeam", "Omegadot2_mid.txt"))
    @test u3_quarter ≈ u3_quarter_ atol=SELFatol
    @test V3_quarter ≈ V3_quarter_ atol=SELFatol
    @test Vdot3_quarter ≈ Vdot3_quarter_ atol=SELFatol
    @test θ2_root ≈ θ2_root_ atol=SELFatol
    @test Ω2_mid ≈ Ω2_mid_ atol=SELFatol
    @test Ωdot2_mid ≈ Ωdot2_mid_ atol=SELFatol
end

@testset "Dynamic analysis of the free response of a beam subjected to an initial velocity profile" begin
    include("examples/initialVelocityBeam.jl")
    # Analytical comparison
    @test u3_quarter ≈ u3_quarter_analytic rtol=1e-2
    @test V3_quarter ≈ V3_quarter_analytic rtol=1e-2
    @test Vdot3_quarter ≈ Vdot3_quarter_analytic rtol=1e-2
    @test θ2_root ≈ θ2_root_analytic rtol=1e-2
    @test Ω2_mid ≈ Ω2_mid_analytic rtol=1e-2
    @test Ωdot2_mid ≈ Ωdot2_mid_analytic rtol=1e-2
    # Self-comparison
    u3_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialVelocityBeam", "u3_quarter.txt"))
    V3_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialVelocityBeam", "V3_quarter.txt"))
    Vdot3_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialVelocityBeam", "Vdot3_quarter.txt"))
    θ2_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialVelocityBeam", "theta2_root.txt"))
    Ω2_mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialVelocityBeam", "Omega2_mid.txt"))
    Ωdot2_mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "initialVelocityBeam", "Omegadot2_mid.txt"))
    @test u3_quarter ≈ u3_quarter_ atol=SELFatol
    @test V3_quarter ≈ V3_quarter_ atol=SELFatol
    @test Vdot3_quarter ≈ Vdot3_quarter_ atol=SELFatol
    @test θ2_root ≈ θ2_root_ atol=SELFatol
    @test Ω2_mid ≈ Ω2_mid_ atol=SELFatol
    @test Ωdot2_mid ≈ Ωdot2_mid_ atol=SELFatol
end

@testset "Dynamic analysis of joined beams under load" begin
    include("examples/joinedBeams.jl")
    # Self-comparison
    u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "joinedBeams", "u1.txt"))
    u2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "joinedBeams", "u2.txt"))
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "joinedBeams", "u3.txt"))
    @test u1 ≈ u1_ atol=SELFatol
    @test u2 ≈ u2_ atol=SELFatol
    @test u3 ≈ u3_ atol=SELFatol
end

@testset "Static analysis of the Lee frame (a right-angled frame) with a dead load" begin
    include("examples/LeeFrameDeadLoad.jl")
    # Self-comparison
    u1_atForce_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "LeeFrameDeadLoad", "u1_atForce.txt"))
    u3_atForce_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "LeeFrameDeadLoad", "u3_atForce.txt"))
    @test u1_atForce ≈ u1_atForce_ atol=SELFatol
    @test u3_atForce ≈ u3_atForce_ atol=SELFatol
end

@testset "Static analysis of the Lee frame (a right-angled frame) with a follower load" begin
    include("examples/LeeFrameFollowerLoad.jl")
    # Self-comparison
    u1_atForce_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "LeeFrameFollowerLoad", "u1_atForce.txt"))
    u3_atForce_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "LeeFrameFollowerLoad", "u3_atForce.txt"))
    @test u1_atForce ≈ u1_atForce_ atol=SELFatol
    @test u3_atForce ≈ u3_atForce_ atol=SELFatol
end

@testset "Trim analysis (reaction loads check) of a simply-supported beam loaded at the middle" begin
    include("examples/midLoadedBeamTrim.jl")
    # Self-comparison
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "midLoadedBeamTrim", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "midLoadedBeamTrim", "M2.txt"))
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end

@testset "Dynamic analysis of a pinned robot arm driven by a couple moment" begin
    include("examples/momentDrivenRobotArm.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "momentDrivenRobotArm", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "momentDrivenRobotArm", "u3_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol
    @test u3_tip ≈ u3_tip_ atol=SELFatol
end

@testset "One-minus-cosine gust response of an airfoil section at several pitch angles" begin
    include("examples/OMCgustTests.jl")
    # Self-comparison
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
            for (k,testCase) in enumerate(1:6)
                if typeof(aeroSolver) == QuasiSteady
                    aeroSolverName = "QS"
                elseif typeof(aeroSolver) == Indicial
                    aeroSolverName = "Indicial"
                elseif typeof(aeroSolver) == Inflow
                    aeroSolverName = "Inflow"
                elseif typeof(aeroSolver) == BLi
                    aeroSolverName = "BLi"
                elseif typeof(aeroSolver) == BLo
                    aeroSolverName = "BLo"
                end
                gustSolverName = gustLoadsSolver.indicialFunctionName
                Δcl_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "OMCgustTests", string("OMCgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dcl",i,j,k,".txt")))
                @test Δcl[i,j,k] ≈ Δcl_ atol=SELFatol
            end
        end
    end
end

@testset "Modal analysis of the Pazy wing with flared folding wing tip (FFWT)" begin
    include("examples/PazyFFWTeigen.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTeigen", "freqs.txt"))
    damps_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyFFWTeigen", "damps.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
    @test hcat(damps...)' ≈ damps_ atol=SELFatol
end

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

@testset "Static analysis of the pure bending test of the Pazy wing" begin
    include("examples/PazyWingBendingTest.jl")
    # Self-comparison
    tip_OOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingBendingTest", "tip_OOP.txt"))
    tip_twist_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingBendingTest", "tip_twist.txt"))
    @test tip_OOP ≈ tip_OOP_ atol=SELFatol
    @test tip_twist ≈ tip_twist_ atol=SELFatol
end

@testset "Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over time" begin
    include("examples/PazyWingContinuous1DGust.jl")
    # Self-comparison
    tipAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DGust", "tipAoA.txt"))
    tipOOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DGust", "tipOOP.txt"))
    tqSpan_cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DGust", "tqSpan_cn.txt"))
    tqSpan_cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DGust", "tqSpan_cm.txt"))
    tqSpan_ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DGust", "tqSpan_ct.txt"))
    @test tipAoA ≈ tipAoA_ atol=SELFatol
    @test tipOOP ≈ tipOOP_ atol=SELFatol
    @test tqSpan_cn ≈ tqSpan_cn_ atol=SELFatol
    @test tqSpan_cm ≈ tqSpan_cm_ atol=SELFatol
    @test tqSpan_ct ≈ tqSpan_ct_ atol=SELFatol
end

@testset "Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over space" begin
    include("examples/PazyWingContinuous1DSpaceGust.jl")
    # Self-comparison
    tipAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DSpaceGust", "tipAoA.txt"))
    tipOOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DSpaceGust", "tipOOP.txt"))
    tqSpan_cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DSpaceGust", "tqSpan_cn.txt"))
    tqSpan_cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DSpaceGust", "tqSpan_cm.txt"))
    tqSpan_ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous1DSpaceGust", "tqSpan_ct.txt"))
    @test tipAoA ≈ tipAoA_ atol=SELFatol
    @test tipOOP ≈ tipOOP_ atol=SELFatol
    @test tqSpan_cn ≈ tqSpan_cn_ atol=SELFatol
    @test tqSpan_cm ≈ tqSpan_cm_ atol=SELFatol
    @test tqSpan_ct ≈ tqSpan_ct_ atol=SELFatol
end

@testset "Dynamic analysis of the Pazy wing encountering a continuous, 2-dimensional gust" begin
    include("examples/PazyWingContinuous2DSpaceGust.jl")
    # Self-comparison
    tipAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous2DSpaceGust", "tipAoA.txt"))
    tipOOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous2DSpaceGust", "tipOOP.txt"))
    tqSpan_cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous2DSpaceGust", "tqSpan_cn.txt"))
    tqSpan_cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous2DSpaceGust", "tqSpan_cm.txt"))
    tqSpan_ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingContinuous2DSpaceGust", "tqSpan_ct.txt"))
    @test tipAoA ≈ tipAoA_ atol=SELFatol
    @test tipOOP ≈ tipOOP_ atol=SELFatol
    @test tqSpan_cn ≈ tqSpan_cn_ atol=SELFatol
    @test tqSpan_cm ≈ tqSpan_cm_ atol=SELFatol
    @test tqSpan_ct ≈ tqSpan_ct_ atol=SELFatol
end

@testset "Dynamic analysis of the Pazy wing encountering a DARPA gust" begin
    include("examples/PazyWingDARPAGust.jl")
    # Self-comparison
    tipAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingDARPAGust", "tipAoA.txt"))
    tipOOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingDARPAGust", "tipOOP.txt"))
    tqSpan_cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingDARPAGust", "tqSpan_cn.txt"))
    tqSpan_cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingDARPAGust", "tqSpan_cm.txt"))
    tqSpan_ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingDARPAGust", "tqSpan_ct.txt"))
    @test tipAoA ≈ tipAoA_ atol=SELFatol
    @test tipOOP ≈ tipOOP_ atol=SELFatol
    @test tqSpan_cn ≈ tqSpan_cn_ atol=SELFatol
    @test tqSpan_cm ≈ tqSpan_cm_ atol=SELFatol
    @test tqSpan_ct ≈ tqSpan_ct_ atol=SELFatol
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

@testset "Modal analyses of the Pazy wing in horizontal and vertical positions" begin
    include("examples/PazyWingModal.jl")
    # Reference comparison
    @test freqs[1][1:4] ≈ modalFreqsGVT[1,1:4] rtol=5e-2
    @test freqs[2][1:4] ≈ modalFreqsGVT[2,1:4] rtol=5e-2
    @test freqs[1] ≈ modalFreqsUMNAST[1,:] rtol=2e-2
    @test freqs[2] ≈ modalFreqsUMNAST[2,:] rtol=2e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingModal", "freqs.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
end

@testset "Dynamic analysis of the Pazy wing encountering a one-minus-cosine gust" begin
    include("examples/PazyWingOMCGust.jl")
    # Self-comparison
    tipAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingOMCGust", "tipAoA.txt"))
    tipOOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingOMCGust", "tipOOP.txt"))
    tqSpan_cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingOMCGust", "tqSpan_cn.txt"))
    tqSpan_cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingOMCGust", "tqSpan_cm.txt"))
    tqSpan_ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingOMCGust", "tqSpan_ct.txt"))
    @test tipAoA ≈ tipAoA_ atol=SELFatol
    @test tipOOP ≈ tipOOP_ atol=SELFatol
    @test tqSpan_cn ≈ tqSpan_cn_ atol=SELFatol
    @test tqSpan_cm ≈ tqSpan_cm_ atol=SELFatol
    @test tqSpan_ct ≈ tqSpan_ct_ atol=SELFatol
end

@testset "Steady analysis of the Pazy wing with varying root pitch angle" begin
    include("examples/PazyWingPitchRange.jl")
    # Self-comparison
    tip_AoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingPitchRange", "tip_AoA.txt"))
    tip_OOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingPitchRange", "tip_OOP.txt"))
    tip_IP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingPitchRange", "tip_IP.txt"))
    tip_twist_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingPitchRange", "tip_twist.txt"))
    @test filter(!isnan,tip_AoA) ≈ filter(!isnan,tip_AoA_) atol=SELFatol
    @test tip_OOP ≈ tip_OOP_ atol=SELFatol
    @test tip_IP ≈ tip_IP_ atol=SELFatol
    @test tip_twist ≈ tip_twist_ atol=SELFatol
end

@testset "Dynamic analysis of the Pazy wing with a tip impulse force" begin
    include("examples/PazyWingTipImpulse.jl")
    # Self-comparison
    tipAoA_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTipImpulse", "tipAoA.txt"))
    tipOOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTipImpulse", "tipOOP.txt"))
    tqSpan_cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTipImpulse", "tqSpan_cn.txt"))
    tqSpan_cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTipImpulse", "tqSpan_cm.txt"))
    tqSpan_ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTipImpulse", "tqSpan_ct.txt"))
    @test tipAoA ≈ tipAoA_ atol=SELFatol
    @test tipOOP ≈ tipOOP_ atol=SELFatol
    @test tqSpan_cn ≈ tqSpan_cn_ atol=SELFatol
    @test tqSpan_cm ≈ tqSpan_cm_ atol=SELFatol
    @test tqSpan_ct ≈ tqSpan_ct_ atol=SELFatol
end

@testset "Static analysis of the coupled torsion-bending test of the Pazy wing" begin
    include("examples/PazyWingTorsionTest.jl")
    # Self-comparison
    tip_OOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTorsionTest", "tip_OOP.txt"))
    tip_twist_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTorsionTest", "tip_twist.txt"))
    @test tip_OOP ≈ tip_OOP_ atol=SELFatol
    @test tip_twist ≈ tip_twist_ atol=SELFatol
end

@testset "Dynamic analysis of a pendulum released from rest" begin
    include("examples/pendulum.jl")
    # Reference (analytical) comparison
    @test u1_tip ≈ u1_tip_analytical rtol=5e-2
    @test u3_tip ≈ u3_tip_analytical rtol=5e-2
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "pendulum", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "pendulum", "u3_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol
    @test u3_tip ≈ u3_tip_ atol=SELFatol
end

@testset "Static analysis of a mid-loaded arch pinned at one end and clamped at the other" begin
    include("examples/pinnedClampedArch.jl")
    # Reference (analytical) comparison
    @test -u1_atForce[end]/R ≈ u1Ref[1,45] rtol=2e-2
    @test -u3_atForce[end]/R ≈ u3Ref[1,57] rtol=2e-2
    # Self-comparison
    u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "pinnedClampedArch", "u1.txt"))
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "pinnedClampedArch", "u3.txt"))
    @test u1_atForce ≈ u1_ atol=SELFatol
    @test u3_atForce ≈ u3_ atol=SELFatol
end

@testset "Modal analysis of a beam pinned at one end and transversely springed at the other" begin
    include("examples/pinnedSpringedBeamEigen.jl")
    # Analytical comparison
    @test freqs ≈ freqsRef rtol=1e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "pinnedSpringedBeamEigen", "freqs.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
end

@testset "Static analysis of a right-angled frame under a tip transverse follower force" begin
    include("examples/rightAngledFrame.jl")
    # Self-comparison
    tip_u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rightAngledFrame", "tip_u1.txt"))
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rightAngledFrame", "tip_u3.txt"))
    tip_angle_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rightAngledFrame", "tip_angle.txt"))
    @test tip_u1 ≈ tip_u1_ atol=SELFatol
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
    @test tip_angle ≈ tip_angle_ atol=SELFatol
end

@testset "Static analysis of a right-angled frame under a 'buckling' load" begin
    include("examples/rightAngledFrameBuckling.jl")
    # Self-comparison
    tip_u2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rightAngledFrameBuckling", "tip_u2.txt"))
    @test tip_u2 ≈ tip_u2_ atol=SELFatol
end

@testset "Trim analysis (reaction loads check) of a right-angled frame" begin
    include("examples/rightAngledFrameTrim.jl")
    # Analytical comparison
    @test balanceVerticalForce ≈ balanceVerticalForceAnalytical rtol=1e-4
    @test balanceHorizontalForce ≈ balanceHorizontalForceAnalytical rtol=1e-4
    # Self-comparison
    balanceHorizontalForce_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rightAngledFrameTrim", "balanceHorizontalForce.txt"))[1]
    balanceVerticalForce_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rightAngledFrameTrim", "balanceVerticalForce.txt"))[1]
    @test balanceHorizontalForce ≈ balanceHorizontalForce_ atol=SELFatol
    @test balanceVerticalForce ≈ balanceVerticalForce_ atol=SELFatol
end

@testset "Dynamic analysis of a clamped beam with root sinusoidal oscillation (ω = 9 Hz ≈ 2nd bending mode)" begin
    include("examples/rootExcitationBeam1.jl")
    # Self-comparison
    u3b_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam1", "u3b_root.txt"))
    u3b_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam1", "u3b_tip.txt"))
    V3_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam1", "V3_root.txt"))
    V3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam1", "V3_tip.txt"))
    @test u3b_root ≈ u3b_root_ atol=SELFatol
    @test u3b_tip ≈ u3b_tip_ atol=SELFatol
    @test V3_root ≈ V3_root_ atol=SELFatol
    @test V3_tip ≈ V3_tip_ atol=SELFatol
end

@testset "Dynamic analysis of a clamped beam with root sinusoidal oscillation (ω = 32 Hz ≈ 3rd bending mode)" begin
    include("examples/rootExcitationBeam2.jl")
    # Self-comparison
    u3b_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam2", "u3b_root.txt"))
    u3b_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam2", "u3b_tip.txt"))
    V3_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam2", "V3_root.txt"))
    V3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam2", "V3_tip.txt"))
    @test u3b_root ≈ u3b_root_ atol=SELFatol
    @test u3b_tip ≈ u3b_tip_ atol=SELFatol
    @test V3_root ≈ V3_root_ atol=SELFatol
    @test V3_tip ≈ V3_tip_ atol=SELFatol
end

@testset "Dynamic analysis of a rotary shaft with specified rotation" begin
    include("examples/rotaryShaft.jl")
    # Analytical comparison
    @test pNum ≈ p.(t) rtol = 1e-2
    @test pdotNum ≈ pdot.(t) rtol = 1e-2
    @test ΩNum ≈ θdot.(t) rtol = 1e-2
    @test ΩdotNum ≈ θddot.(t) rtol = 1e-2
    @test MNum ≈ -θddot.(t)*ρ*Is rtol = 0.1
    # Self-comparison
    pNum_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotaryShaft", "pNum.txt"))
    pdotNum_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotaryShaft", "pdotNum.txt"))
    ΩNum_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotaryShaft", "OmegaNum.txt"))
    ΩdotNum_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotaryShaft", "OmegadotNum.txt"))
    MNum_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotaryShaft", "MNum.txt"))
    @test pNum ≈ pNum_ atol=SELFatol
    @test pdotNum ≈ pdotNum_ atol=SELFatol
    @test ΩNum ≈ ΩNum_ atol=SELFatol
    @test ΩdotNum ≈ ΩdotNum_ atol=SELFatol
    @test MNum ≈ MNum_ atol=SELFatol
end

@testset "Dynamic analysis of an articulated robot arm driven by specified rotation" begin
    include("examples/rotationDrivenArticulatedRobotArm.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotationDrivenArticulatedRobotArm", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotationDrivenArticulatedRobotArm", "u3_tip.txt"))
    u1_hinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotationDrivenArticulatedRobotArm", "u1_hinge.txt"))
    u3_hinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotationDrivenArticulatedRobotArm", "u3_hinge.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol
    @test u3_tip ≈ u3_tip_ atol=SELFatol
    @test u1_hinge ≈ u1_hinge_ atol=SELFatol
    @test u3_hinge ≈ u3_hinge_ atol=SELFatol
end

@testset "Dynamic analysis of a robot arm driven by specified rotation" begin
    include("examples/rotationDrivenRobotArm.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotationDrivenRobotArm", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotationDrivenRobotArm", "u3_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol
    @test u3_tip ≈ u3_tip_ atol=SELFatol
end

@testset "Sharp-edged gust response of an airfoil section at several pitch angles" begin
    include("examples/SEgustTests.jl")
    # Self-comparison
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
            for (k,testCase) in enumerate(1:3)
                if typeof(aeroSolver) == QuasiSteady
                    aeroSolverName = "QS"
                elseif typeof(aeroSolver) == Indicial
                    aeroSolverName = "Indicial"
                elseif typeof(aeroSolver) == Inflow
                    aeroSolverName = "Inflow"
                elseif typeof(aeroSolver) == BLi
                    aeroSolverName = "BLi"
                elseif typeof(aeroSolver) == BLo
                    aeroSolverName = "BLo"
                end
                gustSolverName = gustLoadsSolver.indicialFunctionName
                Δcl_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SEgustTests", string("SEgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dcl",i,j,k,".txt")))
                @test Δcl[i,j,k] ≈ Δcl_ atol=SELFatol
            end
        end
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

@testset "Linear flutter analysis of the sixteen-meter-wing" begin
    include("examples/SMWLinearFlutter.jl")
    # Reference comparison
    @test flutterSpeed ≈ flutterSpeedRef rtol=1e-2
    @test flutterFreq ≈ flutterFreqRef rtol=2e-2
    # Self-comparison
    flutterSpeed_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWLinearFlutter", "flutterSpeed.txt"))[1]
    flutterFreq_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWLinearFlutter", "flutterFreq.txt"))[1]
    @test flutterSpeed ≈ flutterSpeed_ atol=SELFatol
    @test flutterFreq ≈ flutterFreq_ atol=SELFatol
end

@testset "Modal analysis of the sixteen-meter-wing" begin
    include("examples/SMWModal.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=2e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWModal", "freqs.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
end

@testset "Steady aeroelastic analysis of the sixteen-meter-wing" begin
    include("examples/SMWSteady.jl")
    # Self-comparison
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWSteady", "tip_u3.txt"))
    tip_twist_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWSteady", "tip_twist.txt"))
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
    @test tip_twist ≈ tip_twist_ atol=SELFatol
end

@testset "Dynamic analysis of the spin-up maneuver of a robot arm" begin
    include("examples/spinupRobotArm.jl")
    # Reference comparison
    @test θ3_root[end] ≈ θ(t[end]) rtol=1e-2
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "spinupRobotArm", "u1_tip.txt"))
    u2_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "spinupRobotArm", "u2_tip.txt"))
    θ3_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "spinupRobotArm", "theta3_root.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol
    @test u2_tip ≈ u2_tip_ atol=SELFatol
    @test θ3_root ≈ θ3_root_ atol=SELFatol
end

@testset "Static analysis of a L-frame with a doubly-attached spring and tip load" begin
    include("examples/springedLFrame.jl")
    # Self-comparison
    u3_b_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "springedLFrame", "u3_b.txt"))
    @test u3_b ≈ u3_b_ atol=SELFatol
end

@testset "Eigen-analysis of a straight rotor" begin
    include("examples/straightRotor.jl")
    # Reference comparison
    @test numFreqs[end][[1,2,5,6]]/(2π) ≈ expFreqs[:,end] rtol=5e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "straightRotor", "freqs.txt"))
    @test hcat(numFreqs...)' ≈ freqs_ atol=SELFatol
end

@testset "Eigen-analysis of a swept-tip rotor" begin
    include("examples/sweptTipRotor.jl")
    # Reference comparison (3 first bending frequencies @ ω = 750 rpm and tipAngle = 45 deg)
    @test numFreqs[end,end][[1,3,4]]/(2π) ≈ [expFreqs1[end]; expFreqs2[end]; expFreqs3[end]] rtol=5e-2
    # Self-comparison (frequencies @ ω = 750 rpm and tipAngle = 45 deg)
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "sweptTipRotor", "freqs.txt"))
    @test numFreqs[end,end] ≈ freqs_ atol=SELFatol
end

@testset "Static analysis of a semi-circular arch with tangential follower force" begin
    include("examples/tangentiallyForcedArch.jl")
    # Reference comparison
    @test -tip_u1[end]/R ≈ u1Ref[1,end] rtol=5e-2
    @test tip_u3[end]/R ≈ u3Ref[1,end] rtol=5e-2
    @test -tip_angle[end]/(π/2) ≈ θRef[1,end] rtol=5e-2
    # Self-comparison
    tip_u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tangentiallyForcedArch", "tip_u1.txt"))
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tangentiallyForcedArch", "tip_u3.txt"))
    tip_angle_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tangentiallyForcedArch", "tip_angle.txt"))
    @test tip_u1 ≈ tip_u1_ atol=SELFatol
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
    @test tip_angle ≈ tip_angle_ atol=SELFatol
end

@testset "Modal analysis of a cantilevered tapered beam" begin
    include("examples/taperedBeamEigen.jl")
    # Analytical comparison
    @test freqs ≈ freqsRef rtol=1e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "taperedBeamEigen", "freqs.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
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

@testset "Modal analysis of the Tang&Dowell wing at varying pitch angles" begin
    include("examples/TDWingPitchRange.jl")
    # Reference comparison
    @test -tip_u3[1] ≈ u3_exp[2,1] rtol=5e-2
    @test -tip_u2[end] ≈ u2_exp[2,end] rtol=5e-2
    @test -tip_twist[div(length(θRange)-1,2)] ≈ th_exp[2,2] rtol=0.1
    @test freqs[1][[1,2,4]] ≈ freqs_exp[2:4,1] rtol=0.1
    # Self-comparison
    tip_u2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "TDWingPitchRange", "tip_u2.txt"))
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "TDWingPitchRange", "tip_u3.txt"))
    tip_twist_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "TDWingPitchRange", "tip_twist.txt"))
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "TDWingPitchRange", "freqs.txt"))
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
    @test tip_u2 ≈ tip_u2_ atol=SELFatol
    @test tip_twist ≈ tip_twist_ atol=SELFatol
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
end

@testset "Dynamic analysis of an airfoil section sinusoidally surging (facing a time-varying freestream)" begin
    include("examples/timeVaryingFreestream.jl")
    # Loop normalized airspeed amplitude
    for i in eachindex(λᵤRange)
        # Reference comparison
        @test cn[i][end]/(2π*θ) ≈ cnCFD[i][2,end] atol=0.3
        @test cm[i][end] ≈ cmCFD[i][2,end] atol=1e-2
        # Analytical comparison
        @test Vdot2[i][end] ≈ Vdot2Analytical[i][end] rtol=1e-2
        @test Vdot3[i][end] ≈ Vdot3Analytical[i][end] rtol=1e-2
    end
    # Self-comparison (for the last λᵤ)
    cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "timeVaryingFreestream", "cn.txt"))
    cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "timeVaryingFreestream", "cm.txt"))
    @test cn[end] ≈ cn_ atol=SELFatol
    @test cm[end] ≈ cm_ atol=SELFatol
end

@testset "Dynamic analysis of an airfoil section sinusoidally pitching and surging (facing a time-varying freestream)" begin
    include("examples/timeVaryingFreestreamAndPitch.jl")
    # Loop normalized airspeed amplitude
    for i in eachindex(λᵤRange)
        # index at 3/4 of the cycle
        indexTqcycle = rangeLastCycle[i][round(Int,length(rangeLastCycle[i])*3/4)]
        # Reference comparison
        @test cn[i][end]/cn_qs(t[i][end]) ≈ cnCFD[i][2,end] atol=0.5
        @test cm[i][end] ≈ cmCFD[i][2,end] atol=1e-2
        # Analytical comparison
        @test Vdot2[i][end] ≈ Vdot2Analytical[i][end] rtol=1e-2
        @test Vdot3[i][end] ≈ Vdot3Analytical[i][end] rtol=1e-2
        @test Ω1[i][end] ≈ Ω1Analytical(t[i][end]) rtol=1e-2
        @test Ωdot1[i][indexTqcycle] ≈ Ωdot1Analytical(t[i][indexTqcycle]) rtol=0.1
    end
    # Self-comparison (for the last λᵤ)
    cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "timeVaryingFreestreamAndPitch", "cn.txt"))
    cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "timeVaryingFreestreamAndPitch", "cm.txt"))
    @test cn[end] ≈ cn_ atol=SELFatol
    @test cm[end] ≈ cm_ atol=SELFatol
end

@testset "Static analysis cantilever with tip follower transverse force" begin
    include("examples/tipFollowerForceCantilever.jl")
    # Reference comparison
    @test -tip_u1[end]/L ≈ u1Ref[2,end]/L rtol=1e-2
    @test tip_u3[end]/L ≈ u3Ref[2,end]/L rtol=1e-2
    # Self-comparison
    tip_u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipFollowerForceCantilever", "tip_u1.txt"))
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipFollowerForceCantilever", "tip_u3.txt"))
    @test tip_u1 ≈ tip_u1_ atol=SELFatol
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
end

@testset "Static analysis of a cantilever with tip follower transverse force (force split over 2 BCs)" begin
    include("examples/tipFollowerForceCantilever2.jl")
    # Reference comparison
    @test -tip_u1[end]/L ≈ u1Ref[2,end]/L rtol=1e-2
    @test tip_u3[end]/L ≈ u3Ref[2,end]/L rtol=1e-2
    # Self-comparison
    tip_u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipFollowerForceCantilever2", "tip_u1.txt"))
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipFollowerForceCantilever2", "tip_u3.txt"))
    @test tip_u1 ≈ tip_u1_ atol=SELFatol
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
end

@testset "Trim analysis of a cantilever with tip force" begin
    include("examples/tipLoadedCantileverTrim.jl")
    # Analytical comparison
    @test u3[end] ≈ u3_analytical(L) rtol=1e-3
    @test F3[1] ≈ F3_analytical(0) rtol=1e-3
    @test M2[1] ≈ M2_analytical(0) rtol=1e-3
    # Self-comparison
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipLoadedCantileverTrim", "u3.txt"))
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipLoadedCantileverTrim", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipLoadedCantileverTrim", "M2.txt"))
    @test u3 ≈ u3_ atol=SELFatol
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end

@testset "Static analysis of a cantilever with tip moment" begin
    include("examples/tipMomentCantilever.jl")
    # Analytical comparison
    @test tip_u1[end]/L ≈ -1 atol=1e-3
    @test tip_u3[end]/L ≈ 0 atol=1e-4
    @test tip_angle[end]/π ≈ 2 atol=1e-4
    # Self-comparison
    tip_u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipMomentCantilever", "tip_u1.txt"))
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipMomentCantilever", "tip_u3.txt"))
    tip_angle_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipMomentCantilever", "tip_angle.txt"))
    @test tip_u1 ≈ tip_u1_ atol=SELFatol
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
    @test tip_angle ≈ tip_angle_ atol=SELFatol
end

@testset "Dynamic analysis of a pendulum with tip mass" begin
    include("examples/tipPendulum.jl")
    # Analytical comparison
    @test u1_tip[end]/L ≈ u1_tip_analytical[end]/L atol=1e-3
    @test u3_tip[end]/L ≈ u3_tip_analytical[end]/L atol=1e-3
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipPendulum", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipPendulum", "u3_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol
    @test u3_tip ≈ u3_tip_ atol=SELFatol
end

@testset "Dynamic analysis of a cantilever with tip sinusoidal force" begin
    include("examples/tipSineLoadedCantilever.jl")
    # Self-comparison
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipSineLoadedCantilever", "u3_tip.txt"))
    F3_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipSineLoadedCantilever", "F3_root.txt"))
    M2_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipSineLoadedCantilever", "M2_root.txt"))
    @test u3_tip ≈ u3_tip_ atol=SELFatol
    @test F3_root ≈ F3_root_ atol=SELFatol
    @test M2_root ≈ M2_root_ atol=SELFatol
end

@testset "Static analysis of a semi-circular arch with transverse follower force" begin
    include("examples/transverselyForcedArch.jl")
    # Reference comparison
    @test -tip_u1[end]/R ≈ u1Ref[1,end] rtol=2e-2
    @test -tip_u3[end]/R ≈ u3Ref[1,end] atol=5e-2
    @test -tip_angle[end]/(π/2) ≈ θRef[1,end] rtol=5e-2
    # Self-comparison
    tip_u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "transverselyForcedArch", "tip_u1.txt"))
    tip_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "transverselyForcedArch", "tip_u3.txt"))
    tip_angle_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "transverselyForcedArch", "tip_angle.txt"))
    @test tip_u1 ≈ tip_u1_ atol=SELFatol
    @test tip_u3 ≈ tip_u3_ atol=SELFatol
    @test tip_angle ≈ tip_angle_ atol=SELFatol
end

@testset "Static analysis of a beam with triangular distributed load" begin
    include("examples/triangleLoadBeam.jl")
    # Self-comparison
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "triangleLoadBeam", "u3.txt"))
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "triangleLoadBeam", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "triangleLoadBeam", "M2.txt"))
    @test u3 ≈ u3_ atol=SELFatol
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end

@testset "Eigen analysis of a 2-story frame" begin
    include("examples/twoStoryFrame.jl")
    # Reference comparison
    @test freqs[[1,4]]/(2π) ≈ refFreqs rtol=1e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "twoStoryFrame", "freqs.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
end

println("Finished tests")
# Eigenvalue structural (modal) problems

@testset "Modal analysis of the axial vibration of a beam under clamped-clamped boundary conditions" begin
    include("examples/beamAxialVibrationCC.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=5e-3
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamAxialVibrationCC", "freqs.txt"))
    u1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamAxialVibrationCC", "u1_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(u1_modeShapes...)' ≈ u1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of the axial vibration of a beam under clamped-free boundary conditions" begin
    include("examples/beamAxialVibrationCF.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=5e-3
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamAxialVibrationCF", "freqs.txt"))
    u1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamAxialVibrationCF", "u1_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(u1_modeShapes...)' ≈ u1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of the axial vibration of a beam under free-free boundary conditions" begin
    include("examples/beamAxialVibrationFF.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=5e-3
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamAxialVibrationFF", "freqs.txt"))
    u1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamAxialVibrationFF", "u1_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(u1_modeShapes...)' ≈ u1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of the bending vibration of a beam under clamped-clamped boundary conditions" begin
    include("examples/beamBendingVibrationCC.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=1e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCC", "freqs.txt"))
    u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCC", "u3_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of the bending vibration of a beam under clamped-free boundary conditions" begin
    include("examples/beamBendingVibrationCF.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=1e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCF", "freqs.txt"))
    u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCF", "u3_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of the bending vibration of a beam under clamped-pinned boundary conditions" begin
    include("examples/beamBendingVibrationCP.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=1e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCP", "freqs.txt"))
    u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCP", "u3_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of the bending vibration of a beam under clamped-sliding boundary conditions" begin
    include("examples/beamBendingVibrationCS.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=1e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCS", "freqs.txt"))
    u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationCS", "u3_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of the bending vibration of a beam under free-free boundary conditions" begin
    include("examples/beamBendingVibrationFF.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=1e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationFF", "freqs.txt"))
    u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationFF", "u3_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of the bending vibration of a beam under pinned-pinned boundary conditions" begin
    include("examples/beamBendingVibrationPP.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=1e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationPP", "freqs.txt"))
    u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamBendingVibrationPP", "u3_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of the torsional vibration of a beam under clamped-clamped boundary conditions" begin
    include("examples/beamTorsionalVibrationCC.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=5e-3
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamTorsionalVibrationCC", "freqs.txt"))
    p1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamTorsionalVibrationCC", "p1_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(p1_modeShapes...)' ≈ p1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of the torsional vibration of a beam under clamped-free boundary conditions" begin
    include("examples/beamTorsionalVibrationCF.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=5e-3
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamTorsionalVibrationCF", "freqs.txt"))
    p1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamTorsionalVibrationCF", "p1_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(p1_modeShapes...)' ≈ p1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of the torsional vibration of a beam under free-free boundary conditions" begin
    include("examples/beamTorsionalVibrationFF.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=5e-3
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamTorsionalVibrationFF", "freqs.txt"))
    p1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "beamTorsionalVibrationFF", "p1_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(p1_modeShapes...)' ≈ p1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of a cantilever beam with a tip axial inertia" begin
    include("examples/cantileverWithTipAxialMassEigen.jl")
    # Analytical comparison
    @test freqsNorm ≈ freqsNormAnalytical rtol=2e-2
    # Self-comparison
    freqsNorm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipAxialMassEigen", "freqsNorm.txt"))
    u1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipAxialMassEigen", "u1_modeShapes.txt"))
    @test freqsNorm ≈ freqsNorm_ atol=SELFatol
    # @test hcat(u1_modeShapes...)' ≈ u1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of a cantilever beam with a tip axial spring" begin
    include("examples/cantileverWithTipAxialSpringEigen.jl")
    # Analytical comparison
    @test freqsNorm ≈ freqsNormAnalytical rtol=5e-3
    # Self-comparison
    freqsNorm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipAxialSpringEigen", "freqsNorm.txt"))
    u1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipAxialSpringEigen", "u1_modeShapes.txt"))
    @test freqsNorm ≈ freqsNorm_ atol=SELFatol
    # @test hcat(u1_modeShapes...)' ≈ u1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of a cantilever beam with a tip spring in bending" begin
    include("examples/cantileverWithTipSpringEigen.jl")
    # Analytical comparison
    @test freqsNorm ≈ freqsNormRef rtol=1e-2
    # Self-comparison
    freqsNorm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipSpringEigen", "freqsNorm.txt"))
    u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipSpringEigen", "u3_modeShapes.txt"))
    @test freqsNorm ≈ freqsNorm_ atol=SELFatol
    # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of a cantilever beam with a tip torsional inertia" begin
    include("examples/cantileverWithTipTorsionalInertiaEigen.jl")
    # Analytical comparison
    @test freqsNorm ≈ freqsNormAnalytical rtol=2e-2
    # Self-comparison
    freqsNorm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipTorsionalInertiaEigen", "freqsNorm.txt"))
    p1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipTorsionalInertiaEigen", "p1_modeShapes.txt"))
    @test freqsNorm ≈ freqsNorm_ atol=SELFatol
    # @test hcat(p1_modeShapes...)' ≈ p1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of a cantilever beam with a tip torsional spring" begin
    include("examples/cantileverWithTipTorsionalSpringEigen.jl")
    # Analytical comparison
    @test freqsNorm ≈ freqsNormAnalytical rtol=1e-2
    # Self-comparison
    freqsNorm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipTorsionalSpringEigen", "freqsNorm.txt"))
    p1_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "cantileverWithTipTorsionalSpringEigen", "p1_modeShapes.txt"))
    @test freqsNorm ≈ freqsNorm_ atol=SELFatol
    # @test hcat(p1_modeShapes...)' ≈ p1_modeShapes_ atol=SELFatol
end

@testset "Modal analysis a beam clamped at one end, simply-supported at the other and with a tip inertia" begin
    include("examples/clampedSSBeamWIthTipInertia.jl")
    # Reference (analytical) comparison
    @test freqs[1] ≈ freqRef rtol=2e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedSSBeamWIthTipInertia", "freqs.txt"))
    u3_modeShapes_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedSSBeamWIthTipInertia", "u3_modeShapes.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
    # @test hcat(u3_modeShapes...)' ≈ u3_modeShapes_ atol=SELFatol
end

@testset "Modal analysis of the baseline Healy free FFWT wing without gravity" begin
    include("examples/HealyBaselineFFWTModalFree.jl")
    # Reference comparison
    @test freqs[1] ≈ freqsRef[1] atol=1e-2
    @test freqs[order[2:9]] ≈ freqsRef[2:end] rtol=6e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "HealyBaselineFFWTModalFree", "freqs.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
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

@testset "Modal analysis of a beam pinned at one end and transversely springed at the other" begin
    include("examples/pinnedSpringedBeamEigen.jl")
    # Analytical comparison
    @test freqs ≈ freqsRef rtol=1e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "pinnedSpringedBeamEigen", "freqs.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
end

@testset "Modal analysis of the sixteen-meter-wing" begin
    include("examples/SMWModal.jl")
    # Analytical comparison
    @test freqs ≈ freqsAnalytical rtol=2e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SMWModal", "freqs.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
end

@testset "Modal analysis of a straight rotor under varying angular velocities" begin
    include("examples/straightRotor.jl")
    # Reference comparison
    @test numFreqs[end][[1,2,5,6]]/(2π) ≈ expFreqs[:,end] rtol=5e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "straightRotor", "freqs.txt"))
    @test hcat(numFreqs...)' ≈ freqs_ atol=SELFatol
end

@testset "Modal analysis of the undeformed swept Pazy wing" begin
    include("examples/sweptPazyModal.jl")
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "sweptPazyModal", "freqs.txt"))
    @test hcat(freqs...)' ≈ freqs_ atol=SELFatol
end

@testset "Modal analysis of a swept-tip rotor under varying angular velocities" begin
    include("examples/sweptTipRotor.jl")
    # Reference comparison (3 first bending frequencies @ ω = 750 rpm and tipAngle = 45 deg)
    @test numFreqs[end,end][[1,3,4]]/(2π) ≈ [expFreqs1[end]; expFreqs2[end]; expFreqs3[end]] rtol=5e-2
    # Self-comparison (frequencies @ ω = 750 rpm and tipAngle = 45 deg)
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "sweptTipRotor", "freqs.txt"))
    @test numFreqs[end,end] ≈ freqs_ atol=SELFatol
end

@testset "Modal analysis of a cantilevered tapered beam" begin
    include("examples/taperedBeamEigen.jl")
    # Analytical comparison
    @test freqs ≈ freqsRef rtol=1e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "taperedBeamEigen", "freqs.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
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

@testset "Modal analysis of a 2-story frame" begin
    include("examples/twoStoryFrame.jl")
    # Reference comparison
    @test freqs[[1,4]]/(2π) ≈ refFreqs rtol=1e-2
    # Self-comparison
    freqs_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "twoStoryFrame", "freqs.txt"))
    @test freqs ≈ freqs_ atol=SELFatol
end
# Static structural problems

@testset "Static analysis of an arch under a dead pressure load" begin
    include("examples/archUnderDeadPressure.jl")
    # Self-comparison
    mid_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "archUnderDeadPressure", "mid_u3.txt"))
    @test mid_u3 ≈ mid_u3_ atol=SELFatol
end

@testset "Static analysis of an arch under a follower pressure load" begin
    include("examples/archUnderFollowerPressure.jl")
    # Self-comparison
    mid_u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "archUnderFollowerPressure", "mid_u3.txt"))
    @test mid_u3 ≈ mid_u3_ atol=SELFatol
end

@testset "Static analysis of a cantilever beam with an axial spring attached between its middle and tip nodes, subjected to an axial tip force" begin
    include("examples/axialDoublyAttachedSpringCantilever.jl")
    # Reference comparison
    @test u1[end] ≈ uTipAnalytical rtol=1e-10
    @test norm(spring.Fs) ≈ FsAnalytical rtol=1e-10
    # Self-comparison
    u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "axialDoublyAttachedSpringCantilever", "u1.txt"))
    F1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "axialDoublyAttachedSpringCantilever", "F1.txt"))
    @test u1 ≈ u1_ atol=SELFatol
    @test F1 ≈ F1_ atol=SELFatol
end

@testset "Static analysis of a cantilever beam with an axial spring attached between its middle and tip nodes, subjected to an axial tip displacement" begin
    include("examples/axialDoublyAttachedSpringCantilever2.jl")
    # Reference comparison
    @test F1[1] ≈ FrootAnalytical rtol=1e-10
    @test norm(spring.Fs) ≈ FsAnalytical rtol=1e-10
    # Self-comparison
    u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "axialDoublyAttachedSpringCantilever2", "u1.txt"))
    F1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "axialDoublyAttachedSpringCantilever2", "F1.txt"))
    @test u1 ≈ u1_ atol=SELFatol
    @test F1 ≈ F1_ atol=SELFatol
end

@testset "Static analysis of a biclamped, hinged beam under distributed and concentrated loads" begin
    include("examples/biclampedHingedBeam.jl")
    # Reference comparison
    @test u3Mid ≈ u3MidRef rtol=1e-2
    @test p2Left ≈ p2LeftRef rtol=2e-2
    @test p2Right ≈ p2RightRef rtol=1e-2
    # Self-comparison
    u3Mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedHingedBeam", "u3Mid.txt"))[1]
    p2Left_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedHingedBeam", "p2Left.txt"))[1]
    p2Right_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedHingedBeam", "p2Right.txt"))[1]
    @test u3Mid ≈ u3Mid_ atol=SELFatol
    @test p2Left ≈ p2Left_ atol=SELFatol
    @test p2Right ≈ p2Right_ atol=SELFatol
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

@testset "Static analysis of a clamped beam, with a free flared hinge at the middle, under distributed and concentrated loads" begin
    include("examples/clampedFlaredHingedBeam.jl")
    # Self-comparison
    u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedFlaredHingedBeam", "u1.txt"))
    u2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedFlaredHingedBeam", "u2.txt"))
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedFlaredHingedBeam", "u3.txt"))
    p1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedFlaredHingedBeam", "p1.txt"))
    p2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedFlaredHingedBeam", "p2.txt"))
    p3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedFlaredHingedBeam", "p3.txt"))
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedFlaredHingedBeam", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedFlaredHingedBeam", "M2.txt"))
    ϕHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedFlaredHingedBeam", "phiHinge.txt"))[1]
    @test u1 ≈ u1_ atol=SELFatol
    @test u2 ≈ u2_ atol=SELFatol
    @test u3 ≈ u3_ atol=SELFatol
    @test p1 ≈ p1_ atol=SELFatol
    @test p2 ≈ p2_ atol=SELFatol
    @test p3 ≈ p3_ atol=SELFatol
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
    @test ϕHinge ≈ ϕHinge_ atol=SELFatol
end

@testset "Static analysis of a clamped beam, hinged at the middle with imposed hinge angle, under distributed and concentrated loads" begin
    include("examples/clampedHingedBeam.jl")
    # Self-comparison
    u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeam", "u1.txt"))
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeam", "u3.txt"))
    p2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeam", "p2.txt"))
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeam", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeam", "M2.txt"))
    @test u1 ≈ u1_ atol=SELFatol
    @test u3 ≈ u3_ atol=SELFatol
    @test p2 ≈ p2_ atol=SELFatol
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end

@testset "Static analysis of a clamped beam rotated in 3D space, hinged at the middle with imposed hinge angle, under distributed and concentrated loads" begin
    include("examples/clampedHingedBeamRotated.jl")
    # Self-comparison
    u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeamRotated", "u1.txt"))
    u2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeamRotated", "u2.txt"))
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeamRotated", "u3.txt"))
    p2_b_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeamRotated", "p2_b.txt"))
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeamRotated", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeamRotated", "M2.txt"))
    @test u1 ≈ u1_ atol=SELFatol
    @test u2 ≈ u2_ atol=SELFatol
    @test u3 ≈ u3_ atol=SELFatol
    @test p2_b ≈ p2_b_ atol=SELFatol
    @test F3 ≈ F3_ atol=SELFatol
    @test M2 ≈ M2_ atol=SELFatol
end

@testset "Static analysis of a clamped beam, hinged at the middle, under distributed loads, with varying spring stiffnesses around the hinge" begin
    include("examples/clampedHingedBeamSpringRange.jl")
    # Self-comparison
    u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeamSpringRange", "u1.txt"))
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeamSpringRange", "u3.txt"))
    p2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeamSpringRange", "p2.txt"))
    F3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeamSpringRange", "F3.txt"))
    M2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "clampedHingedBeamSpringRange", "M2.txt"))
    @test isapprox(vcat(u1'...), u1_, atol=SELFatol, nans=true)
    @test isapprox(vcat(u3'...), u3_, atol=SELFatol, nans=true)
    @test isapprox(vcat(p2'...), p2_, atol=SELFatol, nans=true)
    @test isapprox(vcat(F3'...), F3_, atol=SELFatol, nans=true)
    @test isapprox(vcat(M2'...), M2_, atol=SELFatol, nans=true)
end

@testset "Static analysis of a clamped beam with a coasting flared hinge at the middle, under distributed loads" begin
    include("examples/coastingFoldingWingtip.jl")
    # Self-comparison
    pHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "coastingFoldingWingtip", "pHinge.txt"))
    @test pHinge ≈ pHinge_ atol=SELFatol
end

@testset "Static analysis of a clamped beam with a driven flared hinge at the middle, under distributed loads" begin
    include("examples/drivenFoldingWingtip.jl")
    # Self-comparison
    pHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "drivenFoldingWingtip", "pHinge.txt"))
    @test pHinge ≈ pHinge_ atol=SELFatol
end

@testset "Static analysis of a clamped beam with a driven, springed, flared hinge at the middle" begin
    include("examples/drivenSpringedHingedBeam.jl")
    # Self-comparison
    pHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "drivenSpringedHingedBeam", "pHinge.txt"))
    @test pHinge ≈ pHinge_ atol=SELFatol
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

@testset "Static analysis of a hinged beam subjected to a distributed load" begin
    include("examples/hingedBeam.jl")
    # Reference comparison
    @test u3Mid ≈ u3MidRef rtol=1e-2
    @test p2Left ≈ p2LeftRef rtol=2e-2
    @test p2Right ≈ p2RightRef rtol=1e-2
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

@testset "Static analysis of the pure bending test of the Pazy wing" begin
    include("examples/PazyWingBendingTest.jl")
    # Self-comparison
    tip_OOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingBendingTest", "tip_OOP.txt"))
    @test tip_OOP ≈ tip_OOP_ atol=SELFatol
end

@testset "Static analysis of the coupled torsion-bending test of the Pazy wing" begin
    include("examples/PazyWingTorsionTest.jl")
    # Self-comparison
    tip_OOP_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTorsionTest", "tip_OOP.txt"))
    tip_twist_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "PazyWingTorsionTest", "tip_twist.txt"))
    @test tip_OOP ≈ tip_OOP_ atol=SELFatol
    @test tip_twist ≈ tip_twist_ atol=SELFatol
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

@testset "Static analysis of a L-frame with a doubly-attached spring and tip load" begin
    include("examples/springedLFrame.jl")
    # Self-comparison
    u3_b_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "springedLFrame", "u3_b.txt"))
    @test u3_b ≈ u3_b_ atol=SELFatol
end

@testset "Static analysis of a swept-back clamped beam with a driven flared folding wingtip" begin
    include("examples/sweptBackFFWTWing.jl")
    # Self-comparison
    ϕ_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "sweptBackFFWTWing", "phi.txt"))[1]
    hingeBalanceM_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "sweptBackFFWTWing", "hingeBalanceM.txt"))[1]
    pHinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "sweptBackFFWTWing", "pHinge.txt"))
    @test hingeBalanceM ≈ hingeBalanceM_ atol=SELFatol
    @test pHinge ≈ pHinge_ atol=SELFatol
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

@testset "Static analysis of a cantilever beam with a torsional spring attached between its middle and tip nodes, subjected to a torsional tip moment" begin
    include("examples/twistDoublyAttachedSpringCantilever.jl")
    # Reference comparison
    @test p1[end] ≈ pTipAnalytical rtol=2e-4
    @test norm(spring.Ms) ≈ MsAnalytical rtol=1e-4
    # Self-comparison
    p1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "twistDoublyAttachedSpringCantilever", "p1.txt"))
    M1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "twistDoublyAttachedSpringCantilever", "M1.txt"))
    @test p1 ≈ p1_ atol=SELFatol
    @test M1 ≈ M1_ atol=SELFatol
end
GC.gc()
sleep(1)
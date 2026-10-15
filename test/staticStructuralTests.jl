@testset "Static analysis of an arch under a dead pressure load" begin
    include_example("examples/staticStructural/archUnderDeadPressure.jl")
    @test matches_test_data(mid_u3, "archUnderDeadPressure", "mid_u3.txt")
end

@testset "Static analysis of an arch under a follower pressure load" begin
    include_example("examples/staticStructural/archUnderFollowerPressure.jl")
    @test matches_test_data(mid_u3, "archUnderFollowerPressure", "mid_u3.txt")
end

@testset "Static analysis of a cantilever beam with an axial spring attached between its middle and tip nodes, subjected to an axial tip force" begin
    include_example("examples/staticStructural/axialDoublyAttachedSpringCantilever.jl")
    @test matches_test_data(u1, "axialDoublyAttachedSpringCantilever", "u1.txt")
    @test matches_test_data(F1, "axialDoublyAttachedSpringCantilever", "F1.txt")
end

@testset "Static analysis of a cantilever beam with an axial spring attached between its middle and tip nodes, subjected to an axial tip displacement" begin
    include_example("examples/staticStructural/axialDoublyAttachedSpringCantilever2.jl")
    @test matches_test_data(u1, "axialDoublyAttachedSpringCantilever2", "u1.txt")
    @test matches_test_data(F1, "axialDoublyAttachedSpringCantilever2", "F1.txt")
end

@testset "Static analysis of a biclamped, hinged beam under distributed and concentrated loads" begin
    include_example("examples/staticStructural/biclampedHingedBeam.jl")
    @test matches_test_data(u3Mid, "biclampedHingedBeam", "u3Mid.txt")
    @test matches_test_data(p2Left, "biclampedHingedBeam", "p2Left.txt")
    @test matches_test_data(p2Right, "biclampedHingedBeam", "p2Right.txt")
end

@testset "Static analysis of a cantilever beam bending under self weight" begin
    include_example("examples/staticStructural/cantileverUnderSelfWeight.jl")
    @test matches_test_data(u1, "cantileverUnderSelfWeight", "u1.txt")
    @test matches_test_data(u3, "cantileverUnderSelfWeight", "u3.txt")
    @test matches_test_data(F3, "cantileverUnderSelfWeight", "F3.txt")
    @test matches_test_data(M2, "cantileverUnderSelfWeight", "M2.txt")
end

@testset "Static analysis of a cantilever beam with a tip spring in bending" begin
    include_example("examples/staticStructural/cantileverWithTipSpring.jl")
    @test matches_test_data(u3, "cantileverWithTipSpring", "u3.txt")
    @test matches_test_data(F3, "cantileverWithTipSpring", "F3.txt")
    @test matches_test_data(M2, "cantileverWithTipSpring", "M2.txt")
end

@testset "Static analysis of a clamped beam, with a free flared hinge at the middle, under distributed and concentrated loads" begin
    include_example("examples/staticStructural/clampedFlaredHingedBeam.jl")
    @test matches_test_data(u1, "clampedFlaredHingedBeam", "u1.txt")
    @test matches_test_data(u2, "clampedFlaredHingedBeam", "u2.txt")
    @test matches_test_data(u3, "clampedFlaredHingedBeam", "u3.txt")
    @test matches_test_data(p1, "clampedFlaredHingedBeam", "p1.txt")
    @test matches_test_data(p2, "clampedFlaredHingedBeam", "p2.txt")
    @test matches_test_data(p3, "clampedFlaredHingedBeam", "p3.txt")
    @test matches_test_data(F3, "clampedFlaredHingedBeam", "F3.txt")
    @test matches_test_data(M2, "clampedFlaredHingedBeam", "M2.txt")
    @test matches_test_data(ϕHinge, "clampedFlaredHingedBeam", "phiHinge.txt")
end

@testset "Static analysis of a clamped beam, hinged at the middle with imposed hinge angle, under distributed and concentrated loads" begin
    include_example("examples/staticStructural/clampedHingedBeam.jl")
    @test matches_test_data(u1, "clampedHingedBeam", "u1.txt")
    @test matches_test_data(u3, "clampedHingedBeam", "u3.txt")
    @test matches_test_data(p2, "clampedHingedBeam", "p2.txt")
    @test matches_test_data(F3, "clampedHingedBeam", "F3.txt")
    @test matches_test_data(M2, "clampedHingedBeam", "M2.txt")
end

@testset "Static analysis of a clamped beam rotated in 3D space, hinged at the middle with imposed hinge angle, under distributed and concentrated loads" begin
    include_example("examples/staticStructural/clampedHingedBeamRotated.jl")
    @test matches_test_data(u1, "clampedHingedBeamRotated", "u1.txt")
    @test matches_test_data(u2, "clampedHingedBeamRotated", "u2.txt")
    @test matches_test_data(u3, "clampedHingedBeamRotated", "u3.txt")
    @test matches_test_data(p2_b, "clampedHingedBeamRotated", "p2_b.txt")
    @test matches_test_data(F3, "clampedHingedBeamRotated", "F3.txt")
    @test matches_test_data(M2, "clampedHingedBeamRotated", "M2.txt")
end

@testset "Static analysis of a clamped beam, hinged at the middle, under distributed loads, with varying spring stiffnesses around the hinge" begin
    include_example("examples/staticStructural/clampedHingedBeamSpringRange.jl")
    @test matches_test_data(u1, "clampedHingedBeamSpringRange", "u1.txt")
    @test matches_test_data(u3, "clampedHingedBeamSpringRange", "u3.txt")
    @test matches_test_data(p2, "clampedHingedBeamSpringRange", "p2.txt")
    @test matches_test_data(F3, "clampedHingedBeamSpringRange", "F3.txt")
    @test matches_test_data(M2, "clampedHingedBeamSpringRange", "M2.txt")
end

@testset "Static analysis of a clamped beam with a coasting flared hinge at the middle, under distributed loads" begin
    include_example("examples/staticStructural/coastingFoldingWingtip.jl")
    @test matches_test_data(pHinge, "coastingFoldingWingtip", "pHinge.txt")
end

@testset "Static analysis of composite laminates subjected to tip loads" begin
    include_example("examples/staticStructural/compositeCantileverMD.jl")
    @test matches_test_data(u1_500mm, "compositeCantileverMD", "u1_500mm.txt")
    @test matches_test_data(u2_500mm, "compositeCantileverMD", "u2_500mm.txt")
    @test matches_test_data(u3_500mm, "compositeCantileverMD", "u3_500mm.txt")
end

@testset "Static analysis of a curved cantilever subjected to a tip dead force" begin
    include_example("examples/staticStructural/curvedCantileverDeadLoad.jl")
    @test matches_test_data(tip_u1, "curvedCantileverDeadLoad", "tip_u1.txt")
    @test matches_test_data(tip_u2, "curvedCantileverDeadLoad", "tip_u2.txt")
    @test matches_test_data(tip_u3, "curvedCantileverDeadLoad", "tip_u3.txt")
end

@testset "Static analysis of a curved cantilever subjected to a tip follower force" begin
    include_example("examples/staticStructural/curvedCantileverStaticFollower.jl")
    @test matches_test_data(tip_u1, "curvedCantileverStaticFollower", "tip_u1.txt")
    @test matches_test_data(tip_u2, "curvedCantileverStaticFollower", "tip_u2.txt")
    @test matches_test_data(tip_u3, "curvedCantileverStaticFollower", "tip_u3.txt")
end

@testset "Static analysis of a cantilever with distributed follower force" begin
    include_example("examples/staticStructural/distributedLoadCantilever.jl")
    @test matches_test_data(tip_u1, "distributedLoadCantilever", "tip_u1.txt")
    @test matches_test_data(tip_u3, "distributedLoadCantilever", "tip_u3.txt")
    @test matches_test_data(tip_angle, "distributedLoadCantilever", "tip_angle.txt")
end

@testset "Static analysis of a clamped beam with a driven flared hinge at the middle, under distributed loads" begin
    include_example("examples/staticStructural/drivenFoldingWingtip.jl")
    @test matches_test_data(pHinge, "drivenFoldingWingtip", "pHinge.txt")
end

@testset "Static analysis of a clamped beam with a driven, springed, flared hinge at the middle" begin
    include_example("examples/staticStructural/drivenSpringedHingedBeam.jl")
    @test matches_test_data(pHinge, "drivenSpringedHingedBeam", "pHinge.txt")
end

@testset "Static analysis of a hinged beam subjected to a distributed load" begin
    include_example("examples/staticStructural/hingedBeam.jl")
    @test matches_test_data(u3, "hingedBeam", "u3.txt")
    @test matches_test_data(p2, "hingedBeam", "p2.txt")
    @test matches_test_data(F3, "hingedBeam", "F3.txt")
    @test matches_test_data(M2, "hingedBeam", "M2.txt")
end

@testset "Static analysis of a hinged beam subjected to a distributed load and a rotational spring" begin
    include_example("examples/staticStructural/hingedSpringedBeam.jl")
    @test matches_test_data(u3, "hingedSpringedBeam", "u3.txt")
    @test matches_test_data(p2, "hingedSpringedBeam", "p2.txt")
    @test matches_test_data(F3, "hingedSpringedBeam", "F3.txt")
    @test matches_test_data(M2, "hingedSpringedBeam", "M2.txt")
end

@testset "Static analysis of a T-frame hinged at the connection subjected to a distributed load" begin
    include_example("examples/staticStructural/hingedTFrame.jl")
    @test matches_test_data(u1_beam1, "hingedTFrame", "u1_beam1.txt")
    @test matches_test_data(u1_beam2, "hingedTFrame", "u1_beam2.txt")
    @test matches_test_data(u3_beam1, "hingedTFrame", "u3_beam1.txt")
    @test matches_test_data(u3_beam2, "hingedTFrame", "u3_beam2.txt")
    @test matches_test_data(F3_beam1, "hingedTFrame", "F3_beam1.txt")
    @test matches_test_data(M2_beam1, "hingedTFrame", "M2_beam1.txt")
end

@testset "Static analysis of the Lee frame (a right-angled frame) with a dead load" begin
    include_example("examples/staticStructural/LeeFrameDeadLoad.jl")
    @test matches_test_data(u1_atForce, "LeeFrameDeadLoad", "u1_atForce.txt")
    @test matches_test_data(u3_atForce, "LeeFrameDeadLoad", "u3_atForce.txt")
end

@testset "Static analysis of the Lee frame (a right-angled frame) with a follower load" begin
    include_example("examples/staticStructural/LeeFrameFollowerLoad.jl")
    @test matches_test_data(u1_atForce, "LeeFrameFollowerLoad", "u1_atForce.txt")
    @test matches_test_data(u3_atForce, "LeeFrameFollowerLoad", "u3_atForce.txt")
end

@testset "Static analysis of the pure bending test of the Pazy wing" begin
    include_example("examples/staticStructural/PazyWingBendingTest.jl")
    @test matches_test_data(tip_OOP, "PazyWingBendingTest", "tip_OOP.txt")
end

@testset "Static analysis of the coupled torsion-bending test of the Pazy wing" begin
    include_example("examples/staticStructural/PazyWingTorsionTest.jl")
    @test matches_test_data(tip_OOP, "PazyWingTorsionTest", "tip_OOP.txt")
    @test matches_test_data(tip_twist, "PazyWingTorsionTest", "tip_twist.txt")
end

@testset "Static analysis of a mid-loaded arch pinned at one end and clamped at the other" begin
    include_example("examples/staticStructural/pinnedClampedArch.jl")
    @test matches_test_data(u1_atForce, "pinnedClampedArch", "u1.txt")
    @test matches_test_data(u3_atForce, "pinnedClampedArch", "u3.txt")
end

@testset "Static analysis of a right-angled frame under a tip transverse follower force" begin
    include_example("examples/staticStructural/rightAngledFrame.jl")
    @test matches_test_data(tip_u1, "rightAngledFrame", "tip_u1.txt")
    @test matches_test_data(tip_u3, "rightAngledFrame", "tip_u3.txt")
    @test matches_test_data(tip_angle, "rightAngledFrame", "tip_angle.txt")
end

@testset "Static analysis of a right-angled frame under a 'buckling' load" begin
    include_example("examples/staticStructural/rightAngledFrameBuckling.jl")
    @test matches_test_data(tip_u2, "rightAngledFrameBuckling", "tip_u2.txt")
end

@testset "Static analysis of a L-frame with a doubly-attached spring and tip load" begin
    include_example("examples/staticStructural/springedLFrame.jl")
    @test matches_test_data(u3_b, "springedLFrame", "u3_b.txt")
end

@testset "Static analysis of a swept-back clamped beam with a driven flared folding wingtip" begin
    include_example("examples/staticStructural/sweptBackFFWTWing.jl")
    @test matches_test_data(ϕHinge, "sweptBackFFWTWing", "phi.txt")
    @test matches_test_data(hingeMoment, "sweptBackFFWTWing", "hingeMoment.txt")
end

@testset "Static analysis of a semi-circular arch with tangential follower force" begin
    include_example("examples/staticStructural/tangentiallyForcedArch.jl")
    @test matches_test_data(tip_u1, "tangentiallyForcedArch", "tip_u1.txt")
    @test matches_test_data(tip_u3, "tangentiallyForcedArch", "tip_u3.txt")
    @test matches_test_data(tip_angle, "tangentiallyForcedArch", "tip_angle.txt")
end

@testset "Static analysis cantilever with tip follower transverse force" begin
    include_example("examples/staticStructural/tipFollowerForceCantilever.jl")
    @test matches_test_data(tip_u1, "tipFollowerForceCantilever", "tip_u1.txt")
    @test matches_test_data(tip_u3, "tipFollowerForceCantilever", "tip_u3.txt")
end

@testset "Static analysis of a cantilever with tip follower transverse force (force split over 2 BCs)" begin
    include_example("examples/staticStructural/tipFollowerForceCantilever2.jl")
    @test matches_test_data(tip_u1, "tipFollowerForceCantilever2", "tip_u1.txt")
    @test matches_test_data(tip_u3, "tipFollowerForceCantilever2", "tip_u3.txt")
end

@testset "Static analysis of a cantilever with tip moment" begin
    include_example("examples/staticStructural/tipMomentCantilever.jl")
    @test matches_test_data(tip_u1, "tipMomentCantilever", "tip_u1.txt")
    @test matches_test_data(tip_u3, "tipMomentCantilever", "tip_u3.txt")
    @test matches_test_data(tip_angle, "tipMomentCantilever", "tip_angle.txt")
end

@testset "Static analysis of a semi-circular arch with transverse follower force" begin
    include_example("examples/staticStructural/transverselyForcedArch.jl")
    @test matches_test_data(tip_u1, "transverselyForcedArch", "tip_u1.txt")
    @test matches_test_data(tip_u3, "transverselyForcedArch", "tip_u3.txt")
    @test matches_test_data(tip_angle, "transverselyForcedArch", "tip_angle.txt")
end

@testset "Static analysis of a beam with triangular distributed load" begin
    include_example("examples/staticStructural/triangleLoadBeam.jl")
    @test matches_test_data(u3, "triangleLoadBeam", "u3.txt")
    @test matches_test_data(F3, "triangleLoadBeam", "F3.txt")
    @test matches_test_data(M2, "triangleLoadBeam", "M2.txt")
end

@testset "Static analysis of a cantilever beam with a torsional spring attached between its middle and tip nodes, subjected to a torsional tip moment" begin
    include_example("examples/staticStructural/twistDoublyAttachedSpringCantilever.jl")
    @test matches_test_data(p1, "twistDoublyAttachedSpringCantilever", "p1.txt")
    @test matches_test_data(M1, "twistDoublyAttachedSpringCantilever", "M1.txt")
end

@testset "Static analysis of a beam with a central hinge driven by opposite end moments" begin
    include_example("examples/staticStructural/drivenSpringedHingedBeam2.jl")
    # Self-comparison
    @test matches_test_data(u3, "drivenSpringedHingedBeam2", "u3.txt")
    @test matches_test_data(p2, "drivenSpringedHingedBeam2", "p2.txt")
    @test matches_test_data(M2, "drivenSpringedHingedBeam2", "M2.txt")
end

@testset "Static analysis of a swept cantilever under tip transverse loading" begin
    include_example("examples/staticStructural/sweptBeamStaticLoading.jl")
    # Self-comparison
    @test matches_test_data(tip_u3, "sweptBeamStaticLoading", "tip_u3.txt")
    @test matches_test_data(tip_twist, "sweptBeamStaticLoading", "tip_twist.txt")
end
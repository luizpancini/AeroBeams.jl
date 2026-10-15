# Static analysis of an arch under a dead pressure load
include("../examples/staticStructural/archUnderDeadPressure.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/archUnderDeadPressure"))
writedlm("test/newTestDataGenerators/archUnderDeadPressure/mid_u3.txt", test_data_matrix(mid_u3))

# Static analysis of an arch under a follower pressure load
include("../examples/staticStructural/archUnderFollowerPressure.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/archUnderFollowerPressure"))
writedlm("test/newTestDataGenerators/archUnderFollowerPressure/mid_u3.txt", test_data_matrix(mid_u3))

# Static analysis of a cantilever beam with an axial spring attached between its middle and tip nodes, subjected to an axial tip force
include("../examples/staticStructural/axialDoublyAttachedSpringCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/axialDoublyAttachedSpringCantilever"))
writedlm("test/newTestDataGenerators/axialDoublyAttachedSpringCantilever/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/axialDoublyAttachedSpringCantilever/F1.txt", test_data_matrix(F1))

# Static analysis of a cantilever beam with an axial spring attached between its middle and tip nodes, subjected to an axial tip displacement
include("../examples/staticStructural/axialDoublyAttachedSpringCantilever2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/axialDoublyAttachedSpringCantilever2"))
writedlm("test/newTestDataGenerators/axialDoublyAttachedSpringCantilever2/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/axialDoublyAttachedSpringCantilever2/F1.txt", test_data_matrix(F1))

# Static analysis of a biclamped, hinged beam under distributed and concentrated loads
include("../examples/staticStructural/biclampedHingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/biclampedHingedBeam"))
writedlm("test/newTestDataGenerators/biclampedHingedBeam/u3Mid.txt", test_data_matrix(u3Mid))
writedlm("test/newTestDataGenerators/biclampedHingedBeam/p2Left.txt", test_data_matrix(p2Left))
writedlm("test/newTestDataGenerators/biclampedHingedBeam/p2Right.txt", test_data_matrix(p2Right))

# Static analysis of a cantilever beam bending under self weight
include("../examples/staticStructural/cantileverUnderSelfWeight.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverUnderSelfWeight"))
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/M2.txt", test_data_matrix(M2))

# Static analysis of a cantilever beam with a tip spring in bending
include("../examples/staticStructural/cantileverWithTipSpring.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipSpring"))
writedlm("test/newTestDataGenerators/cantileverWithTipSpring/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/cantileverWithTipSpring/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/cantileverWithTipSpring/M2.txt", test_data_matrix(M2))

# Static analysis of a clamped beam, with a free flared hinge at the middle, under distributed and concentrated loads
include("../examples/staticStructural/clampedFlaredHingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedFlaredHingedBeam"))
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/u2.txt", test_data_matrix(u2))
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/p1.txt", test_data_matrix(p1))
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/p2.txt", test_data_matrix(p2))
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/p3.txt", test_data_matrix(p3))
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/M2.txt", test_data_matrix(M2))
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/phiHinge.txt", test_data_matrix(ϕHinge))

# Static analysis of a clamped beam, hinged at the middle with imposed hinge angle, under distributed and concentrated loads
include("../examples/staticStructural/clampedHingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedHingedBeam"))
writedlm("test/newTestDataGenerators/clampedHingedBeam/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/clampedHingedBeam/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/clampedHingedBeam/p2.txt", test_data_matrix(p2))
writedlm("test/newTestDataGenerators/clampedHingedBeam/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/clampedHingedBeam/M2.txt", test_data_matrix(M2))

# Static analysis of a clamped beam rotated in 3D space, hinged at the middle with imposed hinge angle, under distributed and concentrated loads
include("../examples/staticStructural/clampedHingedBeamRotated.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedHingedBeamRotated"))
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/u2.txt", test_data_matrix(u2))
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/p2_b.txt", test_data_matrix(p2_b))
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/M2.txt", test_data_matrix(M2))

# Static analysis of a clamped beam, hinged at the middle, under distributed loads, with varying spring stiffnesses around the hinge
include("../examples/staticStructural/clampedHingedBeamSpringRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedHingedBeamSpringRange"))
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/u1.txt", test_data_matrix(u1))
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/p2.txt", test_data_matrix(p2))
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/M2.txt", test_data_matrix(M2))

# Static analysis of a clamped beam with a coasting flared hinge at the middle, under distributed loads
include("../examples/staticStructural/coastingFoldingWingtip.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/coastingFoldingWingtip"))
writedlm("test/newTestDataGenerators/coastingFoldingWingtip/pHinge.txt", test_data_matrix(pHinge))

# Static analysis of a clamped beam with a driven flared hinge at the middle, under distributed loads
include("../examples/staticStructural/drivenFoldingWingtip.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/drivenFoldingWingtip"))
writedlm("test/newTestDataGenerators/drivenFoldingWingtip/pHinge.txt", test_data_matrix(pHinge))

# Static analysis of a clamped beam with a driven, springed, flared hinge at the middle
include("../examples/staticStructural/drivenSpringedHingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/drivenSpringedHingedBeam"))
writedlm("test/newTestDataGenerators/drivenSpringedHingedBeam/pHinge.txt", test_data_matrix(pHinge))

# Static analysis of composite laminates subjected to tip loads
include("../examples/staticStructural/compositeCantileverMD.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/compositeCantileverMD"))
writedlm(string("test/newTestDataGenerators/compositeCantileverMD/u1_500mm.txt"), test_data_matrix(u1_500mm))
writedlm(string("test/newTestDataGenerators/compositeCantileverMD/u2_500mm.txt"), test_data_matrix(u2_500mm))
writedlm(string("test/newTestDataGenerators/compositeCantileverMD/u3_500mm.txt"), test_data_matrix(u3_500mm))

# Static analysis of a curved cantilever subjected to a tip dead force
include("../examples/staticStructural/curvedCantileverDeadLoad.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/curvedCantileverDeadLoad"))
writedlm("test/newTestDataGenerators/curvedCantileverDeadLoad/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/curvedCantileverDeadLoad/tip_u2.txt", test_data_matrix(tip_u2))
writedlm("test/newTestDataGenerators/curvedCantileverDeadLoad/tip_u3.txt", test_data_matrix(tip_u3))

# Static analysis of a curved cantilever subjected to a tip follower force
include("../examples/staticStructural/curvedCantileverStaticFollower.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/curvedCantileverStaticFollower"))
writedlm("test/newTestDataGenerators/curvedCantileverStaticFollower/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/curvedCantileverStaticFollower/tip_u2.txt", test_data_matrix(tip_u2))
writedlm("test/newTestDataGenerators/curvedCantileverStaticFollower/tip_u3.txt", test_data_matrix(tip_u3))

# Static analysis of a cantilever with distributed follower force
include("../examples/staticStructural/distributedLoadCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/distributedLoadCantilever"))
writedlm("test/newTestDataGenerators/distributedLoadCantilever/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/distributedLoadCantilever/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/distributedLoadCantilever/tip_angle.txt", test_data_matrix(tip_angle))

# Static analysis of a hinged beam subjected to a distributed load
include("../examples/staticStructural/hingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/hingedBeam"))
writedlm("test/newTestDataGenerators/hingedBeam/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/hingedBeam/p2.txt", test_data_matrix(p2))
writedlm("test/newTestDataGenerators/hingedBeam/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/hingedBeam/M2.txt", test_data_matrix(M2))

# Static analysis of a hinged beam subjected to a distributed load and a rotational spring
include("../examples/staticStructural/hingedSpringedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/hingedSpringedBeam"))
writedlm("test/newTestDataGenerators/hingedSpringedBeam/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/hingedSpringedBeam/p2.txt", test_data_matrix(p2))
writedlm("test/newTestDataGenerators/hingedSpringedBeam/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/hingedSpringedBeam/M2.txt", test_data_matrix(M2))

# Static analysis of a T-frame hinged at the connection subjected to a distributed load
include("../examples/staticStructural/hingedTFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/hingedTFrame"))
writedlm("test/newTestDataGenerators/hingedTFrame/u1_beam1.txt", test_data_matrix(u1_beam1))
writedlm("test/newTestDataGenerators/hingedTFrame/u1_beam2.txt", test_data_matrix(u1_beam2))
writedlm("test/newTestDataGenerators/hingedTFrame/u3_beam1.txt", test_data_matrix(u3_beam1))
writedlm("test/newTestDataGenerators/hingedTFrame/u3_beam2.txt", test_data_matrix(u3_beam2))
writedlm("test/newTestDataGenerators/hingedTFrame/F3_beam1.txt", test_data_matrix(F3_beam1))
writedlm("test/newTestDataGenerators/hingedTFrame/M2_beam1.txt", test_data_matrix(M2_beam1))

# Static analysis of the Lee frame (a right-angled frame) with a dead load
include("../examples/staticStructural/LeeFrameDeadLoad.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/LeeFrameDeadLoad"))
writedlm("test/newTestDataGenerators/LeeFrameDeadLoad/u1_atForce.txt", test_data_matrix(u1_atForce))
writedlm("test/newTestDataGenerators/LeeFrameDeadLoad/u3_atForce.txt", test_data_matrix(u3_atForce))

# Static analysis of the Lee frame (a right-angled frame) with a follower load
include("../examples/staticStructural/LeeFrameFollowerLoad.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/LeeFrameFollowerLoad"))
writedlm("test/newTestDataGenerators/LeeFrameFollowerLoad/u1_atForce.txt", test_data_matrix(u1_atForce))
writedlm("test/newTestDataGenerators/LeeFrameFollowerLoad/u3_atForce.txt", test_data_matrix(u3_atForce))

# Static analysis of the pure bending test of the Pazy wing
include("../examples/staticStructural/PazyWingBendingTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingBendingTest"))
writedlm("test/newTestDataGenerators/PazyWingBendingTest/tip_OOP.txt", test_data_matrix(tip_OOP))

# Static analysis of the coupled torsion-bending test of the Pazy wing
include("../examples/staticStructural/PazyWingTorsionTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingTorsionTest"))
writedlm("test/newTestDataGenerators/PazyWingTorsionTest/tip_OOP.txt", test_data_matrix(tip_OOP))
writedlm("test/newTestDataGenerators/PazyWingTorsionTest/tip_twist.txt", test_data_matrix(tip_twist))

# Static analysis of a mid-loaded arch pinned at one end and clamped at the other
include("../examples/staticStructural/pinnedClampedArch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pinnedClampedArch"))
writedlm("test/newTestDataGenerators/pinnedClampedArch/u1.txt", test_data_matrix(u1_atForce))
writedlm("test/newTestDataGenerators/pinnedClampedArch/u3.txt", test_data_matrix(u3_atForce))

# Static analysis of a right-angled frame under a tip transverse follower force
include("../examples/staticStructural/rightAngledFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rightAngledFrame"))
writedlm("test/newTestDataGenerators/rightAngledFrame/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/rightAngledFrame/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/rightAngledFrame/tip_angle.txt", test_data_matrix(tip_angle))

# Static analysis of a right-angled frame under a 'buckling' load
include("../examples/staticStructural/rightAngledFrameBuckling.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rightAngledFrameBuckling"))
writedlm("test/newTestDataGenerators/rightAngledFrameBuckling/tip_u2.txt", test_data_matrix(tip_u2))

# Static analysis of an L-frame with a doubly-attached spring and tip load
include("../examples/staticStructural/springedLFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/springedLFrame"))
writedlm("test/newTestDataGenerators/springedLFrame/u3_b.txt", test_data_matrix(u3_b))

# Static analysis of a swept-back clamped beam with a driven flared folding wingtip
include("../examples/staticStructural/sweptBackFFWTWing.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptBackFFWTWing"))
writedlm("test/newTestDataGenerators/sweptBackFFWTWing/phi.txt", test_data_matrix(ϕHinge))
writedlm("test/newTestDataGenerators/sweptBackFFWTWing/hingeMoment.txt", test_data_matrix(hingeMoment))

# Static analysis of a semi-circular arch with tangential follower force
include("../examples/staticStructural/tangentiallyForcedArch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tangentiallyForcedArch"))
writedlm("test/newTestDataGenerators/tangentiallyForcedArch/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/tangentiallyForcedArch/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/tangentiallyForcedArch/tip_angle.txt", test_data_matrix(tip_angle))

# Static analysis of a cantilever with tip follower transverse force 
include("../examples/staticStructural/tipFollowerForceCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipFollowerForceCantilever"))
writedlm("test/newTestDataGenerators/tipFollowerForceCantilever/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/tipFollowerForceCantilever/tip_u3.txt", test_data_matrix(tip_u3))

# Static analysis of a cantilever with tip follower transverse force (force split over 2 BCs)
include("../examples/staticStructural/tipFollowerForceCantilever2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipFollowerForceCantilever2"))
writedlm("test/newTestDataGenerators/tipFollowerForceCantilever2/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/tipFollowerForceCantilever2/tip_u3.txt", test_data_matrix(tip_u3))

# Static analysis of a cantilever with tip moment
include("../examples/staticStructural/tipMomentCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipMomentCantilever"))
writedlm("test/newTestDataGenerators/tipMomentCantilever/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/tipMomentCantilever/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/tipMomentCantilever/tip_angle.txt", test_data_matrix(tip_angle))

# Static analysis of a semi-circular arch with transverse follower force
include("../examples/staticStructural/transverselyForcedArch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/transverselyForcedArch"))
writedlm("test/newTestDataGenerators/transverselyForcedArch/tip_u1.txt", test_data_matrix(tip_u1))
writedlm("test/newTestDataGenerators/transverselyForcedArch/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/transverselyForcedArch/tip_angle.txt", test_data_matrix(tip_angle))

# Static analysis of a beam with triangular distributed load
include("../examples/staticStructural/triangleLoadBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/triangleLoadBeam"))
writedlm("test/newTestDataGenerators/triangleLoadBeam/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/triangleLoadBeam/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/triangleLoadBeam/M2.txt", test_data_matrix(M2))

# Static analysis of a cantilever beam with a torsional spring attached between its middle and tip nodes, subjected to a torsional tip moment
include("../examples/staticStructural/twistDoublyAttachedSpringCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/twistDoublyAttachedSpringCantilever"))
writedlm("test/newTestDataGenerators/twistDoublyAttachedSpringCantilever/p1.txt", test_data_matrix(p1))
writedlm("test/newTestDataGenerators/twistDoublyAttachedSpringCantilever/M1.txt", test_data_matrix(M1))

# Static analysis of a beam with a central hinge driven by opposite end moments
include("../examples/staticStructural/drivenSpringedHingedBeam2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/drivenSpringedHingedBeam2"))
writedlm("test/newTestDataGenerators/drivenSpringedHingedBeam2/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/drivenSpringedHingedBeam2/p2.txt", test_data_matrix(p2))
writedlm("test/newTestDataGenerators/drivenSpringedHingedBeam2/M2.txt", test_data_matrix(M2))

# Static analysis of a swept cantilever under tip transverse loading
include("../examples/staticStructural/sweptBeamStaticLoading.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptBeamStaticLoading"))
writedlm("test/newTestDataGenerators/sweptBeamStaticLoading/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/sweptBeamStaticLoading/tip_twist.txt", test_data_matrix(tip_twist))
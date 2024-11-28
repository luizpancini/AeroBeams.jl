# New reference data for static structural tests

# Static analysis of an arch under a dead pressure load
include("../examples/archUnderDeadPressure.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/archUnderDeadPressure"))
writedlm("test/newTestDataGenerators/archUnderDeadPressure/mid_u3.txt", mid_u3)

# Static analysis of an arch under a follower pressure load
include("../examples/archUnderFollowerPressure.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/archUnderFollowerPressure"))
writedlm("test/newTestDataGenerators/archUnderFollowerPressure/mid_u3.txt", mid_u3)

# Static analysis of a biclamped, hinged beam under distributed and concentrated loads
include("../examples/biclampedHingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/biclampedHingedBeam"))
writedlm("test/newTestDataGenerators/biclampedHingedBeam/u3Mid.txt", u3Mid)
writedlm("test/newTestDataGenerators/biclampedHingedBeam/p2Left.txt", p2Left)
writedlm("test/newTestDataGenerators/biclampedHingedBeam/p2Right.txt", p2Right)

# Static analysis of a cantilever beam bending under self weight
include("../examples/cantileverUnderSelfWeight.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverUnderSelfWeight"))
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/u1.txt", u1)
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/u3.txt", u3)
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/F3.txt", F3)
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/M2.txt", M2)

# Static analysis of a cantilever beam with a tip spring in bending
include("../examples/cantileverWithTipSpring.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipSpring"))
writedlm("test/newTestDataGenerators/cantileverWithTipSpring/u3.txt", u3)
writedlm("test/newTestDataGenerators/cantileverWithTipSpring/F3.txt", F3)
writedlm("test/newTestDataGenerators/cantileverWithTipSpring/M2.txt", M2)

# Static analysis of a clamped beam, hinged at the middle with imposed hinge angle, under distributed and concentrated loads
include("../examples/clampedHingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedHingedBeam"))
writedlm("test/newTestDataGenerators/clampedHingedBeam/u1.txt", u1)
writedlm("test/newTestDataGenerators/clampedHingedBeam/u3.txt", u3)
writedlm("test/newTestDataGenerators/clampedHingedBeam/p2.txt", p2)
writedlm("test/newTestDataGenerators/clampedHingedBeam/F3.txt", F3)
writedlm("test/newTestDataGenerators/clampedHingedBeam/M2.txt", M2)

# Static analysis of a clamped beam rotated in 3D space, hinged at the middle with imposed hinge angle, under distributed and concentrated loads
include("../examples/clampedHingedBeamRotated.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedHingedBeamRotated"))
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/u1.txt", u1)
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/u2.txt", u2)
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/u3.txt", u3)
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/p2_b.txt", p2_b)
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/F3.txt", F3)
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/M2.txt", M2)

# Static analysis of a clamped beam, hinged at the middle, under distributed loads, with varying spring stiffnesses around the hinge
include("../examples/clampedHingedBeamSpringRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedHingedBeamSpringRange"))
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/u1.txt", u1)
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/u3.txt", u3)
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/p2.txt", p2)
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/F3.txt", F3)
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/M2.txt", M2)

# Static analysis of a clamped beam, hinged at the middle with imposed hinge angles about 3 directions, under distributed and concentrated loads
include("../examples/clampedUniversalHingeBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedUniversalHingeBeam"))
writedlm("test/newTestDataGenerators/clampedUniversalHingeBeam/u1.txt", u1)
writedlm("test/newTestDataGenerators/clampedUniversalHingeBeam/u2.txt", u2)
writedlm("test/newTestDataGenerators/clampedUniversalHingeBeam/u3.txt", u3)
writedlm("test/newTestDataGenerators/clampedUniversalHingeBeam/p1.txt", p1)
writedlm("test/newTestDataGenerators/clampedUniversalHingeBeam/p2.txt", p2)
writedlm("test/newTestDataGenerators/clampedUniversalHingeBeam/p3.txt", p3)
writedlm("test/newTestDataGenerators/clampedUniversalHingeBeam/F3.txt", F3)
writedlm("test/newTestDataGenerators/clampedUniversalHingeBeam/M2.txt", M2)

# Static analysis of composite laminates subjected to tip loads
include("../examples/compositeCantileverMD.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/compositeCantileverMD"))
for i in 1:3
    writedlm(string("test/newTestDataGenerators/compositeCantileverMD/u1_500mm_b",i,".txt"), hcat(u1_500mm[i,:]...)')
    writedlm(string("test/newTestDataGenerators/compositeCantileverMD/u2_500mm_b",i,".txt"), hcat(u2_500mm[i,:]...)')
    writedlm(string("test/newTestDataGenerators/compositeCantileverMD/u3_500mm_b",i,".txt"), hcat(u3_500mm[i,:]...)')
end

# Static analysis of a curved cantilever subjected to a tip dead force
include("../examples/curvedCantileverDeadLoad.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/curvedCantileverDeadLoad"))
writedlm("test/newTestDataGenerators/curvedCantileverDeadLoad/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/curvedCantileverDeadLoad/tip_u2.txt", tip_u2)
writedlm("test/newTestDataGenerators/curvedCantileverDeadLoad/tip_u3.txt", tip_u3)

# Static analysis of a curved cantilever subjected to a tip follower force
include("../examples/curvedCantileverStaticFollower.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/curvedCantileverStaticFollower"))
writedlm("test/newTestDataGenerators/curvedCantileverStaticFollower/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/curvedCantileverStaticFollower/tip_u2.txt", tip_u2)
writedlm("test/newTestDataGenerators/curvedCantileverStaticFollower/tip_u3.txt", tip_u3)

# Static analysis of a cantilever with distributed follower force
include("../examples/distributedLoadCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/distributedLoadCantilever"))
writedlm("test/newTestDataGenerators/distributedLoadCantilever/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/distributedLoadCantilever/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/distributedLoadCantilever/tip_angle.txt", tip_angle)

# Static analysis of a hinged beam subjected to a distributed load
include("../examples/hingedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/hingedBeam"))
writedlm("test/newTestDataGenerators/hingedBeam/u3.txt", u3)
writedlm("test/newTestDataGenerators/hingedBeam/p2.txt", p2)
writedlm("test/newTestDataGenerators/hingedBeam/F3.txt", F3)
writedlm("test/newTestDataGenerators/hingedBeam/M2.txt", M2)

# Static analysis of a hinged beam subjected to a distributed load and a rotational spring
include("../examples/hingedSpringedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/hingedSpringedBeam"))
writedlm("test/newTestDataGenerators/hingedSpringedBeam/u3.txt", u3)
writedlm("test/newTestDataGenerators/hingedSpringedBeam/p2.txt", p2)
writedlm("test/newTestDataGenerators/hingedSpringedBeam/F3.txt", F3)
writedlm("test/newTestDataGenerators/hingedSpringedBeam/M2.txt", M2)

# Static analysis of a T-frame hinged at the connection subjected to a distributed load
include("../examples/hingedTFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/hingedTFrame"))
writedlm("test/newTestDataGenerators/hingedTFrame/u1_beam1.txt", u1_beam1)
writedlm("test/newTestDataGenerators/hingedTFrame/u1_beam2.txt", u1_beam2)
writedlm("test/newTestDataGenerators/hingedTFrame/u3_beam1.txt", u3_beam1)
writedlm("test/newTestDataGenerators/hingedTFrame/u3_beam2.txt", u3_beam2)
writedlm("test/newTestDataGenerators/hingedTFrame/F3_beam1.txt", F3_beam1)
writedlm("test/newTestDataGenerators/hingedTFrame/M2_beam1.txt", M2_beam1)

# Static analysis of the Lee frame (a right-angled frame) with a dead load
include("../examples/LeeFrameDeadLoad.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/LeeFrameDeadLoad"))
writedlm("test/newTestDataGenerators/LeeFrameDeadLoad/u1_atForce.txt", u1_atForce)
writedlm("test/newTestDataGenerators/LeeFrameDeadLoad/u3_atForce.txt", u3_atForce)

# Static analysis of the Lee frame (a right-angled frame) with a follower load
include("../examples/LeeFrameFollowerLoad.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/LeeFrameFollowerLoad"))
writedlm("test/newTestDataGenerators/LeeFrameFollowerLoad/u1_atForce.txt", u1_atForce)
writedlm("test/newTestDataGenerators/LeeFrameFollowerLoad/u3_atForce.txt", u3_atForce)

# Static analysis of the pure bending test of the Pazy wing
include("../examples/PazyWingBendingTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingBendingTest"))
writedlm("test/newTestDataGenerators/PazyWingBendingTest/tip_OOP.txt", tip_OOP)

# Static analysis of the coupled torsion-bending test of the Pazy wing
include("../examples/PazyWingTorsionTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingTorsionTest"))
writedlm("test/newTestDataGenerators/PazyWingTorsionTest/tip_OOP.txt", tip_OOP)
writedlm("test/newTestDataGenerators/PazyWingTorsionTest/tip_twist.txt", tip_twist)

# Static analysis of a mid-loaded arch pinned at one end and clamped at the other
include("../examples/pinnedClampedArch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pinnedClampedArch"))
writedlm("test/newTestDataGenerators/pinnedClampedArch/u1.txt", u1_atForce)
writedlm("test/newTestDataGenerators/pinnedClampedArch/u3.txt", u3_atForce)

# Static analysis of a right-angled frame under a tip transverse follower force
include("../examples/rightAngledFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rightAngledFrame"))
writedlm("test/newTestDataGenerators/rightAngledFrame/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/rightAngledFrame/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/rightAngledFrame/tip_angle.txt", tip_angle)

# Static analysis of a right-angled frame under a 'buckling' load
include("../examples/rightAngledFrameBuckling.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rightAngledFrameBuckling"))
writedlm("test/newTestDataGenerators/rightAngledFrameBuckling/tip_u2.txt", tip_u2)

# Static analysis of an L-frame with a doubly-attached spring and tip load
include("../examples/springedLFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/springedLFrame"))
writedlm("test/newTestDataGenerators/springedLFrame/u3_b.txt", u3_b)

# Static analysis of a semi-circular arch with tangential follower force
include("../examples/tangentiallyForcedArch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tangentiallyForcedArch"))
writedlm("test/newTestDataGenerators/tangentiallyForcedArch/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/tangentiallyForcedArch/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/tangentiallyForcedArch/tip_angle.txt", tip_angle)

# Static analysis of a cantilever with tip follower transverse force 
include("../examples/tipFollowerForceCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipFollowerForceCantilever"))
writedlm("test/newTestDataGenerators/tipFollowerForceCantilever/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/tipFollowerForceCantilever/tip_u3.txt", tip_u3)

# Static analysis of a cantilever with tip follower transverse force (force split over 2 BCs)
include("../examples/tipFollowerForceCantilever2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipFollowerForceCantilever2"))
writedlm("test/newTestDataGenerators/tipFollowerForceCantilever2/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/tipFollowerForceCantilever2/tip_u3.txt", tip_u3)

# Static analysis of a cantilever with tip moment
include("../examples/tipMomentCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipMomentCantilever"))
writedlm("test/newTestDataGenerators/tipMomentCantilever/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/tipMomentCantilever/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/tipMomentCantilever/tip_angle.txt", tip_angle)

# Static analysis of a semi-circular arch with transverse follower force
include("../examples/transverselyForcedArch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/transverselyForcedArch"))
writedlm("test/newTestDataGenerators/transverselyForcedArch/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/transverselyForcedArch/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/transverselyForcedArch/tip_angle.txt", tip_angle)

# Static analysis of a beam with triangular distributed load
include("../examples/triangleLoadBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/triangleLoadBeam"))
writedlm("test/newTestDataGenerators/triangleLoadBeam/u3.txt", u3)
writedlm("test/newTestDataGenerators/triangleLoadBeam/F3.txt", F3)
writedlm("test/newTestDataGenerators/triangleLoadBeam/M2.txt", M2)
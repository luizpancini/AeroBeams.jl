# New reference data for dynamic structural tests

# Dynamic analysis of the axial vibration of a beam under a traction force applied suddenly
include("../examples/axialTractionCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/axialTractionCantilever"))
writedlm("test/newTestDataGenerators/axialTractionCantilever/u1_08.txt", u1_08)
writedlm("test/newTestDataGenerators/axialTractionCantilever/u1_10.txt", u1_10)

# Dynamic analysis of the free response of a beam clamped at both ends and subjected to an initial displacement profile
include("../examples/biclampedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/biclampedBeam"))
writedlm("test/newTestDataGenerators/biclampedBeam/u3_mid.txt", u3_mid)
writedlm("test/newTestDataGenerators/biclampedBeam/V3_mid.txt", V3_mid)
writedlm("test/newTestDataGenerators/biclampedBeam/Vdot3_mid.txt", Vdot3_mid)
writedlm("test/newTestDataGenerators/biclampedBeam/theta2_quarter.txt", θ2_quarter)
writedlm("test/newTestDataGenerators/biclampedBeam/Omega2_quarter.txt", Ω2_quarter)
writedlm("test/newTestDataGenerators/biclampedBeam/Omegadot2_quarter.txt", Ωdot2_quarter)

# Dynamic analysis of a composite cantilever beam under a tip sinusoidal load
include("../examples/compositeCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/compositeCantilever"))
writedlm("test/newTestDataGenerators/compositeCantilever/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/compositeCantilever/u2_tip.txt", u2_tip)
writedlm("test/newTestDataGenerators/compositeCantilever/u3_tip.txt", u3_tip)
writedlm("test/newTestDataGenerators/compositeCantilever/p1_tip.txt", p1_tip)
writedlm("test/newTestDataGenerators/compositeCantilever/p2_tip.txt", p2_tip)
writedlm("test/newTestDataGenerators/compositeCantilever/p3_tip.txt", p3_tip)
writedlm("test/newTestDataGenerators/compositeCantilever/F1_root.txt", F1_root)
writedlm("test/newTestDataGenerators/compositeCantilever/F2_root.txt", F2_root)
writedlm("test/newTestDataGenerators/compositeCantilever/F3_root.txt", F3_root)
writedlm("test/newTestDataGenerators/compositeCantilever/M1_root.txt", M1_root)
writedlm("test/newTestDataGenerators/compositeCantilever/M2_root.txt", M2_root)
writedlm("test/newTestDataGenerators/compositeCantilever/M3_root.txt", M3_root)

# Dynamic analysis of a curved cantilever subjected to a tip follower force
include("../examples/curvedCantileverDynamicFollower.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/curvedCantileverDynamicFollower"))
writedlm("test/newTestDataGenerators/curvedCantileverDynamicFollower/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/curvedCantileverDynamicFollower/u2_tip.txt", u2_tip)
writedlm("test/newTestDataGenerators/curvedCantileverDynamicFollower/u3_tip.txt", u3_tip)
writedlm("test/newTestDataGenerators/curvedCantileverDynamicFollower/F1_root.txt", F1_root)
writedlm("test/newTestDataGenerators/curvedCantileverDynamicFollower/F2_root.txt", F2_root)
writedlm("test/newTestDataGenerators/curvedCantileverDynamicFollower/F3_root.txt", F3_root)
writedlm("test/newTestDataGenerators/curvedCantileverDynamicFollower/M1_root.txt", M1_root)
writedlm("test/newTestDataGenerators/curvedCantileverDynamicFollower/M2_root.txt", M2_root)
writedlm("test/newTestDataGenerators/curvedCantileverDynamicFollower/M3_root.txt", M3_root)

# Dynamic analysis of a double pendulum released from rest
include("../examples/doublePendulum.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/doublePendulum"))
writedlm("test/newTestDataGenerators/doublePendulum/u1_hinge.txt", u1_hinge)
writedlm("test/newTestDataGenerators/doublePendulum/u3_hinge.txt", u3_hinge)
writedlm("test/newTestDataGenerators/doublePendulum/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/doublePendulum/u3_tip.txt", u3_tip)

# Dynamic analysis of a right-angled frame subjected to an out-of-plane force
include("../examples/elbowFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/elbowFrame"))
writedlm("test/newTestDataGenerators/elbowFrame/u3_elbow.txt", u3_elbow)
writedlm("test/newTestDataGenerators/elbowFrame/u3_tip.txt", u3_tip)

# Dynamic analysis of a very flexible beam subjected to loads yielding two-dimensional motion
include("../examples/flyingFlexibleBeam2D.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flyingFlexibleBeam2D"))
writedlm("test/newTestDataGenerators/flyingFlexibleBeam2D/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/flyingFlexibleBeam2D/u3_tip.txt", u3_tip)

# Dynamic analysis of a hinged beam in free flight
include("../examples/flyingScissors.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flyingScissors"))
writedlm("test/newTestDataGenerators/flyingScissors/u1_tipA.txt", u1_tipA)
writedlm("test/newTestDataGenerators/flyingScissors/u3_tipA.txt", u3_tipA)
writedlm("test/newTestDataGenerators/flyingScissors/u1_tipB.txt", u1_tipB)
writedlm("test/newTestDataGenerators/flyingScissors/u3_tipB.txt", u3_tipB)
writedlm("test/newTestDataGenerators/flyingScissors/u1_hinge.txt", u1_hinge)
writedlm("test/newTestDataGenerators/flyingScissors/u3_hinge.txt", u3_hinge)

# Dynamic analysis of a very flexible beam subjected to loads yielding two-dimensional motion
include("../examples/flyingSpaghetti2D.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flyingSpaghetti2D"))
writedlm("test/newTestDataGenerators/flyingSpaghetti2D/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/flyingSpaghetti2D/u3_tip.txt", u3_tip)

# Dynamic analysis of a very flexible beam subjected to loads yielding tri-dimensional motion
include("../examples/flyingSpaghetti3D.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flyingSpaghetti3D"))
writedlm("test/newTestDataGenerators/flyingSpaghetti3D/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/flyingSpaghetti3D/u2_tip.txt", u2_tip)
writedlm("test/newTestDataGenerators/flyingSpaghetti3D/u3_tip.txt", u3_tip)
writedlm("test/newTestDataGenerators/flyingSpaghetti3D/theta_tip.txt", θ_tip)

# Dynamic analysis of the free response of a beam subjected to initial displacement and velocity profiles
include("../examples/initialDispAndVelBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/initialDispAndVelBeam"))
writedlm("test/newTestDataGenerators/initialDispAndVelBeam/u3_quarter.txt", u3_quarter)
writedlm("test/newTestDataGenerators/initialDispAndVelBeam/V3_quarter.txt", V3_quarter)
writedlm("test/newTestDataGenerators/initialDispAndVelBeam/Vdot3_quarter.txt", Vdot3_quarter)
writedlm("test/newTestDataGenerators/initialDispAndVelBeam/theta2_root.txt", θ2_root)
writedlm("test/newTestDataGenerators/initialDispAndVelBeam/Omega2_mid.txt", Ω2_mid)
writedlm("test/newTestDataGenerators/initialDispAndVelBeam/Omegadot2_mid.txt", Ωdot2_mid)

# Dynamic analysis of the free response of a beam subjected to an initial displacement profile
include("../examples/initialDisplacementBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/initialDisplacementBeam"))
writedlm("test/newTestDataGenerators/initialDisplacementBeam/u3_quarter.txt", u3_quarter)
writedlm("test/newTestDataGenerators/initialDisplacementBeam/V3_quarter.txt", V3_quarter)
writedlm("test/newTestDataGenerators/initialDisplacementBeam/Vdot3_quarter.txt", Vdot3_quarter)
writedlm("test/newTestDataGenerators/initialDisplacementBeam/theta2_root.txt", θ2_root)
writedlm("test/newTestDataGenerators/initialDisplacementBeam/Omega2_mid.txt", Ω2_mid)
writedlm("test/newTestDataGenerators/initialDisplacementBeam/Omegadot2_mid.txt", Ωdot2_mid)

# Dynamic analysis of the free response of a beam subjected to an initial velocity profile
include("../examples/initialVelocityBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/initialVelocityBeam"))
writedlm("test/newTestDataGenerators/initialVelocityBeam/u3_quarter.txt", u3_quarter)
writedlm("test/newTestDataGenerators/initialVelocityBeam/V3_quarter.txt", V3_quarter)
writedlm("test/newTestDataGenerators/initialVelocityBeam/Vdot3_quarter.txt", Vdot3_quarter)
writedlm("test/newTestDataGenerators/initialVelocityBeam/theta2_root.txt", θ2_root)
writedlm("test/newTestDataGenerators/initialVelocityBeam/Omega2_mid.txt", Ω2_mid)
writedlm("test/newTestDataGenerators/initialVelocityBeam/Omegadot2_mid.txt", Ωdot2_mid)

# Dynamic analysis of joined beams under load
include("../examples/joinedBeams.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/joinedBeams"))
writedlm("test/newTestDataGenerators/joinedBeams/u1.txt", u1)
writedlm("test/newTestDataGenerators/joinedBeams/u2.txt", u2)
writedlm("test/newTestDataGenerators/joinedBeams/u3.txt", u3)

# Dynamic analysis of a pinned robot arm driven by a couple moment
include("../examples/momentDrivenRobotArm.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/momentDrivenRobotArm"))
writedlm("test/newTestDataGenerators/momentDrivenRobotArm/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/momentDrivenRobotArm/u3_tip.txt", u3_tip)

# Dynamic analysis of a pendulum released from rest
include("../examples/pendulum.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pendulum"))
writedlm("test/newTestDataGenerators/pendulum/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/pendulum/u3_tip.txt", u3_tip)

# Dynamic analysis of a clamped beam with root sinusoidal oscillation (ω = 9 Hz ≈ 2nd bending mode)
include("../examples/rootExcitationBeam1.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rootExcitationBeam1"))
writedlm("test/newTestDataGenerators/rootExcitationBeam1/u3b_root.txt", u3b_root)
writedlm("test/newTestDataGenerators/rootExcitationBeam1/u3b_tip.txt", u3b_tip)
writedlm("test/newTestDataGenerators/rootExcitationBeam1/V3_root.txt", V3_root)
writedlm("test/newTestDataGenerators/rootExcitationBeam1/V3_tip.txt", V3_tip)

# Dynamic analysis of a clamped beam with root sinusoidal oscillation (ω = 32 Hz ≈ 3rd bending mode)
include("../examples/rootExcitationBeam2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rootExcitationBeam2"))
writedlm("test/newTestDataGenerators/rootExcitationBeam2/u3b_root.txt", u3b_root)
writedlm("test/newTestDataGenerators/rootExcitationBeam2/u3b_tip.txt", u3b_tip)
writedlm("test/newTestDataGenerators/rootExcitationBeam2/V3_root.txt", V3_root)
writedlm("test/newTestDataGenerators/rootExcitationBeam2/V3_tip.txt", V3_tip)

# Dynamic analysis of a rotary shaft with specified rotation
include("../examples/rotaryShaft.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rotaryShaft"))
writedlm("test/newTestDataGenerators/rotaryShaft/pNum.txt", pNum)
writedlm("test/newTestDataGenerators/rotaryShaft/pdotNum.txt", pdotNum)
writedlm("test/newTestDataGenerators/rotaryShaft/OmegaNum.txt", ΩNum)
writedlm("test/newTestDataGenerators/rotaryShaft/OmegadotNum.txt", ΩdotNum)
writedlm("test/newTestDataGenerators/rotaryShaft/MNum.txt", MNum)

# Dynamic analysis of an articulated robot arm driven by specified rotation
include("../examples/rotationDrivenArticulatedRobotArm.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rotationDrivenArticulatedRobotArm"))
writedlm("test/newTestDataGenerators/rotationDrivenArticulatedRobotArm/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/rotationDrivenArticulatedRobotArm/u3_tip.txt", u3_tip)
writedlm("test/newTestDataGenerators/rotationDrivenArticulatedRobotArm/u1_hinge.txt", u1_hinge)
writedlm("test/newTestDataGenerators/rotationDrivenArticulatedRobotArm/u3_hinge.txt", u3_hinge)

# Dynamic analysis of a robot arm driven by specified rotation
include("../examples/rotationDrivenRobotArm.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rotationDrivenRobotArm"))
writedlm("test/newTestDataGenerators/rotationDrivenRobotArm/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/rotationDrivenRobotArm/u3_tip.txt", u3_tip)

# Dynamic analysis of the spin-up maneuver of a robot arm
include("../examples/spinupRobotArm.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/spinupRobotArm"))
writedlm("test/newTestDataGenerators/spinupRobotArm/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/spinupRobotArm/u2_tip.txt", u2_tip)
writedlm("test/newTestDataGenerators/spinupRobotArm/theta3_root.txt", θ3_root)

# Dynamic analysis of a pendulum with tip mass
include("../examples/tipPendulum.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipPendulum"))
writedlm("test/newTestDataGenerators/tipPendulum/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/tipPendulum/u3_tip.txt", u3_tip)

# Dynamic analysis of a cantilever with tip sinusoidal force
include("../examples/tipSineLoadedCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipSineLoadedCantilever"))
writedlm("test/newTestDataGenerators/tipSineLoadedCantilever/u3_tip.txt", u3_tip)
writedlm("test/newTestDataGenerators/tipSineLoadedCantilever/F3_root.txt", F3_root)
writedlm("test/newTestDataGenerators/tipSineLoadedCantilever/M2_root.txt", M2_root)
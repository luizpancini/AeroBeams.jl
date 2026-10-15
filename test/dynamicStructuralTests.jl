@testset "Dynamic analysis of the axial vibration of a beam under a traction force applied suddenly" begin
    include_example("examples/dynamicStructural/axialTractionCantilever.jl")
    @test matches_test_data(u1_08, "axialTractionCantilever", "u1_08.txt")
    @test matches_test_data(u1_10, "axialTractionCantilever", "u1_10.txt")
end

@testset "Dynamic analysis of the free response of a beam clamped at both ends and subjected to an initial displacement profile" begin
    include_example("examples/dynamicStructural/biclampedBeam.jl")
    @test matches_test_data(u3_mid, "biclampedBeam", "u3_mid.txt")
    @test matches_test_data(V3_mid, "biclampedBeam", "V3_mid.txt")
    @test matches_test_data(Vdot3_mid, "biclampedBeam", "Vdot3_mid.txt")
    @test matches_test_data(θ2_quarter, "biclampedBeam", "theta2_quarter.txt")
    @test matches_test_data(Ω2_quarter, "biclampedBeam", "Omega2_quarter.txt")
    @test matches_test_data(Ωdot2_quarter, "biclampedBeam", "Omegadot2_quarter.txt")
end

@testset "Dynamic analysis of a slanted double pendulum released from rest" begin
    include_example("examples/dynamicStructural/ChengSlantedDoublePendulum.jl")
    @test matches_test_data(u1_tip, "ChengSlantedDoublePendulum", "u1_tip.txt")
    @test matches_test_data(u2_tip, "ChengSlantedDoublePendulum", "u2_tip.txt")
    @test matches_test_data(u3_tip, "ChengSlantedDoublePendulum", "u3_tip.txt")
end

@testset "Dynamic analysis of a composite cantilever beam under a tip sinusoidal load" begin
    include_example("examples/dynamicStructural/compositeCantilever.jl")
    @test matches_test_data(u1_tip, "compositeCantilever", "u1_tip.txt")
    @test matches_test_data(u2_tip, "compositeCantilever", "u2_tip.txt")
    @test matches_test_data(u3_tip, "compositeCantilever", "u3_tip.txt")
    @test matches_test_data(p1_tip, "compositeCantilever", "p1_tip.txt")
    @test matches_test_data(p2_tip, "compositeCantilever", "p2_tip.txt")
    @test matches_test_data(p3_tip, "compositeCantilever", "p3_tip.txt")
    @test matches_test_data(F1_root, "compositeCantilever", "F1_root.txt")
    @test matches_test_data(F2_root, "compositeCantilever", "F2_root.txt")
    @test matches_test_data(F3_root, "compositeCantilever", "F3_root.txt")
    @test matches_test_data(M1_root, "compositeCantilever", "M1_root.txt")
    @test matches_test_data(M2_root, "compositeCantilever", "M2_root.txt")
    @test matches_test_data(M3_root, "compositeCantilever", "M3_root.txt")
end

@testset "Dynamic analysis of a curved cantilever subjected to a tip follower force" begin
    include_example("examples/dynamicStructural/curvedCantileverDynamicFollower.jl")
    @test matches_test_data(u1_tip, "curvedCantileverDynamicFollower", "u1_tip.txt")
    @test matches_test_data(u2_tip, "curvedCantileverDynamicFollower", "u2_tip.txt")
    @test matches_test_data(u3_tip, "curvedCantileverDynamicFollower", "u3_tip.txt")
    @test matches_test_data(F1_root, "curvedCantileverDynamicFollower", "F1_root.txt")
    @test matches_test_data(F2_root, "curvedCantileverDynamicFollower", "F2_root.txt")
    @test matches_test_data(F3_root, "curvedCantileverDynamicFollower", "F3_root.txt")
    @test matches_test_data(M1_root, "curvedCantileverDynamicFollower", "M1_root.txt")
    @test matches_test_data(M2_root, "curvedCantileverDynamicFollower", "M2_root.txt")
    @test matches_test_data(M3_root, "curvedCantileverDynamicFollower", "M3_root.txt")
end

@testset "Dynamic analysis of a double pendulum released from rest" begin
    include_example("examples/dynamicStructural/doublePendulum.jl")
    @test matches_test_data(u1_hinge, "doublePendulum", "u1_hinge.txt")
    @test matches_test_data(u3_hinge, "doublePendulum", "u3_hinge.txt")
    @test matches_test_data(u1_tip, "doublePendulum", "u1_tip.txt")
    @test matches_test_data(u3_tip, "doublePendulum", "u3_tip.txt")
end

@testset "Dynamic analysis of a right-angled frame subjected to an out-of-plane force" begin
    include_example("examples/dynamicStructural/elbowFrame.jl")
    @test matches_test_data(u3_elbow, "elbowFrame", "u3_elbow.txt")
    @test matches_test_data(u3_tip, "elbowFrame", "u3_tip.txt")
end

@testset "Dynamic analysis of a flexible beam subjected to loads yielding two-dimensional motion" begin
    include_example("examples/dynamicStructural/flyingFlexibleBeam2D.jl")
    @test matches_test_data(u1_tip, "flyingFlexibleBeam2D", "u1_tip.txt")
    @test matches_test_data(u3_tip, "flyingFlexibleBeam2D", "u3_tip.txt")
end

@testset "Dynamic analysis of a hinged beam in free flight" begin
    include_example("examples/dynamicStructural/flyingScissors.jl")
    @test matches_test_data(u1_tipA, "flyingScissors", "u1_tipA.txt")
    @test matches_test_data(u3_tipA, "flyingScissors", "u3_tipA.txt")
    @test matches_test_data(u1_tipB, "flyingScissors", "u1_tipB.txt")
    @test matches_test_data(u3_tipB, "flyingScissors", "u3_tipB.txt")
    @test matches_test_data(u1_hinge, "flyingScissors", "u1_hinge.txt")
    @test matches_test_data(u3_hinge, "flyingScissors", "u3_hinge.txt")
end

@testset "Dynamic analysis of a very flexible beam subjected to loads yielding two-dimensional motion" begin
    include_example("examples/dynamicStructural/flyingSpaghetti2D.jl")
    @test matches_test_data(u1_tip, "flyingSpaghetti2D", "u1_tip.txt")
    @test matches_test_data(u3_tip, "flyingSpaghetti2D", "u3_tip.txt")
end

@testset "Dynamic analysis of a very flexible beam subjected to loads yielding tri-dimensional motion" begin
    include_example("examples/dynamicStructural/flyingSpaghetti3D.jl")
    @test matches_test_data(u1_tip, "flyingSpaghetti3D", "u1_tip.txt")
    @test matches_test_data(u2_tip, "flyingSpaghetti3D", "u2_tip.txt")
    @test matches_test_data(u3_tip, "flyingSpaghetti3D", "u3_tip.txt")
end

@testset "Dynamic analysis of the free response of a beam subjected to initial displacement and velocity profiles" begin
    include_example("examples/dynamicStructural/initialDispAndVelBeam.jl")
    @test matches_test_data(u3_quarter, "initialDispAndVelBeam", "u3_quarter.txt")
    @test matches_test_data(V3_quarter, "initialDispAndVelBeam", "V3_quarter.txt")
    @test matches_test_data(Vdot3_quarter, "initialDispAndVelBeam", "Vdot3_quarter.txt")
    @test matches_test_data(θ2_root, "initialDispAndVelBeam", "theta2_root.txt")
    @test matches_test_data(Ω2_mid, "initialDispAndVelBeam", "Omega2_mid.txt")
    @test matches_test_data(Ωdot2_mid, "initialDispAndVelBeam", "Omegadot2_mid.txt")
end

@testset "Dynamic analysis of the free response of a beam subjected to an initial displacement profile" begin
    include_example("examples/dynamicStructural/initialDisplacementBeam.jl")
    @test matches_test_data(u3_quarter, "initialDisplacementBeam", "u3_quarter.txt")
    @test matches_test_data(V3_quarter, "initialDisplacementBeam", "V3_quarter.txt")
    @test matches_test_data(Vdot3_quarter, "initialDisplacementBeam", "Vdot3_quarter.txt")
    @test matches_test_data(θ2_root, "initialDisplacementBeam", "theta2_root.txt")
    @test matches_test_data(Ω2_mid, "initialDisplacementBeam", "Omega2_mid.txt")
    @test matches_test_data(Ωdot2_mid, "initialDisplacementBeam", "Omegadot2_mid.txt")
end

@testset "Dynamic analysis of the free response of a beam subjected to an initial velocity profile" begin
    include_example("examples/dynamicStructural/initialVelocityBeam.jl")
    @test matches_test_data(u3_quarter, "initialVelocityBeam", "u3_quarter.txt")
    @test matches_test_data(V3_quarter, "initialVelocityBeam", "V3_quarter.txt")
    @test matches_test_data(Vdot3_quarter, "initialVelocityBeam", "Vdot3_quarter.txt")
    @test matches_test_data(θ2_root, "initialVelocityBeam", "theta2_root.txt")
    @test matches_test_data(Ω2_mid, "initialVelocityBeam", "Omega2_mid.txt")
    @test matches_test_data(Ωdot2_mid, "initialVelocityBeam", "Omegadot2_mid.txt")
end

@testset "Dynamic analysis of joined beams under load" begin
    include_example("examples/dynamicStructural/joinedBeams.jl")
    @test matches_test_data(u1, "joinedBeams", "u1.txt")
    @test matches_test_data(u2, "joinedBeams", "u2.txt")
    @test matches_test_data(u3, "joinedBeams", "u3.txt")
end

@testset "Dynamic analysis of a pinned robot arm driven by a couple moment" begin
    include_example("examples/dynamicStructural/momentDrivenRobotArm.jl")
    @test matches_test_data(u1_tip, "momentDrivenRobotArm", "u1_tip.txt")
    @test matches_test_data(u3_tip, "momentDrivenRobotArm", "u3_tip.txt")
end

@testset "Dynamic analysis of a pendulum released from rest" begin
    include_example("examples/dynamicStructural/pendulum.jl")
    @test matches_test_data(u1_tip, "pendulum", "u1_tip.txt")
    @test matches_test_data(u3_tip, "pendulum", "u3_tip.txt")
end

@testset "Dynamic analysis of a clamped beam with root sinusoidal oscillation (ω = 9 Hz ≈ 2nd bending mode)" begin
    include_example("examples/dynamicStructural/rootExcitationBeam.jl")
    @test matches_test_data(u3b_root, "rootExcitationBeam", "u3b_root.txt")
    @test matches_test_data(u3b_tip, "rootExcitationBeam", "u3b_tip.txt")
    @test matches_test_data(V3_root, "rootExcitationBeam", "V3_root.txt")
    @test matches_test_data(V3_tip, "rootExcitationBeam", "V3_tip.txt")
end

@testset "Dynamic analysis of a rotary shaft with specified rotation" begin
    include_example("examples/dynamicStructural/rotaryShaft.jl")
    @test matches_test_data(pNum, "rotaryShaft", "pNum.txt")
    @test matches_test_data(pdotNum, "rotaryShaft", "pdotNum.txt")
    @test matches_test_data(ΩNum, "rotaryShaft", "OmegaNum.txt")
    @test matches_test_data(MNum, "rotaryShaft", "MNum.txt")
end

@testset "Dynamic analysis of an articulated robot arm driven by specified rotation" begin
    include_example("examples/dynamicStructural/rotationDrivenArticulatedRobotArm.jl")
    @test matches_test_data(u1_tip, "rotationDrivenArticulatedRobotArm", "u1_tip.txt")
    @test matches_test_data(u3_tip, "rotationDrivenArticulatedRobotArm", "u3_tip.txt")
    @test matches_test_data(u1_hinge, "rotationDrivenArticulatedRobotArm", "u1_hinge.txt")
    @test matches_test_data(u3_hinge, "rotationDrivenArticulatedRobotArm", "u3_hinge.txt")
end

@testset "Dynamic analysis of a robot arm driven by specified rotation" begin
    include_example("examples/dynamicStructural/rotationDrivenRobotArm.jl")
    @test matches_test_data(u1_tip, "rotationDrivenRobotArm", "u1_tip.txt")
    @test matches_test_data(u3_tip, "rotationDrivenRobotArm", "u3_tip.txt")
end

@testset "Dynamic analysis of the spin-up maneuver of a robot arm" begin
    include_example("examples/dynamicStructural/spinupRobotArm.jl")
    @test matches_test_data(u1_tip, "spinupRobotArm", "u1_tip.txt")
    @test matches_test_data(u2_tip, "spinupRobotArm", "u2_tip.txt")
    @test matches_test_data(θ3_root, "spinupRobotArm", "theta3_root.txt")
end

@testset "Dynamic analysis of a pendulum with tip mass" begin
    include_example("examples/dynamicStructural/tipPendulum.jl")
    @test matches_test_data(u1_tip, "tipPendulum", "u1_tip.txt")
    @test matches_test_data(u3_tip, "tipPendulum", "u3_tip.txt")
end

@testset "Dynamic analysis of a cantilever with tip sinusoidal force" begin
    include_example("examples/dynamicStructural/tipSineLoadedCantilever.jl")
    @test matches_test_data(u3_tip, "tipSineLoadedCantilever", "u3_tip.txt")
    @test matches_test_data(F3_root, "tipSineLoadedCantilever", "F3_root.txt")
    @test matches_test_data(M2_root, "tipSineLoadedCantilever", "M2_root.txt")
end

@testset "Dynamic analysis of a cantilever with linearly increasing tip moment" begin
    include_example("examples/dynamicStructural/tipMomentCantileverDynamic.jl")
    @test matches_test_data(κ2, "tipMomentCantileverDynamic", "kappa2.txt")
end
# Dynamic structural problems

@testset "Dynamic analysis of the axial vibration of a beam under a traction force applied suddenly" begin
    include("examples/axialTractionCantilever.jl")
    # Reference comparison
    @test u1_08[end]*1e3 ≈ u1_08_ref[end] rtol=1e-2
    @test u1_10[end]*1e3 ≈ u1_10_ref[end] rtol=2e-2
    # Self-comparison
    u1_08_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "axialTractionCantilever", "u1_08.txt"))
    u1_10_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "axialTractionCantilever", "u1_10.txt"))
    @test u1_08 ≈ u1_08_ atol=SELFatol rtol=SELFrtol
    @test u1_10 ≈ u1_10_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of the free response of a beam clamped at both ends and subjected to an initial displacement profile" begin
    include("examples/biclampedBeam.jl")
    # Analytical comparison
    @test u3_mid ≈ u3_mid_analytic rtol=1e-2
    @test V3_mid ≈ V3_mid_analytic rtol=5e-2
    @test θ2_quarter ≈ θ2_quarter_analytic rtol=2e-2
    @test Ω2_quarter ≈ Ω2_quarter_analytic rtol=10e-2
    # Self-comparison
    u3_mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedBeam", "u3_mid.txt"))
    V3_mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedBeam", "V3_mid.txt"))
    Vdot3_mid_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedBeam", "Vdot3_mid.txt"))
    θ2_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedBeam", "theta2_quarter.txt"))
    Ω2_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedBeam", "Omega2_quarter.txt"))
    Ωdot2_quarter_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "biclampedBeam", "Omegadot2_quarter.txt"))
    @test u3_mid ≈ u3_mid_ atol=SELFatol rtol=SELFrtol
    @test V3_mid ≈ V3_mid_ atol=SELFatol rtol=SELFrtol
    @test Vdot3_mid ≈ Vdot3_mid_ atol=SELFatol rtol=SELFrtol
    @test θ2_quarter ≈ θ2_quarter_ atol=SELFatol rtol=SELFrtol
    @test Ω2_quarter ≈ Ω2_quarter_ atol=SELFatol rtol=SELFrtol
    @test Ωdot2_quarter ≈ Ωdot2_quarter_ atol=SELFatol rtol=SELFrtol
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
    @test u1_tip ≈ u1_tip_ atol=SELFatol rtol=SELFrtol
    @test u2_tip ≈ u2_tip_ atol=SELFatol rtol=SELFrtol
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
    @test p1_tip ≈ p1_tip_ atol=SELFatol rtol=SELFrtol
    @test p2_tip ≈ p2_tip_ atol=SELFatol rtol=SELFrtol
    @test p3_tip ≈ p3_tip_ atol=SELFatol rtol=SELFrtol
    @test F1_root ≈ F1_root_ atol=1
    @test F2_root ≈ F2_root_ atol=1
    @test F3_root ≈ F3_root_ atol=1
    @test M1_root ≈ M1_root_ atol=1
    @test M2_root ≈ M2_root_ atol=1
    @test M3_root ≈ M3_root_ atol=1
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
    @test u1_tip ≈ u1_tip_ atol=SELFatol rtol=SELFrtol
    @test u2_tip ≈ u2_tip_ atol=SELFatol rtol=SELFrtol
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
    @test F1_root ≈ F1_root_ atol=SELFatol rtol=SELFrtol
    @test F2_root ≈ F2_root_ atol=SELFatol rtol=SELFrtol
    @test F3_root ≈ F3_root_ atol=SELFatol rtol=SELFrtol
    @test M1_root ≈ M1_root_ atol=SELFatol rtol=SELFrtol
    @test M2_root ≈ M2_root_ atol=SELFatol rtol=SELFrtol
    @test M3_root ≈ M3_root_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of a double pendulum released from rest" begin
    include("examples/doublePendulum.jl")
    # Self-comparison
    u1_hinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "doublePendulum", "u1_hinge.txt"))
    u3_hinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "doublePendulum", "u3_hinge.txt"))
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "doublePendulum", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "doublePendulum", "u3_tip.txt"))
    @test u1_hinge ≈ u1_hinge_ atol=SELFatol rtol=SELFrtol
    @test u3_hinge ≈ u3_hinge_ atol=SELFatol rtol=SELFrtol
    @test u1_tip ≈ u1_tip_ atol=SELFatol rtol=SELFrtol
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of a right-angled frame subjected to an out-of-plane force" begin
    include("examples/elbowFrame.jl")
    # Analytical comparison
    @test u3_elbow[end] ≈ u3ElbowRef[2,end] atol=0.5
    @test u3_tip[end] ≈ u3TipRef[2,end] atol=1
    # Self-comparison
    u3_elbow_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "elbowFrame", "u3_elbow.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "elbowFrame", "u3_tip.txt"))
    @test u3_elbow ≈ u3_elbow_ atol=SELFatol rtol=SELFrtol
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of a flexible beam subjected to loads yielding two-dimensional motion" begin
    include("examples/flyingFlexibleBeam2D.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingFlexibleBeam2D", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingFlexibleBeam2D", "u3_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol rtol=SELFrtol
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
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
    @test u1_tipA ≈ u1_tipA_ atol=SELFatol rtol=SELFrtol
    @test u3_tipA ≈ u3_tipA_ atol=SELFatol rtol=SELFrtol
    @test u1_tipB ≈ u1_tipB_ atol=SELFatol rtol=SELFrtol
    @test u3_tipB ≈ u3_tipB_ atol=SELFatol rtol=SELFrtol
    @test u1_hinge ≈ u1_hinge_ atol=SELFatol rtol=SELFrtol
    @test u3_hinge ≈ u3_hinge_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of a very flexible beam subjected to loads yielding two-dimensional motion" begin
    include("examples/flyingSpaghetti2D.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingSpaghetti2D", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingSpaghetti2D", "u3_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol rtol=SELFrtol
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of a very flexible beam subjected to loads yielding tri-dimensional motion" begin
    include("examples/flyingSpaghetti3D.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingSpaghetti3D", "u1_tip.txt"))
    u2_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingSpaghetti3D", "u2_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingSpaghetti3D", "u3_tip.txt"))
    θ_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flyingSpaghetti3D", "theta_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol rtol=SELFrtol
    @test u2_tip ≈ u2_tip_ atol=SELFatol rtol=SELFrtol
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
    @test θ_tip ≈ θ_tip_ atol=SELFatol rtol=SELFrtol
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
    @test u3_quarter ≈ u3_quarter_ atol=SELFatol rtol=SELFrtol
    @test V3_quarter ≈ V3_quarter_ atol=SELFatol rtol=SELFrtol
    @test Vdot3_quarter ≈ Vdot3_quarter_ atol=SELFatol rtol=SELFrtol
    @test θ2_root ≈ θ2_root_ atol=SELFatol rtol=SELFrtol
    @test Ω2_mid ≈ Ω2_mid_ atol=SELFatol rtol=SELFrtol
    @test Ωdot2_mid ≈ Ωdot2_mid_ atol=SELFatol rtol=SELFrtol
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
    @test u3_quarter ≈ u3_quarter_ atol=SELFatol rtol=SELFrtol
    @test V3_quarter ≈ V3_quarter_ atol=SELFatol rtol=SELFrtol
    @test Vdot3_quarter ≈ Vdot3_quarter_ atol=SELFatol rtol=SELFrtol
    @test θ2_root ≈ θ2_root_ atol=SELFatol rtol=SELFrtol
    @test Ω2_mid ≈ Ω2_mid_ atol=SELFatol rtol=SELFrtol
    @test Ωdot2_mid ≈ Ωdot2_mid_ atol=SELFatol rtol=SELFrtol
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
    @test u3_quarter ≈ u3_quarter_ atol=SELFatol rtol=SELFrtol
    @test V3_quarter ≈ V3_quarter_ atol=SELFatol rtol=SELFrtol
    @test Vdot3_quarter ≈ Vdot3_quarter_ atol=SELFatol rtol=SELFrtol
    @test θ2_root ≈ θ2_root_ atol=SELFatol rtol=SELFrtol
    @test Ω2_mid ≈ Ω2_mid_ atol=SELFatol rtol=SELFrtol
    @test Ωdot2_mid ≈ Ωdot2_mid_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of joined beams under load" begin
    include("examples/joinedBeams.jl")
    # Self-comparison
    u1_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "joinedBeams", "u1.txt"))
    u2_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "joinedBeams", "u2.txt"))
    u3_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "joinedBeams", "u3.txt"))
    @test u1 ≈ u1_ atol=SELFatol rtol=SELFrtol
    @test u2 ≈ u2_ atol=SELFatol rtol=SELFrtol
    @test u3 ≈ u3_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of a pinned robot arm driven by a couple moment" begin
    include("examples/momentDrivenRobotArm.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "momentDrivenRobotArm", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "momentDrivenRobotArm", "u3_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol rtol=SELFrtol
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of a pendulum released from rest" begin
    include("examples/pendulum.jl")
    # Reference (analytical) comparison
    @test u1_tip ≈ u1_tip_analytical rtol=5e-2
    @test u3_tip ≈ u3_tip_analytical rtol=5e-2
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "pendulum", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "pendulum", "u3_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol rtol=SELFrtol
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of a clamped beam with root sinusoidal oscillation (ω = 9 Hz ≈ 2nd bending mode)" begin
    include("examples/rootExcitationBeam1.jl")
    # Self-comparison
    u3b_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam1", "u3b_root.txt"))
    u3b_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam1", "u3b_tip.txt"))
    V3_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam1", "V3_root.txt"))
    V3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam1", "V3_tip.txt"))
    @test u3b_root ≈ u3b_root_ atol=SELFatol rtol=SELFrtol
    @test u3b_tip ≈ u3b_tip_ atol=SELFatol rtol=SELFrtol
    @test V3_root ≈ V3_root_ atol=SELFatol rtol=SELFrtol
    @test V3_tip ≈ V3_tip_ atol=SELFatol rtol=SELFrtol
end

# Reduce CI time
# @testset "Dynamic analysis of a clamped beam with root sinusoidal oscillation (ω = 32 Hz ≈ 3rd bending mode)" begin
#     include("examples/rootExcitationBeam2.jl")
#     # Self-comparison
#     u3b_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam2", "u3b_root.txt"))
#     u3b_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam2", "u3b_tip.txt"))
#     V3_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam2", "V3_root.txt"))
#     V3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rootExcitationBeam2", "V3_tip.txt"))
#     @test u3b_root ≈ u3b_root_ atol=SELFatol rtol=SELFrtol
#     @test u3b_tip ≈ u3b_tip_ atol=SELFatol rtol=SELFrtol
#     @test V3_root ≈ V3_root_ atol=SELFatol rtol=SELFrtol
#     @test V3_tip ≈ V3_tip_ atol=SELFatol rtol=SELFrtol
# end

@testset "Dynamic analysis of a rotary shaft with specified rotation" begin
    include("examples/rotaryShaft.jl")
    # Analytical comparison
    @test pNum ≈ p.(t) rtol = 1e-2
    @test pdotNum ≈ pdot.(t) rtol = 1e-2
    @test ΩNum ≈ θdot.(t) rtol = 1e-2
    @test ΩdotNum ≈ θddot.(t) rtol = 1e-2
    @test MNum ≈ -θddot.(t)*ρ*Ix rtol = 0.1
    # Self-comparison
    pNum_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotaryShaft", "pNum.txt"))
    pdotNum_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotaryShaft", "pdotNum.txt"))
    ΩNum_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotaryShaft", "OmegaNum.txt"))
    ΩdotNum_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotaryShaft", "OmegadotNum.txt"))
    MNum_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotaryShaft", "MNum.txt"))
    @test pNum ≈ pNum_ atol=SELFatol rtol=SELFrtol
    @test pdotNum ≈ pdotNum_ atol=SELFatol rtol=SELFrtol
    @test ΩNum ≈ ΩNum_ atol=SELFatol rtol=SELFrtol
    # @test ΩdotNum ≈ ΩdotNum_ atol=SELFatol rtol=SELFrtol
    @test MNum ≈ MNum_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of an articulated robot arm driven by specified rotation" begin
    include("examples/rotationDrivenArticulatedRobotArm.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotationDrivenArticulatedRobotArm", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotationDrivenArticulatedRobotArm", "u3_tip.txt"))
    u1_hinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotationDrivenArticulatedRobotArm", "u1_hinge.txt"))
    u3_hinge_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotationDrivenArticulatedRobotArm", "u3_hinge.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol rtol=SELFrtol
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
    @test u1_hinge ≈ u1_hinge_ atol=SELFatol rtol=SELFrtol
    @test u3_hinge ≈ u3_hinge_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of a robot arm driven by specified rotation" begin
    include("examples/rotationDrivenRobotArm.jl")
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotationDrivenRobotArm", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "rotationDrivenRobotArm", "u3_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol rtol=SELFrtol
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of the spin-up maneuver of a robot arm" begin
    include("examples/spinupRobotArm.jl")
    # Reference comparison
    @test θ3_root[end] ≈ θ(t[end]) rtol=1e-2
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "spinupRobotArm", "u1_tip.txt"))
    u2_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "spinupRobotArm", "u2_tip.txt"))
    θ3_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "spinupRobotArm", "theta3_root.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol rtol=SELFrtol
    @test u2_tip ≈ u2_tip_ atol=SELFatol rtol=SELFrtol
    @test θ3_root ≈ θ3_root_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of a pendulum with tip mass" begin
    include("examples/tipPendulum.jl")
    # Analytical comparison
    @test u1_tip[end]/L ≈ u1_tip_analytical[end]/L atol=1e-3
    @test u3_tip[end]/L ≈ u3_tip_analytical[end]/L atol=1e-3
    # Self-comparison
    u1_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipPendulum", "u1_tip.txt"))
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipPendulum", "u3_tip.txt"))
    @test u1_tip ≈ u1_tip_ atol=SELFatol rtol=SELFrtol
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
end

@testset "Dynamic analysis of a cantilever with tip sinusoidal force" begin
    include("examples/tipSineLoadedCantilever.jl")
    # Self-comparison
    u3_tip_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipSineLoadedCantilever", "u3_tip.txt"))
    F3_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipSineLoadedCantilever", "F3_root.txt"))
    M2_root_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "tipSineLoadedCantilever", "M2_root.txt"))
    @test u3_tip ≈ u3_tip_ atol=SELFatol rtol=SELFrtol
    @test F3_root ≈ F3_root_ atol=SELFatol rtol=SELFrtol
    @test M2_root ≈ M2_root_ atol=SELFatol rtol=SELFrtol
end
GC.gc()
sleep(1)
using DelimitedFiles

# Static analysis of an arch under a dead pressure load
include("../examples/archUnderDeadPressure.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/archUnderDeadPressure"))
writedlm("test/newTestDataGenerators/archUnderDeadPressure/mid_u3.txt", mid_u3)

# Static analysis of an arch under a follower pressure load
include("../examples/archUnderFollowerPressure.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/archUnderFollowerPressure"))
writedlm("test/newTestDataGenerators/archUnderFollowerPressure/mid_u3.txt", mid_u3)

# Dynamic analysis of the axial vibration of a beam under a traction force applied suddenly
include("../examples/axialTractionCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/axialTractionCantilever"))
writedlm("test/newTestDataGenerators/axialTractionCantilever/u1_08.txt", u1_08)
writedlm("test/newTestDataGenerators/axialTractionCantilever/u1_10.txt", u1_10)

# Modal analysis of the axial vibration of a beam under clamped-clamped boundary conditions
include("../examples/beamAxialVibrationCC.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamAxialVibrationCC"))
writedlm("test/newTestDataGenerators/beamAxialVibrationCC/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamAxialVibrationCC/u1_modeShapes.txt", u1_modeShapes)

# Modal analysis of the axial vibration of a beam under clamped-free boundary conditions
include("../examples/beamAxialVibrationCF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamAxialVibrationCF"))
writedlm("test/newTestDataGenerators/beamAxialVibrationCF/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamAxialVibrationCF/u1_modeShapes.txt", u1_modeShapes)

# Modal analysis of the axial vibration of a beam under free-free boundary conditions
include("../examples/beamAxialVibrationFF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamAxialVibrationFF"))
writedlm("test/newTestDataGenerators/beamAxialVibrationFF/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamAxialVibrationFF/u1_modeShapes.txt", u1_modeShapes)

# Modal analysis of the bending vibration of a beam under clamped-clamped boundary conditions
include("../examples/beamBendingVibrationCC.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationCC"))
writedlm("test/newTestDataGenerators/beamBendingVibrationCC/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamBendingVibrationCC/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the bending vibration of a beam under clamped-free boundary conditions
include("../examples/beamBendingVibrationCF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationCF"))
writedlm("test/newTestDataGenerators/beamBendingVibrationCF/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamBendingVibrationCF/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the bending vibration of a beam under clamped-pinned boundary conditions
include("../examples/beamBendingVibrationCP.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationCP"))
writedlm("test/newTestDataGenerators/beamBendingVibrationCP/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamBendingVibrationCP/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the bending vibration of a beam under clamped-sliding boundary conditions
include("../examples/beamBendingVibrationCS.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationCS"))
writedlm("test/newTestDataGenerators/beamBendingVibrationCS/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamBendingVibrationCS/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the bending vibration of a beam under free-free boundary conditions
include("../examples/beamBendingVibrationFF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationFF"))
writedlm("test/newTestDataGenerators/beamBendingVibrationFF/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamBendingVibrationFF/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the bending vibration of a beam under pinned-pinned boundary conditions
include("../examples/beamBendingVibrationPP.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationPP"))
writedlm("test/newTestDataGenerators/beamBendingVibrationPP/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamBendingVibrationPP/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the torsional vibration of a beam under clamped-clamped boundary conditions
include("../examples/beamTorsionalVibrationCC.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamTorsionalVibrationCC"))
writedlm("test/newTestDataGenerators/beamTorsionalVibrationCC/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamTorsionalVibrationCC/p1_modeShapes.txt", p1_modeShapes)

# Modal analysis of the torsional vibration of a beam under clamped-free boundary conditions
include("../examples/beamTorsionalVibrationCF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamTorsionalVibrationCF"))
writedlm("test/newTestDataGenerators/beamTorsionalVibrationCF/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamTorsionalVibrationCF/p1_modeShapes.txt", p1_modeShapes)

# Modal analysis of the torsional vibration of a beam under free-free boundary conditions
include("../examples/beamTorsionalVibrationFF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamTorsionalVibrationFF"))
writedlm("test/newTestDataGenerators/beamTorsionalVibrationFF/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamTorsionalVibrationFF/p1_modeShapes.txt", p1_modeShapes)

# Dynamic analysis of the free response of a beam clamped at both ends and subjected to an initial displacement profile
include("../examples/biclampedBeam.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/biclampedBeam"))
writedlm("test/newTestDataGenerators/biclampedBeam/u3_mid.txt", u3_mid)
writedlm("test/newTestDataGenerators/biclampedBeam/V3_mid.txt", V3_mid)
writedlm("test/newTestDataGenerators/biclampedBeam/Vdot3_mid.txt", Vdot3_mid)
writedlm("test/newTestDataGenerators/biclampedBeam/theta2_quarter.txt", θ2_quarter)
writedlm("test/newTestDataGenerators/biclampedBeam/Omega2_quarter.txt", Ω2_quarter)
writedlm("test/newTestDataGenerators/biclampedBeam/Omegadot2_quarter.txt", Ωdot2_quarter)

# Flutter analysis of the Blended-Wing-Body flying wing
include("../examples/BWBflutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/BWBflutter"))
writedlm("test/newTestDataGenerators/BWBflutter/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/BWBflutter/damps.txt", damps)

# Trim analysis of the Blended-Wing-Body flying wing in free flight
include("../examples/BWBtrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/BWBtrim"))
writedlm("test/newTestDataGenerators/BWBtrim/trimAoA.txt", trimAoA)
writedlm("test/newTestDataGenerators/BWBtrim/trimThrust.txt", trimThrust)
writedlm("test/newTestDataGenerators/BWBtrim/trimDelta.txt", trimδ)

# Static analysis of a cantilever beam bending under self weight
include("../examples/cantileverUnderSelfWeight.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverUnderSelfWeight"))
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/u1.txt", u1)
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/u3.txt", u3)
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/F3.txt", F3)
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/M2.txt", M2)

# Modal analysis of a cantilever beam with a tip axial inertia
include("../examples/cantileverWithTipAxialMassEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipAxialMassEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipAxialMassEigen/freqsNorm.txt", freqsNorm)
writedlm("test/newTestDataGenerators/cantileverWithTipAxialMassEigen/u1_modeShapes.txt", u1_modeShapes)

# Modal analysis of a cantilever beam with a tip axial spring
include("../examples/cantileverWithTipAxialSpringEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipAxialSpringEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipAxialSpringEigen/freqsNorm.txt", freqsNorm)
writedlm("test/newTestDataGenerators/cantileverWithTipAxialSpringEigen/u1_modeShapes.txt", u1_modeShapes)

# Static analysis of a cantilever beam with a tip spring in bending
include("../examples/cantileverWithTipSpring.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipSpring"))
writedlm("test/newTestDataGenerators/cantileverWithTipSpring/u3.txt", u3)
writedlm("test/newTestDataGenerators/cantileverWithTipSpring/F3.txt", F3)
writedlm("test/newTestDataGenerators/cantileverWithTipSpring/M2.txt", M2)

# Modal analysis of a cantilever beam with a tip spring in bending
include("../examples/cantileverWithTipSpringEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipSpringEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipSpringEigen/freqsNorm.txt", freqsNorm)
writedlm("test/newTestDataGenerators/cantileverWithTipSpringEigen/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of a cantilever beam with a tip torsional inertia
include("../examples/cantileverWithTipTorsionalInertiaEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipTorsionalInertiaEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipTorsionalInertiaEigen/freqsNorm.txt", freqsNorm)
writedlm("test/newTestDataGenerators/cantileverWithTipTorsionalInertiaEigen/p1_modeShapes.txt", p1_modeShapes)

# Modal analysis of a cantilever beam with a tip torsional spring
include("../examples/cantileverWithTipTorsionalSpringEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipTorsionalSpringEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipTorsionalSpringEigen/freqsNorm.txt", freqsNorm)
writedlm("test/newTestDataGenerators/cantileverWithTipTorsionalSpringEigen/p1_modeShapes.txt", p1_modeShapes)

# Modal analysis a beam clamped at one end, simply-supported at the other and with a tip inertia
include("../examples/clampedSSBeamWIthTipInertia.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedSSBeamWIthTipInertia"))
writedlm("test/newTestDataGenerators/clampedSSBeamWIthTipInertia/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/clampedSSBeamWIthTipInertia/u3_modeShapes.txt", u3_modeShapes)

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

# Static analysis of composite laminates subjected to tip loads
include("../examples/compositeCantileverMD.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/compositeCantileverMD"))
for i in 1:3
    writedlm(string("test/newTestDataGenerators/compositeCantileverMD/u1_500mm_b",i,".txt"), hcat(u1_500mm[i,:]...)')
    writedlm(string("test/newTestDataGenerators/compositeCantileverMD/u2_500mm_b",i,".txt"), hcat(u2_500mm[i,:]...)')
    writedlm(string("test/newTestDataGenerators/compositeCantileverMD/u3_500mm_b",i,".txt"), hcat(u3_500mm[i,:]...)')
end

# Dynamic analysis the conventional HALE aircraft undergoing a checked pitch maneuver
include("../examples/conventionalHALECheckedPitchManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALECheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/conventionalHALECheckedPitchManeuver/rootAoA.txt", rootAoA)
writedlm("test/newTestDataGenerators/conventionalHALECheckedPitchManeuver/Deltau3.txt", Δu3)

# Trim analysis the conventional HALE aircraft in free flight (considering aerodynamics from stabilizers and thrust)
include("../examples/conventionalHALEfullTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALEfullTrim"))
writedlm("test/newTestDataGenerators/conventionalHALEfullTrim/trimAoA.txt", trimAoA)
writedlm("test/newTestDataGenerators/conventionalHALEfullTrim/trimThrust.txt", trimThrust)
writedlm("test/newTestDataGenerators/conventionalHALEfullTrim/trimDelta.txt", trimδ)

# Flutter analysis the conventional HALE aircraft in free flight with structural stiffness as the varying parameter
include("../examples/conventionalHALELambdaRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALELambdaRange"))
writedlm("test/newTestDataGenerators/conventionalHALELambdaRange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/conventionalHALELambdaRange/damps.txt", damps)

# Flutter analysis the conventional HALE aircraft in free flight with airspeed and structural stiffness as the varying parameters
include("../examples/conventionalHALELURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALELURange"))
for (i,λ) in enumerate(λRange)
    writedlm(string("test/newTestDataGenerators/conventionalHALELURange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/conventionalHALELURange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Trim analysis the conventional HALE aircraft in free flight at rigid and flexible configurations (neglecting aerodynamics from stabilizers and thrust)
include("../examples/conventionalHALEtrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALEtrim"))
writedlm("test/newTestDataGenerators/conventionalHALEtrim/trimAoA.txt", trimAoA)

# Flutter analysis the conventional HALE aircraft in free flight with airspeed as the varying parameter
include("../examples/conventionalHALEURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALEURange"))
writedlm("test/newTestDataGenerators/conventionalHALEURange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/conventionalHALEURange/damps.txt", damps)

# Static analysis of a curved cantilever subjected to a tip dead force
include("../examples/curvedCantileverDeadLoad.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/curvedCantileverDeadLoad"))
writedlm("test/newTestDataGenerators/curvedCantileverDeadLoad/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/curvedCantileverDeadLoad/tip_u2.txt", tip_u2)
writedlm("test/newTestDataGenerators/curvedCantileverDeadLoad/tip_u3.txt", tip_u3)

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

# Dynamic analysis of a double pendulum released from rest
include("../examples/doublePendulum.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/doublePendulum"))
writedlm("test/newTestDataGenerators/doublePendulum/u1_hinge.txt", u1_hinge)
writedlm("test/newTestDataGenerators/doublePendulum/u3_hinge.txt", u3_hinge)
writedlm("test/newTestDataGenerators/doublePendulum/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/doublePendulum/u3_tip.txt", u3_tip)

# Dynamic analysis of a harmonically pitching airfoil, using the dynamic stall model
include("../examples/DSModelTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/DSModelTest"))
writedlm("test/newTestDataGenerators/DSModelTest/cn.txt", cn)
writedlm("test/newTestDataGenerators/DSModelTest/cm.txt", cm)
writedlm("test/newTestDataGenerators/DSModelTest/ct.txt", ct)

# Dynamic analysis of a right-angled frame subjected to an out-of-plane force
include("../examples/elbowFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/elbowFrame"))
writedlm("test/newTestDataGenerators/elbowFrame/u3_elbow.txt", u3_elbow)
writedlm("test/newTestDataGenerators/elbowFrame/u3_tip.txt", u3_tip)

# Dynamic analysis of an airfoils with harmonic flap deflection profile
include("../examples/flapOscillation.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flapOscillation"))
writedlm("test/newTestDataGenerators/flapOscillation/cn.txt", cn)
writedlm("test/newTestDataGenerators/flapOscillation/cm.txt", cm)

# Dynamic analysis of two airfoils with linked harmonic flap deflection profiles
include("../examples/flapOscillationLinked.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flapOscillationLinked"))
writedlm("test/newTestDataGenerators/flapOscillationLinked/cnMaster.txt", cnMaster)
writedlm("test/newTestDataGenerators/flapOscillationLinked/cmMaster.txt", cmMaster)
writedlm("test/newTestDataGenerators/flapOscillationLinked/cnSlave.txt", cnSlave)
writedlm("test/newTestDataGenerators/flapOscillationLinked/cmSlave.txt", cmSlave)

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

# Trim analysis (reaction loads check) of a beam loaded at the middle
include("../examples/freeBeamTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/freeBeamTrim"))
writedlm("test/newTestDataGenerators/freeBeamTrim/F3.txt", F3)
writedlm("test/newTestDataGenerators/freeBeamTrim/M2.txt", M2)

# Dynamic analysis the Helios flying-wing undergoing a checked pitch maneuver
include("../examples/heliosCheckedPitchManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosCheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/heliosCheckedPitchManeuver/rootAoA.txt", rootAoA)
writedlm("test/newTestDataGenerators/heliosCheckedPitchManeuver/Deltau3.txt", Δu3)

# Flutter analysis the Helios flying-wing in free flight with payload and structural stiffness as the varying parameters
include("../examples/heliosFlutterPLambdaRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosFlutterPLambdaRange"))
for (i,λ) in enumerate(λRange)
    writedlm(string("test/newTestDataGenerators/heliosFlutterPLambdaRange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/heliosFlutterPLambdaRange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Flutter analysis the Helios flying-wing in free flight with payload as the varying parameter
include("../examples/heliosFlutterPRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosFlutterPRange"))
writedlm("test/newTestDataGenerators/heliosFlutterPRange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/heliosFlutterPRange/damps.txt", damps)

# Flutter analysis the Helios flying-wing in free flight with airspeed as the varying parameter
include("../examples/heliosFlutterURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosFlutterURange"))
writedlm("test/newTestDataGenerators/heliosFlutterURange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/heliosFlutterURange/damps.txt", damps)

# Trim analysis of the Helios flying-wing
include("../examples/heliosTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosTrim"))
writedlm("test/newTestDataGenerators/heliosTrim/trimAoA.txt", trimAoA)
writedlm("test/newTestDataGenerators/heliosTrim/trimThrust.txt", trimThrust)
writedlm("test/newTestDataGenerators/heliosTrim/trimDelta.txt", trimδ)

# Flutter analysis of the wing of the Helios flying-wing
include("../examples/heliosWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosWingFlutter"))
writedlm("test/newTestDataGenerators/heliosWingFlutter/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/heliosWingFlutter/damps.txt", damps)

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

# Trim analysis (reaction loads check) of a simply-supported beam loaded at the middle
include("../examples/midLoadedBeamTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/midLoadedBeamTrim"))
writedlm("test/newTestDataGenerators/midLoadedBeamTrim/F3.txt", F3)
writedlm("test/newTestDataGenerators/midLoadedBeamTrim/M2.txt", M2)

# Dynamic analysis of a pinned robot arm driven by a couple moment
include("../examples/momentDrivenRobotArm.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/momentDrivenRobotArm"))
writedlm("test/newTestDataGenerators/momentDrivenRobotArm/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/momentDrivenRobotArm/u3_tip.txt", u3_tip)

# One-minus-cosine gust response of an airfoil section at several pitch angles
include("../examples/OMCgustTests.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/OMCgustTests"))
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
            relPath = string("/test/newTestDataGenerators/OMCgustTests/OMCgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase)
            absPath = string(pwd(),relPath)
            mkpath(absPath)
            writedlm(string("test/newTestDataGenerators/OMCgustTests/OMCgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/tau",i,j,k,".txt"), τ[i,j,k])
            writedlm(string("test/newTestDataGenerators/OMCgustTests/OMCgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dcl",i,j,k,".txt"), Δcl[i,j,k])
            writedlm(string("test/newTestDataGenerators/OMCgustTests/OMCgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dclRef",i,j,k,".txt"), ΔclRef[i,j,k])
        end
    end
end

# Eigen-analysis of the Pazy wing with flared folding wing tip (FFWT)
include("../examples/PazyFFWTeigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTeigen"))
writedlm("test/newTestDataGenerators/PazyFFWTeigen/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/PazyFFWTeigen/damps.txt", damps)

# Steady analysis of the Pazy wing with flared folding wing tip (FFWT)
include("../examples/PazyFFWTsteady.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteady"))
writedlm("test/newTestDataGenerators/PazyFFWTsteady/u3_of_x1.txt", u3_of_x1)
writedlm("test/newTestDataGenerators/PazyFFWTsteady/p2_of_x1.txt", p2_of_x1)
writedlm("test/newTestDataGenerators/PazyFFWTsteady/M2_of_x1.txt", M2_of_x1)

# Static analysis of the pure bending test of the Pazy wing
include("../examples/PazyWingBendingTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingBendingTest"))
writedlm("test/newTestDataGenerators/PazyWingBendingTest/tip_OOP.txt", tip_OOP)

# Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over time
include("../examples/PazyWingContinuous1DGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingContinuous1DGust"))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tqSpan_ct.txt", tqSpan_ct)

# Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over space
include("../examples/PazyWingContinuous1DSpaceGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingContinuous1DSpaceGust"))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tqSpan_ct.txt", tqSpan_ct)

# Dynamic analysis of the Pazy wing encountering a continuous, 2-dimensional gust
include("../examples/PazyWingContinuous2DSpaceGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingContinuous2DSpaceGust"))
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tqSpan_ct.txt", tqSpan_ct)

# Dynamic analysis of the Pazy wing encountering a DARPA gust
include("../examples/PazyWingDARPAGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingDARPAGust"))
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tqSpan_ct.txt", tqSpan_ct)

# Flutter analysis of the Pazy wing
include("../examples/PazyWingFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutter"))
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOnsetSpeedsOfMode.txt", flutterOnsetSpeedsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOnsetFreqsOfMode.txt", flutterOnsetFreqsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOnsetDispOfMode.txt", flutterOnsetDispOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOffsetSpeedsOfMode.txt", flutterOffsetSpeedsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOffsetFreqsOfMode.txt", flutterOffsetFreqsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOffsetDispOfMode.txt", flutterOffsetDispOfMode)

# Flutter and divergence analysis of the Pazy wing
include("../examples/PazyWingFlutterAndDivergence.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutterAndDivergence"))
writedlm("test/newTestDataGenerators/PazyWingFlutterAndDivergence/flutterOnsetSpeedsOfMode.txt", flutterOnsetSpeedsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutterAndDivergence/flutterOnsetFreqsOfMode.txt", flutterOnsetFreqsOfMode)

# Flutter analysis of the Pazy wing with varying root pitch angle
include("../examples/PazyWingFlutterPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutterPitchRange"))
for (i,θ) in enumerate(θRange)
    writedlm(string("test/newTestDataGenerators/PazyWingFlutterPitchRange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/PazyWingFlutterPitchRange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Flutter analysis of the Pazy wing with varying tip mass positions
include("../examples/PazyWingFlutterTipMassRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutterTipMassRange"))
for i=1:3
    writedlm(string("test/newTestDataGenerators/PazyWingFlutterTipMassRange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/PazyWingFlutterTipMassRange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Modal analyses of the Pazy wing in horizontal and vertical positions
include("../examples/PazyWingModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingModal"))
writedlm("test/newTestDataGenerators/PazyWingModal/freqs.txt", freqs)

# Dynamic analysis of the Pazy wing encountering a one-minus-cosine gust
include("../examples/PazyWingOMCGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingOMCGust"))
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tqSpan_ct.txt", tqSpan_ct)

# Steady analysis of the Pazy wing with varying root pitch angle
include("../examples/PazyWingPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingPitchRange"))
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_AoA.txt", tip_AoA)
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_OOP.txt", tip_OOP)
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_IP.txt", tip_IP)
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_twist.txt", tip_twist)

# Dynamic analysis of the Pazy wing with a tip impulse force
include("../examples/PazyWingTipImpulse.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingTipImpulse"))
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tqSpan_ct.txt", tqSpan_ct)

# Static analysis of the coupled torsion-bending test of the Pazy wing
include("../examples/PazyWingTorsionTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingTorsionTest"))
writedlm("test/newTestDataGenerators/PazyWingTorsionTest/tip_OOP.txt", tip_OOP)
writedlm("test/newTestDataGenerators/PazyWingTorsionTest/tip_twist.txt", tip_twist)

# Dynamic analysis of a pendulum released from rest
include("../examples/pendulum.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pendulum"))
writedlm("test/newTestDataGenerators/pendulum/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/pendulum/u3_tip.txt", u3_tip)

# Static analysis of a mid-loaded arch pinned at one end and clamped at the other
include("../examples/pinnedClampedArch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pinnedClampedArch"))
writedlm("test/newTestDataGenerators/pinnedClampedArch/u1.txt", u1_atForce)
writedlm("test/newTestDataGenerators/pinnedClampedArch/u3.txt", u3_atForce)

# Modal analysis of a beam pinned at one end and transversely springed at the other
include("../examples/pinnedSpringedBeamEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pinnedSpringedBeamEigen"))
writedlm("test/newTestDataGenerators/pinnedSpringedBeamEigen/freqs.txt", freqs)

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

# Trim analysis (reaction loads check) of a right-angled frame
include("../examples/rightAngledFrameTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rightAngledFrameTrim"))
writedlm("test/newTestDataGenerators/rightAngledFrameTrim/balanceHorizontalForce.txt", balanceHorizontalForce)
writedlm("test/newTestDataGenerators/rightAngledFrameTrim/balanceVerticalForce.txt", balanceVerticalForce)

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

# Sharp-edged gust response of an airfoil section at several pitch angles
include("../examples/SEgustTests.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SEgustTests"))
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
            relPath = string("/test/newTestDataGenerators/SEgustTests/SEgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase)
            absPath = string(pwd(),relPath)
            mkpath(absPath)
            writedlm(string("test/newTestDataGenerators/SEgustTests/SEgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/tau",i,j,k,".txt"), τ[i,j,k])
            writedlm(string("test/newTestDataGenerators/SEgustTests/SEgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dcl",i,j,k,".txt"), Δcl[i,j,k])
            writedlm(string("test/newTestDataGenerators/SEgustTests/SEgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dclRef",i,j,k,".txt"), ΔclRef[i,j,k])
        end
    end
end

# Flutter analysis of the sixteen-meter-wing
include("../examples/SMWFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutter"))
writedlm("test/newTestDataGenerators/SMWFlutter/x1_def.txt",x1_def[end]/L)
writedlm("test/newTestDataGenerators/SMWFlutter/x3_def.txt",x3_def[end]/L)
writedlm("test/newTestDataGenerators/SMWFlutter/alpha.txt",α_of_x1[end]*180/pi)
for mode in 1:nModes
    writedlm(string("test/newTestDataGenerators/SMWFlutter/freqsMode",mode,".txt"), modeFrequencies[mode])
    writedlm(string("test/newTestDataGenerators/SMWFlutter/dampsMode",mode,".txt"), modeDampings[mode])
end

# Flutter boundary analysis of the sixteen-meter-wing as a function of the pitch angle
include("../examples/SMWFlutterPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterPitchRange"))
for i in 1:nModes
    writedlm(string("test/newTestDataGenerators/SMWFlutterPitchRange/freqsMode",i,".txt"),modeFrequencies[i])
    writedlm(string("test/newTestDataGenerators/SMWFlutterPitchRange/dampsMode",i,".txt"),modeDampings[i])
end

# Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with tip load as the varying parameter
include("../examples/SMWFlutterPrecurvatureRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterPrecurvatureRange"))
for i in eachindex(kRange)
    writedlm(string("test/newTestDataGenerators/SMWFlutterPrecurvatureRange/flutterSpeedk",i,".txt"),flutterSpeed[i,:])
end

# Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with root angle as the varying parameter
include("../examples/SMWFlutterPrecurvatureRange2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterPrecurvatureRange2"))
for (ki,k) in enumerate(kRange)
    for (i,θ) in enumerate(θRange)
        writedlm(string("test/newTestDataGenerators/SMWFlutterPrecurvatureRange2/freqs_k",ki,"th",i,".txt"),hcat(freqs[ki,i,:]...)')
        writedlm(string("test/newTestDataGenerators/SMWFlutterPrecurvatureRange2/damps_k",ki,"th",i,".txt"),hcat(damps[ki,i,:]...)')
    end
end

# Flutter boundary analysis of the sixteen-meter-wing as a function of the bending-torsion coupling factor
include("../examples/SMWFlutterStructuralCouplingRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterStructuralCouplingRange"))
writedlm("test/newTestDataGenerators/SMWFlutterStructuralCouplingRange/flutterSpeed.txt", flutterSpeed)

# Flutter boundary analysis of the sixteen-meter-wing as a function of the tip displacement
include("../examples/SMWFlutterTipDispRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterTipDispRange"))
writedlm("test/newTestDataGenerators/SMWFlutterTipDispRange/flutterSpeed.txt", flutterSpeed)
writedlm("test/newTestDataGenerators/SMWFlutterTipDispRange/flutterFreq.txt", flutterFreq)

# Linear flutter analysis of the sixteen-meter-wing
include("../examples/SMWLinearFlutter.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWLinearFlutter"))
writedlm("test/newTestDataGenerators/SMWLinearFlutter/flutterSpeed.txt", flutterSpeed)
writedlm("test/newTestDataGenerators/SMWLinearFlutter/flutterFreq.txt", flutterFreq)

# Modal analysis of the sixteen-meter-wing
include("../examples/SMWModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWModal"))
writedlm("test/newTestDataGenerators/SMWModal/freqs.txt", freqs)

# Steady aeroelastic analysis of the sixteen-meter-wing
include("../examples/SMWSteady.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWSteady"))
writedlm("test/newTestDataGenerators/SMWSteady/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/SMWSteady/tip_twist.txt", tip_twist)

# Dynamic analysis of the spin-up maneuver of a robot arm
include("../examples/spinupRobotArm.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/spinupRobotArm"))
writedlm("test/newTestDataGenerators/spinupRobotArm/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/spinupRobotArm/u2_tip.txt", u2_tip)
writedlm("test/newTestDataGenerators/spinupRobotArm/theta3_root.txt", θ3_root)

# Static analysis of a L-frame with a doubly-attached spring and tip load
include("../examples/springedLFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/springedLFrame"))
writedlm("test/newTestDataGenerators/springedLFrame/u3_b.txt", u3_b)

# Eigen-analysis of a straight rotor
include("../examples/straightRotor.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/straightRotor"))
writedlm("test/newTestDataGenerators/straightRotor/freqs.txt", numFreqs)

# Eigen-analysis of a swept-tip rotor
include("../examples/sweptTipRotor.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptTipRotor"))
writedlm("test/newTestDataGenerators/sweptTipRotor/freqs.txt", numFreqs[end,end])

# Static analysis of a semi-circular arch with tangential follower force
include("../examples/tangentiallyForcedArch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tangentiallyForcedArch"))
writedlm("test/newTestDataGenerators/tangentiallyForcedArch/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/tangentiallyForcedArch/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/tangentiallyForcedArch/tip_angle.txt", tip_angle)

# Modal analysis of a cantilevered tapered beam
include("../examples/taperedBeamEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/taperedBeamEigen"))
writedlm("test/newTestDataGenerators/taperedBeamEigen/freqs.txt", freqs)

# Steady aeroelastic analysis of the Tang&Dowell wing at varying airspeed
include("../examples/TDWingAirspeedRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/TDWingAirspeedRange"))
writedlm("test/newTestDataGenerators/TDWingAirspeedRange/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/TDWingAirspeedRange/tip_twist.txt", tip_twist)

# Modal analysis of the Tang&Dowell wing at varying pitch angles
include("../examples/TDWingPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/TDWingPitchRange"))
writedlm("test/newTestDataGenerators/TDWingPitchRange/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/TDWingPitchRange/tip_u2.txt", tip_u2)
writedlm("test/newTestDataGenerators/TDWingPitchRange/tip_twist.txt", tip_twist)
writedlm("test/newTestDataGenerators/TDWingPitchRange/freqs.txt", freqs)

# Dynamic analysis of an airfoil section sinusoidally surging (facing a time-varying freestream)
include("../examples/timeVaryingFreestream.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/timeVaryingFreestream"))
writedlm("test/newTestDataGenerators/timeVaryingFreestream/cn.txt", cn[end])
writedlm("test/newTestDataGenerators/timeVaryingFreestream/cm.txt", cm[end])

# Dynamic analysis of an airfoil section sinusoidally pitching and surging (facing a time-varying freestream)
include("../examples/timeVaryingFreestreamAndPitch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/timeVaryingFreestreamAndPitch"))
writedlm("test/newTestDataGenerators/timeVaryingFreestreamAndPitch/cn.txt", cn[end])
writedlm("test/newTestDataGenerators/timeVaryingFreestreamAndPitch/cm.txt", cm[end])

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

# Trim analysis of a cantilever with tip force
include("../examples/tipLoadedCantileverTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipLoadedCantileverTrim"))
writedlm("test/newTestDataGenerators/tipLoadedCantileverTrim/u3.txt", u3)
writedlm("test/newTestDataGenerators/tipLoadedCantileverTrim/F3.txt", F3)
writedlm("test/newTestDataGenerators/tipLoadedCantileverTrim/M2.txt", M2)

# Static analysis of a cantilever with tip moment
include("../examples/tipMomentCantilever.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipMomentCantilever"))
writedlm("test/newTestDataGenerators/tipMomentCantilever/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/tipMomentCantilever/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/tipMomentCantilever/tip_angle.txt", tip_angle)

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

# Eigen analysis of a 2-story frame
include("../examples/twoStoryFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/twoStoryFrame"))
writedlm("test/newTestDataGenerators/twoStoryFrame/freqs.txt", freqs)
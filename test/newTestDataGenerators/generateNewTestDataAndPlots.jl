using DelimitedFiles

# Dynamic analysis of an airfoils with harmonic flap deflection profile
include("../plotGenerators/flapOscillationPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flapOscillation"))
writedlm("test/newTestDataGenerators/flapOscillation/cn.txt", cn)
writedlm("test/newTestDataGenerators/flapOscillation/cm.txt", cm)

# Dynamic analysis of two airfoils with linked harmonic flap deflection profiles
include("../plotGenerators/flapOscillationLinkedPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flapOscillationLinked"))
writedlm("test/newTestDataGenerators/flapOscillationLinked/cnMaster.txt", cnMaster)
writedlm("test/newTestDataGenerators/flapOscillationLinked/cmMaster.txt", cmMaster)
writedlm("test/newTestDataGenerators/flapOscillationLinked/cnSlave.txt", cnSlave)
writedlm("test/newTestDataGenerators/flapOscillationLinked/cmSlave.txt", cmSlave)

# One-minus-cosine gust response of an airfoil section at several pitch angles
include("../plotGenerators/OMCgustTestsPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/OMCgustTests"))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
        for (k,testCase) in enumerate(tests)
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

# One-minus-cosine gust response of an airfoil section at several pitch angles
include("../examples/OMCgustTests2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/OMCgustTests2"))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
        for (k,testCase) in enumerate(tests)
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
            relPath = string("/test/newTestDataGenerators/OMCgustTests2/OMCgustTests2_",aeroSolverName,"_",gustSolverName,"_test",testCase)
            absPath = string(pwd(),relPath)
            mkpath(absPath)
            writedlm(string("test/newTestDataGenerators/OMCgustTests2/OMCgustTests2_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/tau",i,j,k,".txt"), τ[i,j,k])
            writedlm(string("test/newTestDataGenerators/OMCgustTests2/OMCgustTests2_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dcl",i,j,k,".txt"), Δcl[i,j,k])
            writedlm(string("test/newTestDataGenerators/OMCgustTests2/OMCgustTests2_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dclRef",i,j,k,".txt"), ΔclRef[i,j,k])
        end
    end
end

# Dynamic analysis of a harmonically pitching airfoil, using the dynamic stall model
include("../plotGenerators/pitchingAirfoilDSModelTestPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pitchingAirfoilDSModelTest"))
writedlm("test/newTestDataGenerators/pitchingAirfoilDSModelTest/"*frameString*"_cn.txt", cn)
writedlm("test/newTestDataGenerators/pitchingAirfoilDSModelTest/"*frameString*"_cm.txt", cm)
writedlm("test/newTestDataGenerators/pitchingAirfoilDSModelTest/"*frameString*"_ct.txt", ct)

# Dynamic analysis of a harmonically plunging airfoil, using the dynamic stall model
include("../plotGenerators/plungingAirfoilDSModelTestPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/plungingAirfoilDSModelTest"))
writedlm("test/newTestDataGenerators/plungingAirfoilDSModelTest/cn.txt", cn)
writedlm("test/newTestDataGenerators/plungingAirfoilDSModelTest/cm.txt", cm)

# Sharp-edged gust response of an airfoil section at several pitch angles
include("../plotGenerators/SEgustTestsPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SEgustTests"))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
        for (k,testCase) in enumerate(tests)
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

# Dynamic analysis of an airfoil section sinusoidally surging (facing a time-varying freestream)
include("../plotGenerators/timeVaryingFreestreamPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/timeVaryingFreestream"))
writedlm("test/newTestDataGenerators/timeVaryingFreestream/cn.txt", cn[end])
writedlm("test/newTestDataGenerators/timeVaryingFreestream/cm.txt", cm[end])

# Dynamic analysis of an airfoil section sinusoidally pitching and surging (facing a time-varying freestream)
include("../plotGenerators/timeVaryingFreestreamAndPitchPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/timeVaryingFreestreamAndPitch"))
writedlm("test/newTestDataGenerators/timeVaryingFreestreamAndPitch/cn.txt", cn[end])
writedlm("test/newTestDataGenerators/timeVaryingFreestreamAndPitch/cm.txt", cm[end])

# Dynamic analysis of the Blended-Wing-Body vehicle undergoing a checked pitch maneuver
include("../plotGenerators/BWBcheckedPitchManeuverPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/BWBcheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/BWBcheckedPitchManeuver/rootAoA.txt", rootAoA)
writedlm("test/newTestDataGenerators/BWBcheckedPitchManeuver/Deltau3.txt", Δu3)

# Dynamic analysis of the conventional HALE aircraft undergoing a checked pitch maneuver
include("../plotGenerators/conventionalHALECheckedPitchManeuverPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALECheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/conventionalHALECheckedPitchManeuver/wingAoA.txt", wingAoA)
writedlm("test/newTestDataGenerators/conventionalHALECheckedPitchManeuver/Deltau3.txt", Δu3)

# Dynamic analysis of the conventional HALE aircraft undergoing a coordinated turn maneuver
include("../plotGenerators/conventionalHALECheckedRollManeuverPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALECheckedRollManeuver"))
writedlm("test/newTestDataGenerators/conventionalHALECheckedRollManeuver/wingAoA.txt", wingAoA)
writedlm("test/newTestDataGenerators/conventionalHALECheckedRollManeuver/Deltau3.txt", Δu3)

# Dynamic analysis of the baseline Healy FFWT wing under a series of one-minus-cosine gusts
include("../plotGenerators/HealyBaselineFFWTOMCGustFloatingPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTOMCGustFloating"))
for (i,ω) in enumerate(ωRange)
    writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTOMCGustFloating/M2root_omega",ω,".txt"),M2root[i])
    writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTOMCGustFloating/fold_omega",ω,".txt"),ϕ[i])
end

# Dynamic analysis of the baseline Healy FFWT wing for a specific one-minus-cosine gust
include("../examples/HealyBaselineFFWTOMCGustFloating2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTOMCGustFloating2"))
writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTOMCGustFloating2/M2root.txt"),M2root)
writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTOMCGustFloating2/fold.txt"),ϕ)

# Dynamic analysis of the Helios flying-wing undergoing a checked pitch maneuver
include("../plotGenerators/heliosCheckedPitchManeuverPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosCheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/heliosCheckedPitchManeuver/rootAoA.txt", rootAoA)
writedlm("test/newTestDataGenerators/heliosCheckedPitchManeuver/Deltau3.txt", Δu3)

# Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over time
include("../plotGenerators/PazyWingContinuous1DGustPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingContinuous1DGust"))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tqSpan_ct.txt", tqSpan_ct)

# Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over space
include("../plotGenerators/PazyWingContinuous1DSpaceGustPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingContinuous1DSpaceGust"))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tqSpan_ct.txt", tqSpan_ct)

# Dynamic analysis of the Pazy wing encountering a continuous, 2-dimensional gust
include("../plotGenerators/PazyWingContinuous2DSpaceGustPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingContinuous2DSpaceGust"))
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tqSpan_ct.txt", tqSpan_ct)

# Dynamic analysis of the Pazy wing encountering a DARPA gust
include("../plotGenerators/PazyWingDARPAGustPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingDARPAGust"))
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tqSpan_ct.txt", tqSpan_ct)

# Dynamic analysis of the Pazy wing encountering a one-minus-cosine gust
include("../plotGenerators/PazyWingOMCGustPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingOMCGust"))
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tqSpan_ct.txt", tqSpan_ct)

# Dynamic analysis of the Pazy wing with a tip impulse force
include("../plotGenerators/PazyWingTipImpulsePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingTipImpulse"))
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tqSpan_ct.txt", tqSpan_ct)

# Dynamic analysis of the axial vibration of a beam under a traction force applied suddenly
include("../plotGenerators/axialTractionCantileverPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/axialTractionCantilever"))
writedlm("test/newTestDataGenerators/axialTractionCantilever/u1_08.txt", u1_08)
writedlm("test/newTestDataGenerators/axialTractionCantilever/u1_10.txt", u1_10)

# Dynamic analysis of the free response of a beam clamped at both ends and subjected to an initial displacement profile
include("../plotGenerators/biclampedBeamPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/biclampedBeam"))
writedlm("test/newTestDataGenerators/biclampedBeam/u3_mid.txt", u3_mid)
writedlm("test/newTestDataGenerators/biclampedBeam/V3_mid.txt", V3_mid)
writedlm("test/newTestDataGenerators/biclampedBeam/Vdot3_mid.txt", Vdot3_mid)
writedlm("test/newTestDataGenerators/biclampedBeam/theta2_quarter.txt", θ2_quarter)
writedlm("test/newTestDataGenerators/biclampedBeam/Omega2_quarter.txt", Ω2_quarter)
writedlm("test/newTestDataGenerators/biclampedBeam/Omegadot2_quarter.txt", Ωdot2_quarter)

# Dynamic analysis of a composite cantilever beam under a tip sinusoidal load
include("../plotGenerators/compositeCantileverPlotGenerator.jl")
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
include("../plotGenerators/curvedCantileverDynamicFollowerPlotGenerator.jl")
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
include("../plotGenerators/doublePendulumPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/doublePendulum"))
writedlm("test/newTestDataGenerators/doublePendulum/u1_hinge.txt", u1_hinge)
writedlm("test/newTestDataGenerators/doublePendulum/u3_hinge.txt", u3_hinge)
writedlm("test/newTestDataGenerators/doublePendulum/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/doublePendulum/u3_tip.txt", u3_tip)

# Dynamic analysis of a right-angled frame subjected to an out-of-plane force
include("../plotGenerators/elbowFramePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/elbowFrame"))
writedlm("test/newTestDataGenerators/elbowFrame/u3_elbow.txt", u3_elbow)
writedlm("test/newTestDataGenerators/elbowFrame/u3_tip.txt", u3_tip)

# Dynamic analysis of a very flexible beam subjected to loads yielding two-dimensional motion
include("../plotGenerators/flyingFlexibleBeam2DPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flyingFlexibleBeam2D"))
writedlm("test/newTestDataGenerators/flyingFlexibleBeam2D/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/flyingFlexibleBeam2D/u3_tip.txt", u3_tip)

# Dynamic analysis of a hinged beam in free flight
include("../plotGenerators/flyingScissorsPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flyingScissors"))
writedlm("test/newTestDataGenerators/flyingScissors/u1_tipA.txt", u1_tipA)
writedlm("test/newTestDataGenerators/flyingScissors/u3_tipA.txt", u3_tipA)
writedlm("test/newTestDataGenerators/flyingScissors/u1_tipB.txt", u1_tipB)
writedlm("test/newTestDataGenerators/flyingScissors/u3_tipB.txt", u3_tipB)
writedlm("test/newTestDataGenerators/flyingScissors/u1_hinge.txt", u1_hinge)
writedlm("test/newTestDataGenerators/flyingScissors/u3_hinge.txt", u3_hinge)

# Dynamic analysis of a very flexible beam subjected to loads yielding two-dimensional motion
include("../plotGenerators/flyingSpaghetti2DPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flyingSpaghetti2D"))
writedlm("test/newTestDataGenerators/flyingSpaghetti2D/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/flyingSpaghetti2D/u3_tip.txt", u3_tip)

# Dynamic analysis of a very flexible beam subjected to loads yielding tri-dimensional motion
include("../plotGenerators/flyingSpaghetti3DPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flyingSpaghetti3D"))
writedlm("test/newTestDataGenerators/flyingSpaghetti3D/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/flyingSpaghetti3D/u2_tip.txt", u2_tip)
writedlm("test/newTestDataGenerators/flyingSpaghetti3D/u3_tip.txt", u3_tip)
writedlm("test/newTestDataGenerators/flyingSpaghetti3D/theta_tip.txt", θ_tip)

# Dynamic analysis of the free response of a beam subjected to initial displacement and velocity profiles
include("../plotGenerators/initialDispAndVelBeamPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/initialDispAndVelBeam"))
writedlm("test/newTestDataGenerators/initialDispAndVelBeam/u3_quarter.txt", u3_quarter)
writedlm("test/newTestDataGenerators/initialDispAndVelBeam/V3_quarter.txt", V3_quarter)
writedlm("test/newTestDataGenerators/initialDispAndVelBeam/Vdot3_quarter.txt", Vdot3_quarter)
writedlm("test/newTestDataGenerators/initialDispAndVelBeam/theta2_root.txt", θ2_root)
writedlm("test/newTestDataGenerators/initialDispAndVelBeam/Omega2_mid.txt", Ω2_mid)
writedlm("test/newTestDataGenerators/initialDispAndVelBeam/Omegadot2_mid.txt", Ωdot2_mid)

# Dynamic analysis of the free response of a beam subjected to an initial displacement profile
include("../plotGenerators/initialDisplacementBeamPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/initialDisplacementBeam"))
writedlm("test/newTestDataGenerators/initialDisplacementBeam/u3_quarter.txt", u3_quarter)
writedlm("test/newTestDataGenerators/initialDisplacementBeam/V3_quarter.txt", V3_quarter)
writedlm("test/newTestDataGenerators/initialDisplacementBeam/Vdot3_quarter.txt", Vdot3_quarter)
writedlm("test/newTestDataGenerators/initialDisplacementBeam/theta2_root.txt", θ2_root)
writedlm("test/newTestDataGenerators/initialDisplacementBeam/Omega2_mid.txt", Ω2_mid)
writedlm("test/newTestDataGenerators/initialDisplacementBeam/Omegadot2_mid.txt", Ωdot2_mid)

# Dynamic analysis of the free response of a beam subjected to an initial velocity profile
include("../plotGenerators/initialVelocityBeamPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/initialVelocityBeam"))
writedlm("test/newTestDataGenerators/initialVelocityBeam/u3_quarter.txt", u3_quarter)
writedlm("test/newTestDataGenerators/initialVelocityBeam/V3_quarter.txt", V3_quarter)
writedlm("test/newTestDataGenerators/initialVelocityBeam/Vdot3_quarter.txt", Vdot3_quarter)
writedlm("test/newTestDataGenerators/initialVelocityBeam/theta2_root.txt", θ2_root)
writedlm("test/newTestDataGenerators/initialVelocityBeam/Omega2_mid.txt", Ω2_mid)
writedlm("test/newTestDataGenerators/initialVelocityBeam/Omegadot2_mid.txt", Ωdot2_mid)

# Dynamic analysis of joined beams under load
include("../plotGenerators/joinedBeamsPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/joinedBeams"))
writedlm("test/newTestDataGenerators/joinedBeams/u1.txt", u1)
writedlm("test/newTestDataGenerators/joinedBeams/u2.txt", u2)
writedlm("test/newTestDataGenerators/joinedBeams/u3.txt", u3)

# Dynamic analysis of a pinned robot arm driven by a couple moment
include("../plotGenerators/momentDrivenRobotArmPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/momentDrivenRobotArm"))
writedlm("test/newTestDataGenerators/momentDrivenRobotArm/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/momentDrivenRobotArm/u3_tip.txt", u3_tip)

# Dynamic analysis of a pendulum released from rest
include("../plotGenerators/pendulumPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pendulum"))
writedlm("test/newTestDataGenerators/pendulum/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/pendulum/u3_tip.txt", u3_tip)

# Dynamic analysis of a clamped beam with root sinusoidal oscillation (ω = 9 Hz ≈ 2nd bending mode)
include("../plotGenerators/rootExcitationBeam1PlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rootExcitationBeam1"))
writedlm("test/newTestDataGenerators/rootExcitationBeam1/u3b_root.txt", u3b_root)
writedlm("test/newTestDataGenerators/rootExcitationBeam1/u3b_tip.txt", u3b_tip)
writedlm("test/newTestDataGenerators/rootExcitationBeam1/V3_root.txt", V3_root)
writedlm("test/newTestDataGenerators/rootExcitationBeam1/V3_tip.txt", V3_tip)

# Dynamic analysis of a clamped beam with root sinusoidal oscillation (ω = 32 Hz ≈ 3rd bending mode)
include("../plotGenerators/rootExcitationBeam2PlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rootExcitationBeam2"))
writedlm("test/newTestDataGenerators/rootExcitationBeam2/u3b_root.txt", u3b_root)
writedlm("test/newTestDataGenerators/rootExcitationBeam2/u3b_tip.txt", u3b_tip)
writedlm("test/newTestDataGenerators/rootExcitationBeam2/V3_root.txt", V3_root)
writedlm("test/newTestDataGenerators/rootExcitationBeam2/V3_tip.txt", V3_tip)

# Dynamic analysis of a rotary shaft with specified rotation
include("../plotGenerators/rotaryShaftPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rotaryShaft"))
writedlm("test/newTestDataGenerators/rotaryShaft/pNum.txt", pNum)
writedlm("test/newTestDataGenerators/rotaryShaft/pdotNum.txt", pdotNum)
writedlm("test/newTestDataGenerators/rotaryShaft/OmegaNum.txt", ΩNum)
writedlm("test/newTestDataGenerators/rotaryShaft/OmegadotNum.txt", ΩdotNum)
writedlm("test/newTestDataGenerators/rotaryShaft/MNum.txt", MNum)

# Dynamic analysis of an articulated robot arm driven by specified rotation
include("../plotGenerators/rotationDrivenArticulatedRobotArmPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rotationDrivenArticulatedRobotArm"))
writedlm("test/newTestDataGenerators/rotationDrivenArticulatedRobotArm/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/rotationDrivenArticulatedRobotArm/u3_tip.txt", u3_tip)
writedlm("test/newTestDataGenerators/rotationDrivenArticulatedRobotArm/u1_hinge.txt", u1_hinge)
writedlm("test/newTestDataGenerators/rotationDrivenArticulatedRobotArm/u3_hinge.txt", u3_hinge)

# Dynamic analysis of a robot arm driven by specified rotation
include("../plotGenerators/rotationDrivenRobotArmPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rotationDrivenRobotArm"))
writedlm("test/newTestDataGenerators/rotationDrivenRobotArm/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/rotationDrivenRobotArm/u3_tip.txt", u3_tip)

# Dynamic analysis of the spin-up maneuver of a robot arm
include("../plotGenerators/spinupRobotArmPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/spinupRobotArm"))
writedlm("test/newTestDataGenerators/spinupRobotArm/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/spinupRobotArm/u2_tip.txt", u2_tip)
writedlm("test/newTestDataGenerators/spinupRobotArm/theta3_root.txt", θ3_root)

# Dynamic analysis of a pendulum with tip mass
include("../plotGenerators/tipPendulumPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipPendulum"))
writedlm("test/newTestDataGenerators/tipPendulum/u1_tip.txt", u1_tip)
writedlm("test/newTestDataGenerators/tipPendulum/u3_tip.txt", u3_tip)

# Dynamic analysis of a cantilever with tip sinusoidal force
include("../plotGenerators/tipSineLoadedCantileverPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipSineLoadedCantilever"))
writedlm("test/newTestDataGenerators/tipSineLoadedCantilever/u3_tip.txt", u3_tip)
writedlm("test/newTestDataGenerators/tipSineLoadedCantilever/F3_root.txt", F3_root)
writedlm("test/newTestDataGenerators/tipSineLoadedCantilever/M2_root.txt", M2_root)

# Flutter analysis of the Blended-Wing-Body flying wing
include("../plotGenerators/BWBflutterPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/BWBflutter"))
writedlm("test/newTestDataGenerators/BWBflutter/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/BWBflutter/damps.txt", damps)

# Flutter analysis of the conventional HALE aircraft in free flight with structural stiffness as the varying parameter
include("../plotGenerators/conventionalHALELambdaRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALELambdaRange"))
writedlm("test/newTestDataGenerators/conventionalHALELambdaRange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/conventionalHALELambdaRange/damps.txt", damps)

# Flutter analysis of the conventional HALE aircraft in free flight with airspeed and structural stiffness as the varying parameters
include("../plotGenerators/conventionalHALELURangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALELURange"))
for (i,λ) in enumerate(λRange)
    writedlm(string("test/newTestDataGenerators/conventionalHALELURange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/conventionalHALELURange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Flutter analysis of the conventional HALE aircraft in free flight with airspeed as the varying parameter
include("../plotGenerators/conventionalHALEURangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALEURange"))
writedlm("test/newTestDataGenerators/conventionalHALEURange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/conventionalHALEURange/damps.txt", damps)

# Flutter analysis of the Goland wing
include("../plotGenerators/GolandWingFlutterPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/GolandWingFlutter"))
writedlm("test/newTestDataGenerators/GolandWingFlutter/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/GolandWingFlutter/damps.txt", damps)

# Flutter analysis of the baseline Healy free FFWT wing with tip loss, root pitch and airspeed as the varying parameters
include("../plotGenerators/HealyBaselineFFWTfreeFlutterAoARangeURangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTfreeFlutterAoARangeURange"))
for (i,τ) in enumerate(τRange)
    for (j,θ) in enumerate(θRange)
        writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTfreeFlutterAoARangeURange/freqs",i,j,".txt"), hcat(freqs[i,j,:]...)')
        writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTfreeFlutterAoARangeURange/damps",i,j,".txt"), hcat(damps[i,j,:]...)')
    end
end

# Flutter analysis of the baseline Healy free FFWT wing with flare angle and airspeed as the varying parameters
include("../plotGenerators/HealyBaselineFFWTfreeFlutterFlareRangeURangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTfreeFlutterFlareRangeURange"))
for (i,Λ) in enumerate(ΛRange)
    writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTfreeFlutterFlareRangeURange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTfreeFlutterFlareRangeURange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Flutter analysis of the baseline Healy locked FFWT wing with airspeed as the varying parameters
include("../plotGenerators/HealyBaselineFFWTlockedFlutterURangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTlockedFlutterURange"))
writedlm("test/newTestDataGenerators/HealyBaselineFFWTlockedFlutterURange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/HealyBaselineFFWTlockedFlutterURange/damps.txt", damps)

# Flutter analysis of the Healy FFWT wing with sideslip angle as the varying parameter
include("../plotGenerators/HealySideslipFFWTflutterSideslipRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealySideslipFFWTflutterSideslipRange"))
writedlm("test/newTestDataGenerators/HealySideslipFFWTflutterSideslipRange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/HealySideslipFFWTflutterSideslipRange/damps.txt", damps)

# Flutter analysis of the Healy FFWT wing with airspeed as the varying parameter
include("../plotGenerators/HealySideslipFFWTflutterURangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealySideslipFFWTflutterURange"))
writedlm("test/newTestDataGenerators/HealySideslipFFWTflutterURange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/HealySideslipFFWTflutterURange/damps.txt", damps)

# Flutter analysis of the Helios flying-wing in free flight with payload and structural stiffness as the varying parameters
include("../plotGenerators/heliosFlutterPLambdaRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosFlutterPLambdaRange"))
for (i,λ) in enumerate(λRange)
    writedlm(string("test/newTestDataGenerators/heliosFlutterPLambdaRange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/heliosFlutterPLambdaRange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Flutter analysis of the Helios flying-wing in free flight with payload as the varying parameter
include("../plotGenerators/heliosFlutterPRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosFlutterPRange"))
writedlm("test/newTestDataGenerators/heliosFlutterPRange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/heliosFlutterPRange/damps.txt", damps)

# Flutter analysis of the Helios flying-wing in free flight with airspeed as the varying parameter
include("../plotGenerators/heliosFlutterURangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosFlutterURange"))
writedlm("test/newTestDataGenerators/heliosFlutterURange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/heliosFlutterURange/damps.txt", damps)

# Flutter analysis of the wing of the Helios flying-wing
include("../plotGenerators/heliosWingFlutterPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosWingFlutter"))
writedlm("test/newTestDataGenerators/heliosWingFlutter/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/heliosWingFlutter/damps.txt", damps)

# Flutter analysis of the Pazy FFWT wing with airspeed as the varying parameter
include("../plotGenerators/PazyFFWTflutterURangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTflutterURange"))
writedlm("test/newTestDataGenerators/PazyFFWTflutterURange/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/PazyFFWTflutterURange/damps.txt", damps)
writedlm("test/newTestDataGenerators/PazyFFWTflutterURange/phiHinge.txt", ϕHinge)

# Flutter analysis of the Pazy wing
include("../plotGenerators/PazyWingFlutterPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutter"))
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOnsetSpeedsOfMode.txt", flutterOnsetSpeedsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOnsetFreqsOfMode.txt", flutterOnsetFreqsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOnsetDispOfMode.txt", flutterOnsetDispOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOffsetSpeedsOfMode.txt", flutterOffsetSpeedsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOffsetFreqsOfMode.txt", flutterOffsetFreqsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutter/flutterOffsetDispOfMode.txt", flutterOffsetDispOfMode)

# Flutter and divergence analysis of the Pazy wing
include("../plotGenerators/PazyWingFlutterAndDivergencePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutterAndDivergence"))
writedlm("test/newTestDataGenerators/PazyWingFlutterAndDivergence/flutterOnsetSpeedsOfMode.txt", flutterOnsetSpeedsOfMode)
writedlm("test/newTestDataGenerators/PazyWingFlutterAndDivergence/flutterOnsetFreqsOfMode.txt", flutterOnsetFreqsOfMode)

# Flutter analysis of the Pazy wing with varying root pitch angle
include("../plotGenerators/PazyWingFlutterPitchRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutterPitchRange"))
for (i,θ) in enumerate(θRange)
    writedlm(string("test/newTestDataGenerators/PazyWingFlutterPitchRange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/PazyWingFlutterPitchRange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Flutter analysis of the Pazy wing with varying tip mass positions
include("../plotGenerators/PazyWingFlutterTipMassRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFlutterTipMassRange"))
for i=1:3
    writedlm(string("test/newTestDataGenerators/PazyWingFlutterTipMassRange/freqs",i,".txt"), hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/PazyWingFlutterTipMassRange/damps",i,".txt"), hcat(damps[i,:]...)')
end

# Frequency analysis of the Pazy wing
include("../plotGenerators/PazyWingFreqsEvolutionPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingFreqsEvolution"))
writedlm("test/newTestDataGenerators/PazyWingFreqsEvolution/freqs.txt", freqs)

# Flutter analysis of the sixteen-meter-wing
include("../plotGenerators/SMWFlutterPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutter"))
writedlm("test/newTestDataGenerators/SMWFlutter/x1_def.txt",x1_def[end]/L)
writedlm("test/newTestDataGenerators/SMWFlutter/x3_def.txt",x3_def[end]/L)
writedlm("test/newTestDataGenerators/SMWFlutter/alpha.txt",α_of_x1[end]*180/pi)
for mode in 1:nModes
    writedlm(string("test/newTestDataGenerators/SMWFlutter/freqsMode",mode,".txt"), modeFrequencies[mode])
    writedlm(string("test/newTestDataGenerators/SMWFlutter/dampsMode",mode,".txt"), modeDampings[mode])
end

# Flutter boundary analysis of the sixteen-meter-wing as a function of the pitch angle
include("../plotGenerators/SMWFlutterPitchRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterPitchRange"))
for i in 1:nModes
    writedlm(string("test/newTestDataGenerators/SMWFlutterPitchRange/freqsMode",i,".txt"),modeFrequencies[i])
    writedlm(string("test/newTestDataGenerators/SMWFlutterPitchRange/dampsMode",i,".txt"),modeDampings[i])
end

# Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with tip load as the varying parameter
include("../plotGenerators/SMWFlutterPrecurvatureRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterPrecurvatureRange"))
for i in eachindex(kRange)
    writedlm(string("test/newTestDataGenerators/SMWFlutterPrecurvatureRange/flutterSpeedk",i,".txt"),flutterSpeed[i,:])
end

# Flutter boundary analysis of the sixteen-meter-wing as a function of the wing curvature with root angle as the varying parameter
include("../plotGenerators/SMWFlutterPrecurvatureRange2PlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterPrecurvatureRange2"))
for (ki,k) in enumerate(kRange)
    for (i,θ) in enumerate(θRange)
        writedlm(string("test/newTestDataGenerators/SMWFlutterPrecurvatureRange2/freqs_k",ki,"th",i,".txt"),hcat(freqs[ki,i,:]...)')
        writedlm(string("test/newTestDataGenerators/SMWFlutterPrecurvatureRange2/damps_k",ki,"th",i,".txt"),hcat(damps[ki,i,:]...)')
    end
end

# Flutter boundary analysis of the sixteen-meter-wing as a function of the bending-torsion coupling factor
include("../plotGenerators/SMWFlutterStructuralCouplingRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterStructuralCouplingRange"))
writedlm("test/newTestDataGenerators/SMWFlutterStructuralCouplingRange/flutterSpeed.txt", flutterSpeed)

# Flutter boundary analysis of the sixteen-meter-wing as a function of the tip displacement
include("../plotGenerators/SMWFlutterTipDispRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWFlutterTipDispRange"))
writedlm("test/newTestDataGenerators/SMWFlutterTipDispRange/flutterSpeed.txt", flutterSpeed)
writedlm("test/newTestDataGenerators/SMWFlutterTipDispRange/flutterFreq.txt", flutterFreq)

# Linear flutter analysis of the sixteen-meter-wing
include("../plotGenerators/SMWLinearFlutterPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWLinearFlutter"))
writedlm("test/newTestDataGenerators/SMWLinearFlutter/flutterSpeed.txt", flutterSpeed)
writedlm("test/newTestDataGenerators/SMWLinearFlutter/flutterFreq.txt", flutterFreq)

# Torsional divergence analysis of a straight wing
include("../plotGenerators/straightWingTorsionalDivergencePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/straightWingTorsionalDivergence"))
writedlm("test/newTestDataGenerators/straightWingTorsionalDivergence/damps.txt", dampingsNonOscillatory)

# Flutter analysis of the swept Pazy wing
include("../plotGenerators/sweptPazyFlutterPitchRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptPazyFlutterPitchRange"))
for (i,θ) in enumerate(θRange)
    writedlm(string("test/newTestDataGenerators/sweptPazyFlutterPitchRange/freqs_",i,".txt"),hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/sweptPazyFlutterPitchRange/damps_",i,".txt"),hcat(damps[i,:]...)')
end

# Flutter analysis of the undeformed swept Pazy wing
include("../plotGenerators/sweptPazyUndeformedFlutterSweepRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptPazyUndeformedFlutterSweepRange"))
for (i,Λ) in enumerate(ΛRange)
    writedlm(string("test/newTestDataGenerators/sweptPazyUndeformedFlutterSweepRange/freqs_",i,".txt"),hcat(freqs[i,:]...)')
    writedlm(string("test/newTestDataGenerators/sweptPazyUndeformedFlutterSweepRange/damps_",i,".txt"),hcat(damps[i,:]...)')
end

# Bending divergence analyses of swept wings
include("../plotGenerators/sweptWingBendingDivergenceSweepRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptWingBendingDivergenceSweepRange"))
for (i,Λ) in enumerate(ΛRange)
    writedlm(string("test/newTestDataGenerators/sweptWingBendingDivergenceSweepRange/damps_",i,".txt"),hcat(dampingsNonOscillatory[i,:]...)')
end

# Coupled bending-torsion divergence analyses of swept wings
include("../plotGenerators/sweptWingBendingTorsionalDivergenceSweepRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptWingBendingTorsionalDivergenceSweepRange"))
for (i,Λ) in enumerate(ΛRange)
    writedlm(string("test/newTestDataGenerators/sweptWingBendingTorsionalDivergenceSweepRange/damps_",i,".txt"),hcat(dampingsNonOscillatory[i,:]...)')
end

# Torsional divergence analyses of swept wings
include("../plotGenerators/sweptWingTorsionalDivergenceSweepRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptWingTorsionalDivergenceSweepRange"))
for (i,Λ) in enumerate(ΛRange)
    writedlm(string("test/newTestDataGenerators/sweptWingTorsionalDivergenceSweepRange/damps_",i,".txt"),hcat(dampingsNonOscillatory[i,:]...)')
end

# Torsional divergence analysis of a typical section
include("../plotGenerators/typicalSectionDivergencePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/typicalSectionDivergence"))
writedlm("test/newTestDataGenerators/typicalSectionDivergence/damps.txt", dampingsNonOscillatory)

# Flutter and divergence analysis of a typical section
include("../plotGenerators/typicalSectionFlutterAndDivergencePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/typicalSectionFlutterAndDivergence"))
writedlm("test/newTestDataGenerators/typicalSectionFlutterAndDivergence/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/typicalSectionFlutterAndDivergence/damps.txt", damps)

# Modal analysis of the axial vibration of a beam under clamped-clamped boundary conditions
include("../plotGenerators/beamAxialVibrationCCPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamAxialVibrationCC"))
writedlm("test/newTestDataGenerators/beamAxialVibrationCC/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamAxialVibrationCC/u1_modeShapes.txt", u1_modeShapes)

# Modal analysis of the axial vibration of a beam under clamped-free boundary conditions
include("../plotGenerators/beamAxialVibrationCFPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamAxialVibrationCF"))
writedlm("test/newTestDataGenerators/beamAxialVibrationCF/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamAxialVibrationCF/u1_modeShapes.txt", u1_modeShapes)

# Modal analysis of the axial vibration of a beam under free-free boundary conditions
include("../plotGenerators/beamAxialVibrationFFPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamAxialVibrationFF"))
writedlm("test/newTestDataGenerators/beamAxialVibrationFF/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamAxialVibrationFF/u1_modeShapes.txt", u1_modeShapes)

# Modal analysis of the bending vibration of a beam under clamped-clamped boundary conditions
include("../plotGenerators/beamBendingVibrationCCPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationCC"))
writedlm("test/newTestDataGenerators/beamBendingVibrationCC/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamBendingVibrationCC/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the bending vibration of a beam under clamped-free boundary conditions
include("../plotGenerators/beamBendingVibrationCFPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationCF"))
writedlm("test/newTestDataGenerators/beamBendingVibrationCF/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamBendingVibrationCF/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the bending vibration of a beam under clamped-pinned boundary conditions
include("../plotGenerators/beamBendingVibrationCPPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationCP"))
writedlm("test/newTestDataGenerators/beamBendingVibrationCP/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamBendingVibrationCP/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the bending vibration of a beam under clamped-sliding boundary conditions
include("../plotGenerators/beamBendingVibrationCSPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationCS"))
writedlm("test/newTestDataGenerators/beamBendingVibrationCS/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamBendingVibrationCS/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the bending vibration of a beam under free-free boundary conditions
include("../plotGenerators/beamBendingVibrationFFPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationFF"))
writedlm("test/newTestDataGenerators/beamBendingVibrationFF/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamBendingVibrationFF/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the bending vibration of a beam under pinned-pinned boundary conditions
include("../plotGenerators/beamBendingVibrationPPPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationPP"))
writedlm("test/newTestDataGenerators/beamBendingVibrationPP/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamBendingVibrationPP/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the torsional vibration of a beam under clamped-clamped boundary conditions
include("../plotGenerators/beamTorsionalVibrationCCPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamTorsionalVibrationCC"))
writedlm("test/newTestDataGenerators/beamTorsionalVibrationCC/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamTorsionalVibrationCC/p1_modeShapes.txt", p1_modeShapes)

# Modal analysis of the torsional vibration of a beam under clamped-free boundary conditions
include("../plotGenerators/beamTorsionalVibrationCFPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamTorsionalVibrationCF"))
writedlm("test/newTestDataGenerators/beamTorsionalVibrationCF/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamTorsionalVibrationCF/p1_modeShapes.txt", p1_modeShapes)

# Modal analysis of the torsional vibration of a beam under free-free boundary conditions
include("../plotGenerators/beamTorsionalVibrationFFPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamTorsionalVibrationFF"))
writedlm("test/newTestDataGenerators/beamTorsionalVibrationFF/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/beamTorsionalVibrationFF/p1_modeShapes.txt", p1_modeShapes)

# Modal analysis of a cantilever beam with a tip axial inertia
include("../plotGenerators/cantileverWithTipAxialMassEigenPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipAxialMassEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipAxialMassEigen/freqsNorm.txt", freqsNorm)
writedlm("test/newTestDataGenerators/cantileverWithTipAxialMassEigen/u1_modeShapes.txt", u1_modeShapes)

# Modal analysis of a cantilever beam with a tip axial spring
include("../plotGenerators/cantileverWithTipAxialSpringEigenPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipAxialSpringEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipAxialSpringEigen/freqsNorm.txt", freqsNorm)
writedlm("test/newTestDataGenerators/cantileverWithTipAxialSpringEigen/u1_modeShapes.txt", u1_modeShapes)

# Modal analysis of a cantilever beam with a tip spring in bending
include("../plotGenerators/cantileverWithTipSpringEigenPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipSpringEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipSpringEigen/freqsNorm.txt", freqsNorm)
writedlm("test/newTestDataGenerators/cantileverWithTipSpringEigen/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of a cantilever beam with a tip torsional inertia
include("../plotGenerators/cantileverWithTipTorsionalInertiaEigenPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipTorsionalInertiaEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipTorsionalInertiaEigen/freqsNorm.txt", freqsNorm)
writedlm("test/newTestDataGenerators/cantileverWithTipTorsionalInertiaEigen/p1_modeShapes.txt", p1_modeShapes)

# Modal analysis of a cantilever beam with a tip torsional spring
include("../plotGenerators/cantileverWithTipTorsionalSpringEigenPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipTorsionalSpringEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipTorsionalSpringEigen/freqsNorm.txt", freqsNorm)
writedlm("test/newTestDataGenerators/cantileverWithTipTorsionalSpringEigen/p1_modeShapes.txt", p1_modeShapes)

# Modal analysis a beam clamped at one end, simply-supported at the other and with a tip inertia
include("../plotGenerators/clampedSSBeamWIthTipInertiaPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedSSBeamWIthTipInertia"))
writedlm("test/newTestDataGenerators/clampedSSBeamWIthTipInertia/freqs.txt", freqs)
writedlm("test/newTestDataGenerators/clampedSSBeamWIthTipInertia/u3_modeShapes.txt", u3_modeShapes)

# Modal analysis of the baseline Healy free FFWT wing without gravity
include("../plotGenerators/HealyBaselineFFWTModalFreePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTModalFree"))
writedlm("test/newTestDataGenerators/HealyBaselineFFWTModalFree/freqs.txt", freqs)

# Modal analyses of the Pazy wing in horizontal and vertical positions
include("../plotGenerators/PazyWingModalPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingModal"))
writedlm("test/newTestDataGenerators/PazyWingModal/freqs.txt", freqs)

# Modal analysis of a beam pinned at one end and transversely springed at the other
include("../plotGenerators/pinnedSpringedBeamEigenPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pinnedSpringedBeamEigen"))
writedlm("test/newTestDataGenerators/pinnedSpringedBeamEigen/freqs.txt", freqs)

# Modal analysis of the sixteen-meter-wing
include("../plotGenerators/SMWModalPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWModal"))
writedlm("test/newTestDataGenerators/SMWModal/freqs.txt", freqs)

# Modal analysis of a straight rotor under varying angular velocities
include("../plotGenerators/straightRotorPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/straightRotor"))
writedlm("test/newTestDataGenerators/straightRotor/freqs.txt", numFreqs)

# Modal analysis of the undeformed swept Pazy wing
include("../plotGenerators/sweptPazyModalPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptPazyModal"))
writedlm("test/newTestDataGenerators/sweptPazyModal/freqs.txt", freqs)

# Modal analysis of a swept-tip rotor under varying angular velocities
include("../plotGenerators/sweptTipRotorPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptTipRotor"))
writedlm("test/newTestDataGenerators/sweptTipRotor/freqs.txt", numFreqs[end,end])

# Modal analysis of a cantilevered tapered beam
include("../plotGenerators/taperedBeamEigenPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/taperedBeamEigen"))
writedlm("test/newTestDataGenerators/taperedBeamEigen/freqs.txt", freqs)

# Modal analysis of the Tang&Dowell wing at varying pitch angles
include("../plotGenerators/TDWingPitchRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/TDWingPitchRange"))
writedlm("test/newTestDataGenerators/TDWingPitchRange/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/TDWingPitchRange/tip_u2.txt", tip_u2)
writedlm("test/newTestDataGenerators/TDWingPitchRange/tip_twist.txt", tip_twist)
writedlm("test/newTestDataGenerators/TDWingPitchRange/freqs.txt", freqs)

# Modal analysis of a 2-story frame
include("../plotGenerators/twoStoryFramePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/twoStoryFrame"))
writedlm("test/newTestDataGenerators/twoStoryFrame/freqs.txt", freqs)

# Static analysis of an arch under a dead pressure load
include("../plotGenerators/archUnderDeadPressurePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/archUnderDeadPressure"))
writedlm("test/newTestDataGenerators/archUnderDeadPressure/mid_u3.txt", mid_u3)

# Static analysis of an arch under a follower pressure load
include("../plotGenerators/archUnderFollowerPressurePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/archUnderFollowerPressure"))
writedlm("test/newTestDataGenerators/archUnderFollowerPressure/mid_u3.txt", mid_u3)

# Static analysis of a cantilever beam with an axial spring attached between its middle and tip nodes, subjected to an axial tip force
include("../plotGenerators/axialDoublyAttachedSpringCantileverPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/axialDoublyAttachedSpringCantilever"))
writedlm("test/newTestDataGenerators/axialDoublyAttachedSpringCantilever/u1.txt", u1)
writedlm("test/newTestDataGenerators/axialDoublyAttachedSpringCantilever/F1.txt", F1)

# Static analysis of a cantilever beam with an axial spring attached between its middle and tip nodes, subjected to an axial tip displacement
include("../plotGenerators/axialDoublyAttachedSpringCantilever2PlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/axialDoublyAttachedSpringCantilever2"))
writedlm("test/newTestDataGenerators/axialDoublyAttachedSpringCantilever2/u1.txt", u1)
writedlm("test/newTestDataGenerators/axialDoublyAttachedSpringCantilever2/F1.txt", F1)

# Static analysis of a biclamped, hinged beam under distributed and concentrated loads
include("../plotGenerators/biclampedHingedBeamPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/biclampedHingedBeam"))
writedlm("test/newTestDataGenerators/biclampedHingedBeam/u3Mid.txt", u3Mid)
writedlm("test/newTestDataGenerators/biclampedHingedBeam/p2Left.txt", p2Left)
writedlm("test/newTestDataGenerators/biclampedHingedBeam/p2Right.txt", p2Right)

# Static analysis of a cantilever beam bending under self weight
include("../plotGenerators/cantileverUnderSelfWeightPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverUnderSelfWeight"))
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/u1.txt", u1)
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/u3.txt", u3)
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/F3.txt", F3)
writedlm("test/newTestDataGenerators/cantileverUnderSelfWeight/M2.txt", M2)

# Static analysis of a cantilever beam with a tip spring in bending
include("../plotGenerators/cantileverWithTipSpringPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipSpring"))
writedlm("test/newTestDataGenerators/cantileverWithTipSpring/u3.txt", u3)
writedlm("test/newTestDataGenerators/cantileverWithTipSpring/F3.txt", F3)
writedlm("test/newTestDataGenerators/cantileverWithTipSpring/M2.txt", M2)

# Static analysis of a clamped beam, with a free flared hinge at the middle, under distributed and concentrated loads
include("../plotGenerators/clampedFlaredHingedBeamPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedFlaredHingedBeam"))
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/u1.txt", u1)
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/u2.txt", u2)
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/u3.txt", u3)
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/p1.txt", p1)
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/p2.txt", p2)
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/p3.txt", p3)
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/F3.txt", F3)
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/M2.txt", M2)
writedlm("test/newTestDataGenerators/clampedFlaredHingedBeam/phiHinge.txt", ϕHinge)

# Static analysis of a clamped beam, hinged at the middle with imposed hinge angle, under distributed and concentrated loads
include("../plotGenerators/clampedHingedBeamPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedHingedBeam"))
writedlm("test/newTestDataGenerators/clampedHingedBeam/u1.txt", u1)
writedlm("test/newTestDataGenerators/clampedHingedBeam/u3.txt", u3)
writedlm("test/newTestDataGenerators/clampedHingedBeam/p2.txt", p2)
writedlm("test/newTestDataGenerators/clampedHingedBeam/F3.txt", F3)
writedlm("test/newTestDataGenerators/clampedHingedBeam/M2.txt", M2)

# Static analysis of a clamped beam rotated in 3D space, hinged at the middle with imposed hinge angle, under distributed and concentrated loads
include("../plotGenerators/clampedHingedBeamRotatedPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedHingedBeamRotated"))
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/u1.txt", u1)
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/u2.txt", u2)
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/u3.txt", u3)
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/p2_b.txt", p2_b)
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/F3.txt", F3)
writedlm("test/newTestDataGenerators/clampedHingedBeamRotated/M2.txt", M2)

# Static analysis of a clamped beam, hinged at the middle, under distributed loads, with varying spring stiffnesses around the hinge
include("../plotGenerators/clampedHingedBeamSpringRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedHingedBeamSpringRange"))
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/u1.txt", u1)
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/u3.txt", u3)
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/p2.txt", p2)
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/F3.txt", F3)
writedlm("test/newTestDataGenerators/clampedHingedBeamSpringRange/M2.txt", M2)

# Static analysis of a clamped beam with a coasting flared hinge at the middle, under distributed loads
include("../plotGenerators/coastingFoldingWingtipPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/coastingFoldingWingtip"))
writedlm("test/newTestDataGenerators/coastingFoldingWingtip/pHinge.txt", pHinge)

# Static analysis of a clamped beam with a driven flared hinge at the middle, under distributed loads
include("../plotGenerators/drivenFoldingWingtipPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/drivenFoldingWingtip"))
writedlm("test/newTestDataGenerators/drivenFoldingWingtip/pHinge.txt", pHinge)

# Static analysis of a clamped beam with a driven, springed, flared hinge at the middle
include("../plotGenerators/drivenSpringedHingedBeamPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/drivenSpringedHingedBeam"))
writedlm("test/newTestDataGenerators/drivenSpringedHingedBeam/pHinge.txt", pHinge)

# Static analysis of composite laminates subjected to tip loads
include("../plotGenerators/compositeCantileverMDPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/compositeCantileverMD"))
for i in 1:3
    writedlm(string("test/newTestDataGenerators/compositeCantileverMD/u1_500mm_b",i,".txt"), hcat(u1_500mm[i,:]...)')
    writedlm(string("test/newTestDataGenerators/compositeCantileverMD/u2_500mm_b",i,".txt"), hcat(u2_500mm[i,:]...)')
    writedlm(string("test/newTestDataGenerators/compositeCantileverMD/u3_500mm_b",i,".txt"), hcat(u3_500mm[i,:]...)')
end

# Static analysis of a curved cantilever subjected to a tip dead force
include("../plotGenerators/curvedCantileverDeadLoadPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/curvedCantileverDeadLoad"))
writedlm("test/newTestDataGenerators/curvedCantileverDeadLoad/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/curvedCantileverDeadLoad/tip_u2.txt", tip_u2)
writedlm("test/newTestDataGenerators/curvedCantileverDeadLoad/tip_u3.txt", tip_u3)

# Static analysis of a curved cantilever subjected to a tip follower force
include("../plotGenerators/curvedCantileverStaticFollowerPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/curvedCantileverStaticFollower"))
writedlm("test/newTestDataGenerators/curvedCantileverStaticFollower/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/curvedCantileverStaticFollower/tip_u2.txt", tip_u2)
writedlm("test/newTestDataGenerators/curvedCantileverStaticFollower/tip_u3.txt", tip_u3)

# Static analysis of a cantilever with distributed follower force
include("../plotGenerators/distributedLoadCantileverPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/distributedLoadCantilever"))
writedlm("test/newTestDataGenerators/distributedLoadCantilever/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/distributedLoadCantilever/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/distributedLoadCantilever/tip_angle.txt", tip_angle)

# Static analysis of a hinged beam subjected to a distributed load
include("../plotGenerators/hingedBeamPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/hingedBeam"))
writedlm("test/newTestDataGenerators/hingedBeam/u3.txt", u3)
writedlm("test/newTestDataGenerators/hingedBeam/p2.txt", p2)
writedlm("test/newTestDataGenerators/hingedBeam/F3.txt", F3)
writedlm("test/newTestDataGenerators/hingedBeam/M2.txt", M2)

# Static analysis of a hinged beam subjected to a distributed load and a rotational spring
include("../plotGenerators/hingedSpringedBeamPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/hingedSpringedBeam"))
writedlm("test/newTestDataGenerators/hingedSpringedBeam/u3.txt", u3)
writedlm("test/newTestDataGenerators/hingedSpringedBeam/p2.txt", p2)
writedlm("test/newTestDataGenerators/hingedSpringedBeam/F3.txt", F3)
writedlm("test/newTestDataGenerators/hingedSpringedBeam/M2.txt", M2)

# Static analysis of a T-frame hinged at the connection subjected to a distributed load
include("../plotGenerators/hingedTFramePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/hingedTFrame"))
writedlm("test/newTestDataGenerators/hingedTFrame/u1_beam1.txt", u1_beam1)
writedlm("test/newTestDataGenerators/hingedTFrame/u1_beam2.txt", u1_beam2)
writedlm("test/newTestDataGenerators/hingedTFrame/u3_beam1.txt", u3_beam1)
writedlm("test/newTestDataGenerators/hingedTFrame/u3_beam2.txt", u3_beam2)
writedlm("test/newTestDataGenerators/hingedTFrame/F3_beam1.txt", F3_beam1)
writedlm("test/newTestDataGenerators/hingedTFrame/M2_beam1.txt", M2_beam1)

# Static analysis of the Lee frame (a right-angled frame) with a dead load
include("../plotGenerators/LeeFrameDeadLoadPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/LeeFrameDeadLoad"))
writedlm("test/newTestDataGenerators/LeeFrameDeadLoad/u1_atForce.txt", u1_atForce)
writedlm("test/newTestDataGenerators/LeeFrameDeadLoad/u3_atForce.txt", u3_atForce)

# Static analysis of the Lee frame (a right-angled frame) with a follower load
include("../plotGenerators/LeeFrameFollowerLoadPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/LeeFrameFollowerLoad"))
writedlm("test/newTestDataGenerators/LeeFrameFollowerLoad/u1_atForce.txt", u1_atForce)
writedlm("test/newTestDataGenerators/LeeFrameFollowerLoad/u3_atForce.txt", u3_atForce)

# Static analysis of the pure bending test of the Pazy wing
include("../plotGenerators/PazyWingBendingTestPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingBendingTest"))
writedlm("test/newTestDataGenerators/PazyWingBendingTest/tip_OOP.txt", tip_OOP)

# Static analysis of the coupled torsion-bending test of the Pazy wing
include("../plotGenerators/PazyWingTorsionTestPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingTorsionTest"))
writedlm("test/newTestDataGenerators/PazyWingTorsionTest/tip_OOP.txt", tip_OOP)
writedlm("test/newTestDataGenerators/PazyWingTorsionTest/tip_twist.txt", tip_twist)

# Static analysis of a mid-loaded arch pinned at one end and clamped at the other
include("../plotGenerators/pinnedClampedArchPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pinnedClampedArch"))
writedlm("test/newTestDataGenerators/pinnedClampedArch/u1.txt", u1_atForce)
writedlm("test/newTestDataGenerators/pinnedClampedArch/u3.txt", u3_atForce)

# Static analysis of a right-angled frame under a tip transverse follower force
include("../plotGenerators/rightAngledFramePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rightAngledFrame"))
writedlm("test/newTestDataGenerators/rightAngledFrame/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/rightAngledFrame/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/rightAngledFrame/tip_angle.txt", tip_angle)

# Static analysis of a right-angled frame under a 'buckling' load
include("../plotGenerators/rightAngledFrameBucklingPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rightAngledFrameBuckling"))
writedlm("test/newTestDataGenerators/rightAngledFrameBuckling/tip_u2.txt", tip_u2)

# Static analysis of an L-frame with a doubly-attached spring and tip load
include("../plotGenerators/springedLFramePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/springedLFrame"))
writedlm("test/newTestDataGenerators/springedLFrame/u3_b.txt", u3_b)

# Static analysis of a swept-back clamped beam with a driven flared folding wingtip
include("../plotGenerators/sweptBackFFWTWingPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptBackFFWTWing"))
writedlm("test/newTestDataGenerators/sweptBackFFWTWing/phi.txt", ϕ)
writedlm("test/newTestDataGenerators/sweptBackFFWTWing/hingeBalanceM.txt", hingeBalanceM)
writedlm("test/newTestDataGenerators/sweptBackFFWTWing/pHinge.txt", pHinge)

# Static analysis of a semi-circular arch with tangential follower force
include("../plotGenerators/tangentiallyForcedArchPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tangentiallyForcedArch"))
writedlm("test/newTestDataGenerators/tangentiallyForcedArch/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/tangentiallyForcedArch/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/tangentiallyForcedArch/tip_angle.txt", tip_angle)

# Static analysis of a cantilever with tip follower transverse force 
include("../plotGenerators/tipFollowerForceCantileverPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipFollowerForceCantilever"))
writedlm("test/newTestDataGenerators/tipFollowerForceCantilever/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/tipFollowerForceCantilever/tip_u3.txt", tip_u3)

# Static analysis of a cantilever with tip follower transverse force (force split over 2 BCs)
include("../plotGenerators/tipFollowerForceCantilever2PlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipFollowerForceCantilever2"))
writedlm("test/newTestDataGenerators/tipFollowerForceCantilever2/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/tipFollowerForceCantilever2/tip_u3.txt", tip_u3)

# Static analysis of a cantilever with tip moment
include("../plotGenerators/tipMomentCantileverPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipMomentCantilever"))
writedlm("test/newTestDataGenerators/tipMomentCantilever/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/tipMomentCantilever/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/tipMomentCantilever/tip_angle.txt", tip_angle)

# Static analysis of a semi-circular arch with transverse follower force
include("../plotGenerators/transverselyForcedArchPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/transverselyForcedArch"))
writedlm("test/newTestDataGenerators/transverselyForcedArch/tip_u1.txt", tip_u1)
writedlm("test/newTestDataGenerators/transverselyForcedArch/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/transverselyForcedArch/tip_angle.txt", tip_angle)

# Static analysis of a beam with triangular distributed load
include("../plotGenerators/triangleLoadBeamPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/triangleLoadBeam"))
writedlm("test/newTestDataGenerators/triangleLoadBeam/u3.txt", u3)
writedlm("test/newTestDataGenerators/triangleLoadBeam/F3.txt", F3)
writedlm("test/newTestDataGenerators/triangleLoadBeam/M2.txt", M2)

# Static analysis of a cantilever beam with a torsional spring attached between its middle and tip nodes, subjected to a torsional tip moment
include("../plotGenerators/twistDoublyAttachedSpringCantileverPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/twistDoublyAttachedSpringCantilever"))
writedlm("test/newTestDataGenerators/twistDoublyAttachedSpringCantilever/p1.txt", p1)
writedlm("test/newTestDataGenerators/twistDoublyAttachedSpringCantilever/M1.txt", M1)

# Steady aeroelastic analysis of the clamped conventional HALE
include("../plotGenerators/conventionalHALEclampedSteadyPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALEclampedSteady"))
writedlm("test/newTestDataGenerators/conventionalHALEclampedSteady/x1_def.txt", x1_def)
writedlm("test/newTestDataGenerators/conventionalHALEclampedSteady/x3_def.txt", x3_def)

# Steady analysis of the baseline Healy free FFWT wing with varying root pitch and airspeed
include("../plotGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoastPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoast"))
writedlm("test/newTestDataGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoast/phiHinge.txt", ϕHinge)
writedlm("test/newTestDataGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoast/u3Hinge.txt", u3Hinge)
writedlm("test/newTestDataGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoast/M2root.txt", M2root)

# Steady analysis of the Healy sideslip FFWT wing with varying flare angle and root pitch angle
include("../plotGenerators/HealySideslipFFWTsteadyFlareRangeAoARangeCoastPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealySideslipFFWTsteadyFlareRangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/phiHinge.txt", ϕHinge)

# Steady analysis of the Healy sideslip FFWT wing with varying flare angle, root pitch angle and airspeed
include("../plotGenerators/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoastPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/phiHinge.txt", ϕHinge)

# Steady analysis of the Healy sideslip FFWT wing with varying wingtip twist, root pitch angle and sideslip angle
include("../plotGenerators/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoastPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast"))
writedlm("test/newTestDataGenerators/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast/phiHinge.txt", ϕHinge)

# Steady analysis of the Pazy wing with a coasting FFWT
include("../plotGenerators/PazyFFWTsteadyCoastPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyCoast"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyCoast/phiHinge.txt", ϕHinge)

# Steady analysis of the Pazy wing with a FFWT at a fixed fold angle
include("../plotGenerators/PazyFFWTsteadyFixedFoldPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyFixedFold"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyFixedFold/phiHinge.txt", ϕHinge)
writedlm("test/newTestDataGenerators/PazyFFWTsteadyFixedFold/hingeBalanceM.txt", hingeBalanceM)

# Steady analysis of the Pazy wing with a coasting FFWT, at varying airspeed and root pitch angle
include("../plotGenerators/PazyFFWTsteadyURangeAoARangeCoastPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyURangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURangeAoARangeCoast/phiHinge.txt", ϕHinge)

# Steady analysis of the Pazy wing with a coasting FFWT, at varying airspeed
include("../plotGenerators/PazyFFWTsteadyURangeCoastPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyURangeCoast"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURangeCoast/phiHinge.txt", ϕHinge)

# Steady analysis of the Pazy wing with a FFWT at a fixed fold angle, at varying airspeed
include("../plotGenerators/PazyFFWTsteadyURangeFixedFoldPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyURangeFixedFold"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURangeFixedFold/phiHinge.txt", ϕHinge)
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURangeFixedFold/hingeMoment.txt", hingeMoment)

# Steady analysis of the Pazy wing with varying root pitch angle
include("../plotGenerators/PazyWingPitchRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingPitchRange"))
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_AoA.txt", tip_AoA)
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_OOP.txt", tip_OOP)
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_IP.txt", tip_IP)
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_twist.txt", tip_twist)

# Steady aeroelastic analysis of the sixteen-meter-wing
include("../plotGenerators/SMWSteadyPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWSteady"))
writedlm("test/newTestDataGenerators/SMWSteady/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/SMWSteady/tip_twist.txt", tip_twist)

# Steady aeroelastic analysis of the Tang&Dowell wing at varying airspeed
include("../plotGenerators/TDWingAirspeedRangePlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/TDWingAirspeedRange"))
writedlm("test/newTestDataGenerators/TDWingAirspeedRange/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/TDWingAirspeedRange/tip_twist.txt", tip_twist)

# Trim analysis of the Blended-Wing-Body flying wing in free flight
include("../plotGenerators/BWBtrimPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/BWBtrim"))
writedlm("test/newTestDataGenerators/BWBtrim/trimAoA.txt", trimAoA)
writedlm("test/newTestDataGenerators/BWBtrim/trimThrust.txt", trimThrust)
writedlm("test/newTestDataGenerators/BWBtrim/trimDelta.txt", trimδ)

# Trim analysis of the conventional HALE aircraft in free flight (considering aerodynamics from stabilizers and thrust)
include("../plotGenerators/conventionalHALEfullTrimPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALEfullTrim"))
writedlm("test/newTestDataGenerators/conventionalHALEfullTrim/trimAoA.txt", trimAoA)
writedlm("test/newTestDataGenerators/conventionalHALEfullTrim/trimThrust.txt", trimThrust)
writedlm("test/newTestDataGenerators/conventionalHALEfullTrim/trimDelta.txt", trimδ)

# Trim analysis of the conventional HALE aircraft in free flight at rigid and flexible configurations (neglecting aerodynamics from stabilizers and thrust)
include("../plotGenerators/conventionalHALEtrimPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALEtrim"))
writedlm("test/newTestDataGenerators/conventionalHALEtrim/trimAoA.txt", trimAoA)

# Trim analysis of the Helios flying-wing
include("../plotGenerators/heliosTrimPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosTrim"))
writedlm("test/newTestDataGenerators/heliosTrim/trimAoA.txt", trimAoA)
writedlm("test/newTestDataGenerators/heliosTrim/trimThrust.txt", trimThrust)
writedlm("test/newTestDataGenerators/heliosTrim/trimDelta.txt", trimδ)

# Trim analysis (reaction loads check) of a beam loaded at the middle
include("../plotGenerators/freeBeamTrimPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/freeBeamTrim"))
writedlm("test/newTestDataGenerators/freeBeamTrim/F3.txt", F3)
writedlm("test/newTestDataGenerators/freeBeamTrim/M2.txt", M2)

# Trim analysis (reaction loads check) of a simply-supported beam loaded at the middle
include("../plotGenerators/midLoadedBeamTrimPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/midLoadedBeamTrim"))
writedlm("test/newTestDataGenerators/midLoadedBeamTrim/F3.txt", F3)
writedlm("test/newTestDataGenerators/midLoadedBeamTrim/M2.txt", M2)

# Trim analysis (reaction loads check) of a right-angled frame
include("../plotGenerators/rightAngledFrameTrimPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rightAngledFrameTrim"))
writedlm("test/newTestDataGenerators/rightAngledFrameTrim/balanceHorizontalForce.txt", balanceHorizontalForce)
writedlm("test/newTestDataGenerators/rightAngledFrameTrim/balanceVerticalForce.txt", balanceVerticalForce)

# Trim analysis of a cantilever with tip force
include("../plotGenerators/tipLoadedCantileverTrimPlotGenerator.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipLoadedCantileverTrim"))
writedlm("test/newTestDataGenerators/tipLoadedCantileverTrim/u3.txt", u3)
writedlm("test/newTestDataGenerators/tipLoadedCantileverTrim/F3.txt", F3)
writedlm("test/newTestDataGenerators/tipLoadedCantileverTrim/M2.txt", M2)

println("Finished generating new test data and plots")
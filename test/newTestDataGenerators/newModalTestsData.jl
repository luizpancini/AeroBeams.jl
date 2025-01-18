# New reference data for modal tests

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

# Modal analysis of the baseline Healy free FFWT wing without gravity
include("../examples/HealyBaselineFFWTModalFree.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTModalFree"))
writedlm("test/newTestDataGenerators/HealyBaselineFFWTModalFree/freqs.txt", freqs)

# # Modal analyses of the Pazy wing in horizontal and vertical positions
# include("../examples/PazyWingModal.jl")
# mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingModal"))
# writedlm("test/newTestDataGenerators/PazyWingModal/freqs.txt", freqs)

# # Modal analysis of a beam pinned at one end and transversely springed at the other
# include("../examples/pinnedSpringedBeamEigen.jl")
# mkpath(string(pwd(),"/test/newTestDataGenerators/pinnedSpringedBeamEigen"))
# writedlm("test/newTestDataGenerators/pinnedSpringedBeamEigen/freqs.txt", freqs)

# # Modal analysis of the sixteen-meter-wing
# include("../examples/SMWModal.jl")
# mkpath(string(pwd(),"/test/newTestDataGenerators/SMWModal"))
# writedlm("test/newTestDataGenerators/SMWModal/freqs.txt", freqs)

# # Modal analysis of a straight rotor under varying angular velocities
# include("../examples/straightRotor.jl")
# mkpath(string(pwd(),"/test/newTestDataGenerators/straightRotor"))
# writedlm("test/newTestDataGenerators/straightRotor/freqs.txt", numFreqs)

# # Modal analysis of a swept-tip rotor under varying angular velocities
# include("../examples/sweptTipRotor.jl")
# mkpath(string(pwd(),"/test/newTestDataGenerators/sweptTipRotor"))
# writedlm("test/newTestDataGenerators/sweptTipRotor/freqs.txt", numFreqs[end,end])

# # Modal analysis of a cantilevered tapered beam
# include("../examples/taperedBeamEigen.jl")
# mkpath(string(pwd(),"/test/newTestDataGenerators/taperedBeamEigen"))
# writedlm("test/newTestDataGenerators/taperedBeamEigen/freqs.txt", freqs)

# # Modal analysis of the Tang&Dowell wing at varying pitch angles
# include("../examples/TDWingPitchRange.jl")
# mkpath(string(pwd(),"/test/newTestDataGenerators/TDWingPitchRange"))
# writedlm("test/newTestDataGenerators/TDWingPitchRange/tip_u3.txt", tip_u3)
# writedlm("test/newTestDataGenerators/TDWingPitchRange/tip_u2.txt", tip_u2)
# writedlm("test/newTestDataGenerators/TDWingPitchRange/tip_twist.txt", tip_twist)
# writedlm("test/newTestDataGenerators/TDWingPitchRange/freqs.txt", freqs)

# # Modal analysis of a 2-story frame
# include("../examples/twoStoryFrame.jl")
# mkpath(string(pwd(),"/test/newTestDataGenerators/twoStoryFrame"))
# writedlm("test/newTestDataGenerators/twoStoryFrame/freqs.txt", freqs)
# Modal analysis of the axial vibration of a beam under clamped-clamped boundary conditions
include("../examples/eigenStructural/beamAxialVibrationCC.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamAxialVibrationCC"))
writedlm("test/newTestDataGenerators/beamAxialVibrationCC/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/beamAxialVibrationCC/u1_modeShapes.txt", test_data_matrix(u1_modeShapes))

# Modal analysis of the axial vibration of a beam under clamped-free boundary conditions
include("../examples/eigenStructural/beamAxialVibrationCF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamAxialVibrationCF"))
writedlm("test/newTestDataGenerators/beamAxialVibrationCF/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/beamAxialVibrationCF/u1_modeShapes.txt", test_data_matrix(u1_modeShapes))

# Modal analysis of the axial vibration of a beam under free-free boundary conditions
include("../examples/eigenStructural/beamAxialVibrationFF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamAxialVibrationFF"))
writedlm("test/newTestDataGenerators/beamAxialVibrationFF/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/beamAxialVibrationFF/u1_modeShapes.txt", test_data_matrix(u1_modeShapes))

# Modal analysis of the bending vibration of a beam under clamped-clamped boundary conditions
include("../examples/eigenStructural/beamBendingVibrationCC.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationCC"))
writedlm("test/newTestDataGenerators/beamBendingVibrationCC/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/beamBendingVibrationCC/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the bending vibration of a beam under clamped-free boundary conditions
include("../examples/eigenStructural/beamBendingVibrationCF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationCF"))
writedlm("test/newTestDataGenerators/beamBendingVibrationCF/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/beamBendingVibrationCF/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the bending vibration of a beam under clamped-pinned boundary conditions
include("../examples/eigenStructural/beamBendingVibrationCP.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationCP"))
writedlm("test/newTestDataGenerators/beamBendingVibrationCP/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/beamBendingVibrationCP/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the bending vibration of a beam under clamped-sliding boundary conditions
include("../examples/eigenStructural/beamBendingVibrationCS.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationCS"))
writedlm("test/newTestDataGenerators/beamBendingVibrationCS/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/beamBendingVibrationCS/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the bending vibration of a beam under free-free boundary conditions
include("../examples/eigenStructural/beamBendingVibrationFF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationFF"))
writedlm("test/newTestDataGenerators/beamBendingVibrationFF/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/beamBendingVibrationFF/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the bending vibration of a beam under pinned-pinned boundary conditions
include("../examples/eigenStructural/beamBendingVibrationPP.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamBendingVibrationPP"))
writedlm("test/newTestDataGenerators/beamBendingVibrationPP/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/beamBendingVibrationPP/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the torsional vibration of a beam under clamped-clamped boundary conditions
include("../examples/eigenStructural/beamTorsionalVibrationCC.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamTorsionalVibrationCC"))
writedlm("test/newTestDataGenerators/beamTorsionalVibrationCC/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/beamTorsionalVibrationCC/p1_modeShapes.txt", test_data_matrix(p1_modeShapes))

# Modal analysis of the torsional vibration of a beam under clamped-free boundary conditions
include("../examples/eigenStructural/beamTorsionalVibrationCF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamTorsionalVibrationCF"))
writedlm("test/newTestDataGenerators/beamTorsionalVibrationCF/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/beamTorsionalVibrationCF/p1_modeShapes.txt", test_data_matrix(p1_modeShapes))

# Modal analysis of the torsional vibration of a beam under free-free boundary conditions
include("../examples/eigenStructural/beamTorsionalVibrationFF.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/beamTorsionalVibrationFF"))
writedlm("test/newTestDataGenerators/beamTorsionalVibrationFF/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/beamTorsionalVibrationFF/p1_modeShapes.txt", test_data_matrix(p1_modeShapes))

# Modal analysis of a cantilever beam with a tip axial inertia
include("../examples/eigenStructural/cantileverWithTipAxialMassEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipAxialMassEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipAxialMassEigen/freqsNorm.txt", test_data_matrix(freqsNorm))
writedlm("test/newTestDataGenerators/cantileverWithTipAxialMassEigen/u1_modeShapes.txt", test_data_matrix(u1_modeShapes))

# Modal analysis of a cantilever beam with a tip axial spring
include("../examples/eigenStructural/cantileverWithTipAxialSpringEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipAxialSpringEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipAxialSpringEigen/freqsNorm.txt", test_data_matrix(freqsNorm))
writedlm("test/newTestDataGenerators/cantileverWithTipAxialSpringEigen/u1_modeShapes.txt", test_data_matrix(u1_modeShapes))

# Modal analysis of a cantilever beam with a tip spring in bending
include("../examples/eigenStructural/cantileverWithTipSpringEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipSpringEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipSpringEigen/freqsNorm.txt", test_data_matrix(freqsNorm))
writedlm("test/newTestDataGenerators/cantileverWithTipSpringEigen/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of a cantilever beam with a tip torsional inertia
include("../examples/eigenStructural/cantileverWithTipTorsionalInertiaEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipTorsionalInertiaEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipTorsionalInertiaEigen/freqsNorm.txt", test_data_matrix(freqsNorm))
writedlm("test/newTestDataGenerators/cantileverWithTipTorsionalInertiaEigen/p1_modeShapes.txt", test_data_matrix(p1_modeShapes))

# Modal analysis of a cantilever beam with a tip torsional spring
include("../examples/eigenStructural/cantileverWithTipTorsionalSpringEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cantileverWithTipTorsionalSpringEigen"))
writedlm("test/newTestDataGenerators/cantileverWithTipTorsionalSpringEigen/freqsNorm.txt", test_data_matrix(freqsNorm))
writedlm("test/newTestDataGenerators/cantileverWithTipTorsionalSpringEigen/p1_modeShapes.txt", test_data_matrix(p1_modeShapes))

# Modal analysis a beam clamped at one end, simply-supported at the other and with a tip inertia
include("../examples/eigenStructural/clampedSSBeamWIthTipInertiaEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/clampedSSBeamWIthTipInertiaEigen"))
writedlm("test/newTestDataGenerators/clampedSSBeamWIthTipInertiaEigen/freqs.txt", test_data_matrix(freqs))
writedlm("test/newTestDataGenerators/clampedSSBeamWIthTipInertiaEigen/u3_modeShapes.txt", test_data_matrix(u3_modeShapes))

# Modal analysis of the baseline Healy free FFWT wing without gravity
include("../examples/eigenStructural/HealyBaselineFFWTModalFree.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTModalFree"))
writedlm("test/newTestDataGenerators/HealyBaselineFFWTModalFree/freqs.txt", test_data_matrix(freqs))

# Modal analyses of the Pazy wing in horizontal and vertical positions
include("../examples/eigenStructural/PazyWingModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingModal"))
writedlm("test/newTestDataGenerators/PazyWingModal/freqs.txt", test_data_matrix(freqs))

# Modal analysis of a beam pinned at one end and transversely springed at the other
include("../examples/eigenStructural/pinnedSpringedBeamEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pinnedSpringedBeamEigen"))
writedlm("test/newTestDataGenerators/pinnedSpringedBeamEigen/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the sixteen-meter-wing
include("../examples/eigenStructural/SMWModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWModal"))
writedlm("test/newTestDataGenerators/SMWModal/freqs.txt", test_data_matrix(freqs))

# Modal analysis of a straight rotor under varying angular velocities
include("../examples/eigenStructural/straightRotor.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/straightRotor"))
writedlm("test/newTestDataGenerators/straightRotor/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the undeformed swept Pazy wing
include("../examples/eigenStructural/sweptPazyModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptPazyModal"))
writedlm("test/newTestDataGenerators/sweptPazyModal/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the baseline Healy locked FFWT wing
include("../examples/eigenStructural/HealyLCOFFWTLockedModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyLCOFFWTLockedModal"))
writedlm("test/newTestDataGenerators/HealyLCOFFWTLockedModal/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the baseline Healy FFWT wing with wingtip removed
include("../examples/eigenStructural/HealyLCOFFWTRemovedModal.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyLCOFFWTRemovedModal"))
writedlm("test/newTestDataGenerators/HealyLCOFFWTRemovedModal/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the swept Pazy wing with varying tip displacement
include("../examples/eigenStructural/sweptPazyModalTipDispRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptPazyModalTipDispRange"))
writedlm("test/newTestDataGenerators/sweptPazyModalTipDispRange/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/sweptPazyModalTipDispRange/tip_twist.txt", test_data_matrix(tip_twist))
writedlm("test/newTestDataGenerators/sweptPazyModalTipDispRange/freqs.txt", test_data_matrix(freqsOsc))

# Modal analysis of a swept-tip rotor under varying angular velocities
include("../examples/eigenStructural/sweptTipRotor.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptTipRotor"))
writedlm("test/newTestDataGenerators/sweptTipRotor/freqs.txt", test_data_matrix(freqs))

# Modal analysis of a cantilevered tapered beam
include("../examples/eigenStructural/taperedBeamEigen.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/taperedBeamEigen"))
writedlm("test/newTestDataGenerators/taperedBeamEigen/freqs.txt", test_data_matrix(freqs))

# Modal analysis of the Tang&Dowell wing at varying pitch angles
include("../examples/eigenStructural/TDWingPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/TDWingPitchRange"))
writedlm("test/newTestDataGenerators/TDWingPitchRange/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/TDWingPitchRange/tip_u2.txt", test_data_matrix(tip_u2))
writedlm("test/newTestDataGenerators/TDWingPitchRange/tip_twist.txt", test_data_matrix(tip_twist))
writedlm("test/newTestDataGenerators/TDWingPitchRange/freqs.txt", test_data_matrix(freqs))

# Modal analysis of a 2-story frame
include("../examples/eigenStructural/twoStoryFrame.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/twoStoryFrame"))
writedlm("test/newTestDataGenerators/twoStoryFrame/freqs.txt", test_data_matrix(freqs))
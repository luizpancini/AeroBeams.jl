# Dynamic analysis of the Blended-Wing-Body vehicle undergoing a checked pitch maneuver
include("../examples/dynamicAeroelastic/BWBcheckedPitchManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/BWBcheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/BWBcheckedPitchManeuver/rootAoA.txt", test_data_matrix(rootAoA))
writedlm("test/newTestDataGenerators/BWBcheckedPitchManeuver/Deltau3.txt", test_data_matrix(Δu3))

# Dynamic analysis of the conventional HALE aircraft undergoing a checked pitch maneuver
include("../examples/dynamicAeroelastic/cHALECheckedPitchManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cHALECheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/cHALECheckedPitchManeuver/wingAoA.txt", test_data_matrix(wingAoA))
writedlm("test/newTestDataGenerators/cHALECheckedPitchManeuver/Deltau3.txt", test_data_matrix(Δu3))

# Dynamic analysis of the conventional HALE aircraft undergoing a coordinated turn maneuver
include("../examples/dynamicAeroelastic/cHALECheckedRollManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cHALECheckedRollManeuver"))
writedlm("test/newTestDataGenerators/cHALECheckedRollManeuver/wingAoA.txt", test_data_matrix(wingAoA))
writedlm("test/newTestDataGenerators/cHALECheckedRollManeuver/Deltau3.txt", test_data_matrix(Δu3))

# Dynamic analysis of the baseline Healy FFWT wing under a series of one-minus-cosine gusts
include("../examples/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTOMCGustFloating"))
writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTOMCGustFloating/M2root.txt"),test_data_matrix(M2root))
writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTOMCGustFloating/fold.txt"),test_data_matrix(ϕ))

# Dynamic analysis of the baseline Healy FFWT wing for a specific one-minus-cosine gust
include("../examples/dynamicAeroelastic/HealyBaselineFFWTOMCGustFloating2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTOMCGustFloating2"))
writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTOMCGustFloating2/M2root.txt"),test_data_matrix(M2root))
writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTOMCGustFloating2/fold.txt"),test_data_matrix(ϕ))

# Dynamic analysis of the Helios flying-wing undergoing a checked pitch maneuver
include("../examples/dynamicAeroelastic/heliosCheckedPitchManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosCheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/heliosCheckedPitchManeuver/rootAoA.txt", test_data_matrix(rootAoA))
writedlm("test/newTestDataGenerators/heliosCheckedPitchManeuver/Deltau3.txt", test_data_matrix(Δu3))

# Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over time
include("../examples/dynamicAeroelastic/PazyWingContinuous1DGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingContinuous1DGust"))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tipAoA.txt", test_data_matrix(tipAoA))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tqSpan_cm.txt", test_data_matrix(tqSpan_cm))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DGust/tqSpan_ct.txt", test_data_matrix(tqSpan_ct))

# Dynamic analysis of the Pazy wing encountering a continuous, 1-dimensional gust defined over space
include("../examples/dynamicAeroelastic/PazyWingContinuous1DSpaceGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingContinuous1DSpaceGust"))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tipAoA.txt", test_data_matrix(tipAoA))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tqSpan_cm.txt", test_data_matrix(tqSpan_cm))
writedlm("test/newTestDataGenerators/PazyWingContinuous1DSpaceGust/tqSpan_ct.txt", test_data_matrix(tqSpan_ct))

# Dynamic analysis of the Pazy wing encountering a continuous, 2-dimensional gust
include("../examples/dynamicAeroelastic/PazyWingContinuous2DSpaceGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingContinuous2DSpaceGust"))
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tipAoA.txt", test_data_matrix(tipAoA))
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tqSpan_cm.txt", test_data_matrix(tqSpan_cm))
writedlm("test/newTestDataGenerators/PazyWingContinuous2DSpaceGust/tqSpan_ct.txt", test_data_matrix(tqSpan_ct))

# Dynamic analysis of the Pazy wing encountering a DARPA gust
include("../examples/dynamicAeroelastic/PazyWingDARPAGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingDARPAGust"))
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tipAoA.txt", test_data_matrix(tipAoA))
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tqSpan_cm.txt", test_data_matrix(tqSpan_cm))
writedlm("test/newTestDataGenerators/PazyWingDARPAGust/tqSpan_ct.txt", test_data_matrix(tqSpan_ct))

# Dynamic analysis of the Pazy wing encountering a one-minus-cosine gust
include("../examples/dynamicAeroelastic/PazyWingOMCGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingOMCGust"))
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tipAoA.txt", test_data_matrix(tipAoA))
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tqSpan_cm.txt", test_data_matrix(tqSpan_cm))
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tqSpan_ct.txt", test_data_matrix(tqSpan_ct))

# Dynamic analysis of the Pazy wing with a tip impulse force
include("../examples/dynamicAeroelastic/PazyWingTipImpulse.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingTipImpulse"))
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tipAoA.txt", test_data_matrix(tipAoA))
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tqSpan_cm.txt", test_data_matrix(tqSpan_cm))
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tqSpan_ct.txt", test_data_matrix(tqSpan_ct))

# Dynamic analysis of the baseline Healy locked FFWT wing under a series of one-minus-cosine gusts
include("../examples/dynamicAeroelastic/HealyBaselineFFWTlockedOMCGustFloating.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTlockedOMCGustFloating"))
writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTlockedOMCGustFloating/M2root.txt"),test_data_matrix(M2root))

# Dynamic analysis of the baseline Healy free and locked FFWT wings under one-minus-cosine gusts
include("../examples/dynamicAeroelastic/HealyBaselineFFWTlockedVfreeOMCGustFloating.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTlockedVfreeOMCGustFloating"))
writedlm(string("test/newTestDataGenerators/HealyBaselineFFWTlockedVfreeOMCGustFloating/dM2peaks.txt"),test_data_matrix(ΔM2peaks))

# Dynamic analysis of the Healy sideslip FFWT wing under a one-minus-cosine gust
include("../examples/dynamicAeroelastic/HealySideslipFFWTOMCGustFloating.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealySideslipFFWTOMCGustFloating"))
writedlm("test/newTestDataGenerators/HealySideslipFFWTOMCGustFloating/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/HealySideslipFFWTOMCGustFloating/M2_root.txt", test_data_matrix(M2_root))
writedlm("test/newTestDataGenerators/HealySideslipFFWTOMCGustFloating/root_alphae.txt", test_data_matrix(root_αₑ))
writedlm("test/newTestDataGenerators/HealySideslipFFWTOMCGustFloating/root_cn.txt", test_data_matrix(root_cn))

# Dynamic analysis of the Pazy FFWT wing under a one-minus-cosine gust
include("../examples/dynamicAeroelastic/PazyFFWTOMCGustFloating.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTOMCGustFloating"))
writedlm("test/newTestDataGenerators/PazyFFWTOMCGustFloating/fold.txt", test_data_matrix(ϕ))
writedlm("test/newTestDataGenerators/PazyFFWTOMCGustFloating/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/PazyFFWTOMCGustFloating/M2_root.txt", test_data_matrix(M2_root))
writedlm("test/newTestDataGenerators/PazyFFWTOMCGustFloating/tqSpan_cn.txt", test_data_matrix(tqSpan_cn))

# Dynamic analysis of a rotating wing using the dynamic stall model
include("../examples/dynamicAeroelastic/rotorDStest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rotorDStest"))
writedlm("test/newTestDataGenerators/rotorDStest/cl.txt", test_data_matrix(cl))
writedlm("test/newTestDataGenerators/rotorDStest/cm.txt", test_data_matrix(cm))
writedlm("test/newTestDataGenerators/rotorDStest/cdrag.txt", test_data_matrix(cdrag))

# Dynamic analysis of a pitching wing with dynamic stall models
include("../examples/dynamicAeroelastic/wingDStest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/wingDStest"))
writedlm("test/newTestDataGenerators/wingDStest/cl.txt", test_data_matrix(cl))
writedlm("test/newTestDataGenerators/wingDStest/cm.txt", test_data_matrix(cm))
writedlm("test/newTestDataGenerators/wingDStest/cdrag.txt", test_data_matrix(cdrag))
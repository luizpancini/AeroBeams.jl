# New reference data for dynamic aeroelastic tests

# Dynamic analysis of the Blended-Wing-Body vehicle undergoing a checked pitch maneuver
include("../examples/BWBcheckedPitchManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/BWBcheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/BWBcheckedPitchManeuver/rootAoA.txt", rootAoA)
writedlm("test/newTestDataGenerators/BWBcheckedPitchManeuver/Deltau3.txt", Δu3)

# Dynamic analysis of the conventional HALE aircraft undergoing a checked pitch maneuver
include("../examples/conventionalHALECheckedPitchManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALECheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/conventionalHALECheckedPitchManeuver/wingAoA.txt", wingAoA)
writedlm("test/newTestDataGenerators/conventionalHALECheckedPitchManeuver/Deltau3.txt", Δu3)

# Dynamic analysis of the conventional HALE aircraft undergoing a coordinated turn maneuver
include("../examples/conventionalHALECheckedRollManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALECheckedRollManeuver"))
writedlm("test/newTestDataGenerators/conventionalHALECheckedRollManeuver/wingAoA.txt", wingAoA)
writedlm("test/newTestDataGenerators/conventionalHALECheckedRollManeuver/Deltau3.txt", Δu3)

# Dynamic analysis of the Helios flying-wing undergoing a checked pitch maneuver
include("../examples/heliosCheckedPitchManeuver.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosCheckedPitchManeuver"))
writedlm("test/newTestDataGenerators/heliosCheckedPitchManeuver/rootAoA.txt", rootAoA)
writedlm("test/newTestDataGenerators/heliosCheckedPitchManeuver/Deltau3.txt", Δu3)

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

# Dynamic analysis of the Pazy wing encountering a one-minus-cosine gust
include("../examples/PazyWingOMCGust.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingOMCGust"))
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingOMCGust/tqSpan_ct.txt", tqSpan_ct)

# Dynamic analysis of the Pazy wing with a tip impulse force
include("../examples/PazyWingTipImpulse.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingTipImpulse"))
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tipAoA.txt", tipAoA)
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tipOOP.txt", tipOOP)
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tqSpan_cn.txt", tqSpan_cn)
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tqSpan_cm.txt", tqSpan_cm)
writedlm("test/newTestDataGenerators/PazyWingTipImpulse/tqSpan_ct.txt", tqSpan_ct)
# New reference data for steady aeroelastic tests

# Steady aeroelastic analysis of the clamped conventional HALE
include("../examples/conventionalHALEclampedSteady.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALEclampedSteady"))
writedlm("test/newTestDataGenerators/conventionalHALEclampedSteady/x1_def.txt", x1_def)
writedlm("test/newTestDataGenerators/conventionalHALEclampedSteady/x3_def.txt", x3_def)

# Steady analysis of the baseline Healy free FFWT wing with varying root pitch and airspeed
include("../examples/HealyBaselineFFWTsteadyAoARangeURangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoast"))
writedlm("test/newTestDataGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoast/phiHinge.txt", ϕHinge)
writedlm("test/newTestDataGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoast/u3Hinge.txt", u3Hinge)
writedlm("test/newTestDataGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoast/M2root.txt", M2root)

# Steady analysis of the Healy FFWT wing with varying flare angle and root pitch angle
include("../examples/HealyFFWTsteadyFlareRangeAoARangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyFFWTsteadyFlareRangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/HealyFFWTsteadyFlareRangeAoARangeCoast/phiHinge.txt", ϕHinge)

# Steady analysis of the Healy FFWT wing with varying flare angle, root pitch angle and airspeed
include("../examples/HealyFFWTsteadyFlareRangeURangeAoARangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyFFWTsteadyFlareRangeURangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/HealyFFWTsteadyFlareRangeURangeAoARangeCoast/phiHinge.txt", ϕHinge)

# Steady analysis of the Healy FFWT wing with varying wingtip twist, root pitch angle and sideslip angle
include("../examples/HealyFFWTsteadyTwistRangeAoARangeSideslipRangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyFFWTsteadyTwistRangeAoARangeSideslipRangeCoast"))
writedlm("test/newTestDataGenerators/HealyFFWTsteadyTwistRangeAoARangeSideslipRangeCoast/phiHinge.txt", ϕHinge)

# Steady analysis of the Pazy wing with a coasting FFWT
include("../examples/PazyFFWTsteadyCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyCoast"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyCoast/phiHinge.txt", ϕHinge)

# Steady analysis of the Pazy wing with a FFWT at a fixed fold angle
include("../examples/PazyFFWTsteadyFixedFold.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyFixedFold"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyFixedFold/phiHinge.txt", ϕHinge)
writedlm("test/newTestDataGenerators/PazyFFWTsteadyFixedFold/hingeBalanceM.txt", hingeBalanceM)

# Steady analysis of the Pazy wing with a coasting FFWT, at varying airspeed and root pitch angle
include("../examples/PazyFFWTsteadyURangeAoARangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyURangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURangeAoARangeCoast/phiHinge.txt", ϕHinge)

# Steady analysis of the Pazy wing with a coasting FFWT, at varying airspeed
include("../examples/PazyFFWTsteadyURangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyURangeCoast"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURangeCoast/phiHinge.txt", ϕHinge)

# Steady analysis of the Pazy wing with a FFWT at a fixed fold angle, at varying airspeed
include("../examples/PazyFFWTsteadyURangeFixedFold.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyURangeFixedFold"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURangeFixedFold/phiHinge.txt", ϕHinge)
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURangeFixedFold/hingeMoment.txt", hingeMoment)

# Steady analysis of the Pazy wing with varying root pitch angle
include("../examples/PazyWingPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingPitchRange"))
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_AoA.txt", tip_AoA)
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_OOP.txt", tip_OOP)
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_IP.txt", tip_IP)
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_twist.txt", tip_twist)

# Steady aeroelastic analysis of the sixteen-meter-wing
include("../examples/SMWSteady.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWSteady"))
writedlm("test/newTestDataGenerators/SMWSteady/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/SMWSteady/tip_twist.txt", tip_twist)

# Steady aeroelastic analysis of the Tang&Dowell wing at varying airspeed
include("../examples/TDWingAirspeedRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/TDWingAirspeedRange"))
writedlm("test/newTestDataGenerators/TDWingAirspeedRange/tip_u3.txt", tip_u3)
writedlm("test/newTestDataGenerators/TDWingAirspeedRange/tip_twist.txt", tip_twist)
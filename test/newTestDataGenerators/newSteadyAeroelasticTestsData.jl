# New reference data for steady aeroelastic tests

# Steady aeroelastic analysis of the clamped conventional HALE
include("../examples/conventionalHALEclampedSteady.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALEclampedSteady"))
writedlm("test/newTestDataGenerators/conventionalHALEclampedSteady/x1_def.txt", x1_def)
writedlm("test/newTestDataGenerators/conventionalHALEclampedSteady/x3_def.txt", x3_def)

# Steady analysis of the Pazy wing with flared folding wing tip (FFWT)
include("../examples/PazyFFWTsteady.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteady"))
writedlm("test/newTestDataGenerators/PazyFFWTsteady/u1.txt", u1)
writedlm("test/newTestDataGenerators/PazyFFWTsteady/u3.txt", u3)
writedlm("test/newTestDataGenerators/PazyFFWTsteady/p1.txt", p1)
writedlm("test/newTestDataGenerators/PazyFFWTsteady/p2.txt", p2)
writedlm("test/newTestDataGenerators/PazyFFWTsteady/M2.txt", M2)
writedlm("test/newTestDataGenerators/PazyFFWTsteady/alpha.txt", α)
writedlm("test/newTestDataGenerators/PazyFFWTsteady/cn.txt", cn)
writedlm("test/newTestDataGenerators/PazyFFWTsteady/hingeBalanceM.txt", hingeBalanceM)

# Steady analysis of the Pazy wing with flared folding wing tip (FFWT) and varying airspeed
include("../examples/PazyFFWTsteadyURange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyURange"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURange/u1.txt", u1)
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURange/u3.txt", u3)
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURange/p1.txt", p1)
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURange/p2.txt", p2)
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURange/M2.txt", M2)
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURange/alpha.txt", α)
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURange/cn.txt", cn)

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
# New reference data for steady aeroelastic tests

# Steady analysis of the Pazy wing with flared folding wing tip (FFWT)
include("../examples/PazyFFWTsteady.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteady"))
writedlm("test/newTestDataGenerators/PazyFFWTsteady/u3_of_x1.txt", u3_of_x1)
writedlm("test/newTestDataGenerators/PazyFFWTsteady/p2_of_x1.txt", p2_of_x1)
writedlm("test/newTestDataGenerators/PazyFFWTsteady/M2_of_x1.txt", M2_of_x1)

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
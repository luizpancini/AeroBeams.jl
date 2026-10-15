# Steady aeroelastic analysis of the clamped conventional HALE
include("../examples/steadyAeroelastic/cHALEclampedSteady.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cHALEclampedSteady"))
writedlm("test/newTestDataGenerators/cHALEclampedSteady/x1_def.txt", test_data_matrix(x1_def))
writedlm("test/newTestDataGenerators/cHALEclampedSteady/x3_def.txt", test_data_matrix(x3_def))

# Steady analysis of the baseline Healy free FFWT wing with varying root pitch and airspeed
include("../examples/steadyAeroelastic/HealyBaselineFFWTsteadyAoARangeURangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoast"))
writedlm("test/newTestDataGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))
writedlm("test/newTestDataGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoast/u3Hinge.txt", test_data_matrix(u3Hinge))
writedlm("test/newTestDataGenerators/HealyBaselineFFWTsteadyAoARangeURangeCoast/M2root.txt", test_data_matrix(M2root))

# Steady analysis of the Healy sideslip FFWT wing with varying flare angle and root pitch angle
include("../examples/steadyAeroelastic/HealySideslipFFWTsteadyFlareRangeAoARangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealySideslipFFWTsteadyFlareRangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/HealySideslipFFWTsteadyFlareRangeAoARangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the Healy sideslip FFWT wing with varying flare angle, root pitch angle and airspeed
include("../examples/steadyAeroelastic/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/HealySideslipFFWTsteadyFlareRangeURangeAoARangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the Healy sideslip FFWT wing with varying wingtip twist, root pitch angle and sideslip angle
include("../examples/steadyAeroelastic/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast"))
writedlm("test/newTestDataGenerators/HealySideslipFFWTsteadyTwistRangeAoARangeSideslipRangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the Pazy wing with a coasting FFWT
include("../examples/steadyAeroelastic/PazyFFWTsteadyCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyCoast"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the Pazy wing with a FFWT at a fixed fold angle
include("../examples/steadyAeroelastic/PazyFFWTsteadyFixedFold.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyFixedFold"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyFixedFold/phiHinge.txt", test_data_matrix(ϕHinge))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyFixedFold/hingeMoment.txt", test_data_matrix(hingeMoment))

# Steady analysis of the Pazy wing with a coasting FFWT, at varying airspeed and root pitch angle
include("../examples/steadyAeroelastic/PazyFFWTsteadyURangeAoARangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyURangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURangeAoARangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the Pazy wing with a coasting FFWT, at varying airspeed
include("../examples/steadyAeroelastic/PazyFFWTsteadyURangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyURangeCoast"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Steady analysis of the Pazy wing with a FFWT at a fixed fold angle, at varying airspeed
include("../examples/steadyAeroelastic/PazyFFWTsteadyURangeFixedFold.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyFFWTsteadyURangeFixedFold"))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURangeFixedFold/phiHinge.txt", test_data_matrix(ϕHinge))
writedlm("test/newTestDataGenerators/PazyFFWTsteadyURangeFixedFold/hingeMoment.txt", test_data_matrix(hingeMoment))

# Steady analysis of the Pazy wing with varying root pitch angle
include("../examples/steadyAeroelastic/PazyWingPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/PazyWingPitchRange"))
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_AoA.txt", test_data_matrix(tip_AoA))
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_OOP.txt", test_data_matrix(tip_OOP))
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_IP.txt", test_data_matrix(tip_IP))
writedlm("test/newTestDataGenerators/PazyWingPitchRange/tip_twist.txt", test_data_matrix(tip_twist))

# Steady aeroelastic analysis of the sixteen-meter-wing
include("../examples/steadyAeroelastic/SMWSteady.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SMWSteady"))
writedlm("test/newTestDataGenerators/SMWSteady/tip_u3.txt", test_data_matrix(tip_u3))
writedlm("test/newTestDataGenerators/SMWSteady/tip_twist.txt", test_data_matrix(tip_twist))

# Steady analysis of the Healy LCO FFWT wing with varying airspeed and root pitch angle
include("../examples/steadyAeroelastic/HealyLCOFFWTsteadyURangeAoARangeCoast.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/HealyLCOFFWTsteadyURangeAoARangeCoast"))
writedlm("test/newTestDataGenerators/HealyLCOFFWTsteadyURangeAoARangeCoast/phiHinge.txt", test_data_matrix(ϕHinge))

# Creating the conventional HALE model
include("../examples/steadyAeroelastic/cHALEmodel.jl")

# Steady analysis of the swept Pazy wing with varying root pitch angle
include("../examples/steadyAeroelastic/sweptPazySteadyPitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptPazySteadyPitchRange"))
writedlm("test/newTestDataGenerators/sweptPazySteadyPitchRange/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/sweptPazySteadyPitchRange/tipTwist.txt", test_data_matrix(tipTwist))
writedlm("test/newTestDataGenerators/sweptPazySteadyPitchRange/tipAoA.txt", test_data_matrix(tipAoA))

# Steady analysis of the swept Pazy wing with varying sweep and root pitch angles
include("../examples/steadyAeroelastic/sweptPazySteadySweepRangePitchRange.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptPazySteadySweepRangePitchRange"))
writedlm("test/newTestDataGenerators/sweptPazySteadySweepRangePitchRange/tipOOP.txt", test_data_matrix(tipOOP))
writedlm("test/newTestDataGenerators/sweptPazySteadySweepRangePitchRange/tipTwist.txt", test_data_matrix(tipTwist))

# Steady lift distribution of swept wings
include("../examples/steadyAeroelastic/sweptWingsLiftDistribution.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/sweptWingsLiftDistribution"))
writedlm("test/newTestDataGenerators/sweptWingsLiftDistribution/cn.txt", test_data_matrix(cn))
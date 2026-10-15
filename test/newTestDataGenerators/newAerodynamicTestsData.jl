# Dynamic analysis of an airfoils with harmonic flap deflection profile
include("../examples/aerodynamic/flapOscillation.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flapOscillation"))
writedlm("test/newTestDataGenerators/flapOscillation/cn.txt", test_data_matrix(cn))
writedlm("test/newTestDataGenerators/flapOscillation/cm.txt", test_data_matrix(cm))

# Dynamic analysis of two airfoils with linked harmonic flap deflection profiles
include("../examples/aerodynamic/flapOscillationLinked.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flapOscillationLinked"))
writedlm("test/newTestDataGenerators/flapOscillationLinked/cnMaster.txt", test_data_matrix(cnMaster))
writedlm("test/newTestDataGenerators/flapOscillationLinked/cmMaster.txt", test_data_matrix(cmMaster))
writedlm("test/newTestDataGenerators/flapOscillationLinked/cnSlave.txt", test_data_matrix(cnSlave))
writedlm("test/newTestDataGenerators/flapOscillationLinked/cmSlave.txt", test_data_matrix(cmSlave))

# One-minus-cosine gust response of an airfoil section at several pitch angles
include("../examples/aerodynamic/OMCgustTests.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/OMCgustTests"))
writedlm("test/newTestDataGenerators/OMCgustTests/dcl.txt", test_data_matrix(Δcl))

# One-minus-cosine gust response of an airfoil section at several pitch angles
include("../examples/aerodynamic/OMCgustTests2.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/OMCgustTests2"))
writedlm("test/newTestDataGenerators/OMCgustTests2/dcl.txt", test_data_matrix(Δcl))

# Dynamic analysis of a harmonically pitching airfoil, using the dynamic stall model
include("../examples/aerodynamic/pitchingAirfoilDSModelTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pitchingAirfoilDSModelTest"))
writedlm("test/newTestDataGenerators/pitchingAirfoilDSModelTest/"*frameString*"_cn.txt", test_data_matrix(cn))
writedlm("test/newTestDataGenerators/pitchingAirfoilDSModelTest/"*frameString*"_cm.txt", test_data_matrix(cm))
writedlm("test/newTestDataGenerators/pitchingAirfoilDSModelTest/"*frameString*"_ct.txt", test_data_matrix(ct))

# Dynamic analysis of a harmonically plunging airfoil, using the dynamic stall model
include("../examples/aerodynamic/plungingAirfoilDSModelTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/plungingAirfoilDSModelTest"))
writedlm("test/newTestDataGenerators/plungingAirfoilDSModelTest/cn.txt", test_data_matrix(cn))
writedlm("test/newTestDataGenerators/plungingAirfoilDSModelTest/cm.txt", test_data_matrix(cm))

# Sharp-edged gust response of an airfoil section at several pitch angles
include("../examples/aerodynamic/SEgustTests.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SEgustTests"))
writedlm("test/newTestDataGenerators/SEgustTests/dcl.txt", test_data_matrix(Δcl))

# Dynamic analysis of an airfoil section sinusoidally surging (facing a time-varying freestream)
include("../examples/aerodynamic/timeVaryingFreestream.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/timeVaryingFreestream"))
writedlm("test/newTestDataGenerators/timeVaryingFreestream/cn.txt", test_data_matrix(cn[end]))
writedlm("test/newTestDataGenerators/timeVaryingFreestream/cm.txt", test_data_matrix(cm[end]))

# Dynamic analysis of an airfoil section sinusoidally pitching and surging (facing a time-varying freestream)
include("../examples/aerodynamic/timeVaryingFreestreamAndPitch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/timeVaryingFreestreamAndPitch"))
writedlm("test/newTestDataGenerators/timeVaryingFreestreamAndPitch/cn.txt", test_data_matrix(cn[end]))
writedlm("test/newTestDataGenerators/timeVaryingFreestreamAndPitch/cm.txt", test_data_matrix(cm[end]))
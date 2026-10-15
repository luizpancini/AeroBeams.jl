# Trim analysis of the Blended-Wing-Body flying wing in free flight
include("../examples/trimAeroelastic/BWBtrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/BWBtrim"))
writedlm("test/newTestDataGenerators/BWBtrim/trimAoA.txt", test_data_matrix(trimAoA))
writedlm("test/newTestDataGenerators/BWBtrim/trimThrust.txt", test_data_matrix(trimThrust))
writedlm("test/newTestDataGenerators/BWBtrim/trimDelta.txt", test_data_matrix(trimδ))

# Trim analysis of the conventional HALE aircraft in free flight (considering aerodynamics from stabilizers and thrust)
include("../examples/trimAeroelastic/cHALEfullTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cHALEfullTrim"))
writedlm("test/newTestDataGenerators/cHALEfullTrim/trimAoA.txt", test_data_matrix(trimAoA))
writedlm("test/newTestDataGenerators/cHALEfullTrim/trimThrust.txt", test_data_matrix(trimThrust))
writedlm("test/newTestDataGenerators/cHALEfullTrim/trimDelta.txt", test_data_matrix(trimδ))

# Trim analysis of the conventional HALE aircraft in free flight at rigid and flexible configurations (neglecting aerodynamics from stabilizers and thrust)
include("../examples/trimAeroelastic/cHALEtrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/cHALEtrim"))
writedlm("test/newTestDataGenerators/cHALEtrim/trimAoA.txt", test_data_matrix(trimAoA))

# Trim analysis of the Helios flying-wing
include("../examples/trimAeroelastic/heliosTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosTrim"))
writedlm("test/newTestDataGenerators/heliosTrim/trimAoA.txt", test_data_matrix(trimAoA))
writedlm("test/newTestDataGenerators/heliosTrim/trimThrust.txt", test_data_matrix(trimThrust))
writedlm("test/newTestDataGenerators/heliosTrim/trimDelta.txt", test_data_matrix(trimδ))
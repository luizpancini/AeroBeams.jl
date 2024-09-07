# New reference data for trim aeroelastic tests

# Trim analysis of the Blended-Wing-Body flying wing in free flight
include("../examples/BWBtrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/BWBtrim"))
writedlm("test/newTestDataGenerators/BWBtrim/trimAoA.txt", trimAoA)
writedlm("test/newTestDataGenerators/BWBtrim/trimThrust.txt", trimThrust)
writedlm("test/newTestDataGenerators/BWBtrim/trimDelta.txt", trimδ)

# Trim analysis of the conventional HALE aircraft in free flight (considering aerodynamics from stabilizers and thrust)
include("../examples/conventionalHALEfullTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALEfullTrim"))
writedlm("test/newTestDataGenerators/conventionalHALEfullTrim/trimAoA.txt", trimAoA)
writedlm("test/newTestDataGenerators/conventionalHALEfullTrim/trimThrust.txt", trimThrust)
writedlm("test/newTestDataGenerators/conventionalHALEfullTrim/trimDelta.txt", trimδ)

# Trim analysis of the conventional HALE aircraft in free flight at rigid and flexible configurations (neglecting aerodynamics from stabilizers and thrust)
include("../examples/conventionalHALEtrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/conventionalHALEtrim"))
writedlm("test/newTestDataGenerators/conventionalHALEtrim/trimAoA.txt", trimAoA)

# Trim analysis of the Helios flying-wing
include("../examples/heliosTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/heliosTrim"))
writedlm("test/newTestDataGenerators/heliosTrim/trimAoA.txt", trimAoA)
writedlm("test/newTestDataGenerators/heliosTrim/trimThrust.txt", trimThrust)
writedlm("test/newTestDataGenerators/heliosTrim/trimDelta.txt", trimδ)
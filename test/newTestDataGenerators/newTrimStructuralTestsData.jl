# New reference data for trim structural tests

# Trim analysis (reaction loads check) of a beam loaded at the middle
include("../examples/freeBeamTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/freeBeamTrim"))
writedlm("test/newTestDataGenerators/freeBeamTrim/F3.txt", F3)
writedlm("test/newTestDataGenerators/freeBeamTrim/M2.txt", M2)

# Trim analysis (reaction loads check) of a simply-supported beam loaded at the middle
include("../examples/midLoadedBeamTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/midLoadedBeamTrim"))
writedlm("test/newTestDataGenerators/midLoadedBeamTrim/F3.txt", F3)
writedlm("test/newTestDataGenerators/midLoadedBeamTrim/M2.txt", M2)

# Trim analysis (reaction loads check) of a right-angled frame
include("../examples/rightAngledFrameTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rightAngledFrameTrim"))
writedlm("test/newTestDataGenerators/rightAngledFrameTrim/balanceHorizontalForce.txt", balanceHorizontalForce)
writedlm("test/newTestDataGenerators/rightAngledFrameTrim/balanceVerticalForce.txt", balanceVerticalForce)

# Trim analysis of a cantilever with tip force
include("../examples/tipLoadedCantileverTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipLoadedCantileverTrim"))
writedlm("test/newTestDataGenerators/tipLoadedCantileverTrim/u3.txt", u3)
writedlm("test/newTestDataGenerators/tipLoadedCantileverTrim/F3.txt", F3)
writedlm("test/newTestDataGenerators/tipLoadedCantileverTrim/M2.txt", M2)
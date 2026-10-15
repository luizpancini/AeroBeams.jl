# Trim analysis (reaction loads check) of a beam loaded at the middle
include("../examples/trimStructural/freeBeamTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/freeBeamTrim"))
writedlm("test/newTestDataGenerators/freeBeamTrim/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/freeBeamTrim/M2.txt", test_data_matrix(M2))

# Trim analysis (reaction loads check) of a simply-supported beam loaded at the middle
include("../examples/trimStructural/midDispBeamTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/midDispBeamTrim"))
writedlm("test/newTestDataGenerators/midDispBeamTrim/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/midDispBeamTrim/M2.txt", test_data_matrix(M2))

# Trim analysis (reaction loads check) of a right-angled frame
include("../examples/trimStructural/rightAngledFrameTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/rightAngledFrameTrim"))
writedlm("test/newTestDataGenerators/rightAngledFrameTrim/balanceHorizontalForce.txt", test_data_matrix(balanceHorizontalForce))
writedlm("test/newTestDataGenerators/rightAngledFrameTrim/balanceVerticalForce.txt", test_data_matrix(balanceVerticalForce))

# Trim analysis of a cantilever with tip force
include("../examples/trimStructural/tipLoadedCantileverTrim.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/tipLoadedCantileverTrim"))
writedlm("test/newTestDataGenerators/tipLoadedCantileverTrim/u3.txt", test_data_matrix(u3))
writedlm("test/newTestDataGenerators/tipLoadedCantileverTrim/F3.txt", test_data_matrix(F3))
writedlm("test/newTestDataGenerators/tipLoadedCantileverTrim/M2.txt", test_data_matrix(M2))
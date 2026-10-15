using Test, DelimitedFiles

# Flag to run plot/visualization sections inside example scripts
runExamplePlots = false

# Flag for CI
is_ci = get(ENV, "CI", "false") == "true" || get(ENV, "GITHUB_ACTIONS", "false") == "true"

# Default tolerances for self-comparison
SELFatol = is_ci ? 1e-4 : 1e-6
SELFrtol = is_ci ? 1e-2 : 1e-4

# Load test data helpers
include("testData.jl")

# Tests
include("testDataTests.jl")
include("dynamicAeroelasticTests.jl")
include("eigenAeroelasticTests.jl")
include("aerodynamicTests.jl")
include("trimAeroelasticTests.jl")
include("steadyAeroelasticTests.jl")
include("dynamicStructuralTests.jl")
include("eigenStructuralTests.jl")
include("staticStructuralTests.jl")
include("trimStructuralTests.jl")

println("Finished tests")
using Test, BenchmarkTools, DelimitedFiles

# Option to locally run tests failing in CI
runCIfailingTests = false

# Flag for CI
is_ci = get(ENV, "CI", "false") == "true" || get(ENV, "GITHUB_ACTIONS", "false") == "true"

# Default tolerances for self-comparison
SELFatol = is_ci ? 1e-3 : 1e-4
SELFrtol = is_ci ? 0 : 1e-4

include("dynamicAeroelasticTests.jl")
include("eigenAeroelasticTests.jl")
include("trimAeroelasticTests.jl")
include("steadyAeroelasticTests.jl")
include("aerodynamicTests.jl")
include("dynamicStructuralTests.jl")
include("trimStructuralTests.jl")
include("modalTests.jl")
include("staticStructuralTests.jl")

println("Finished tests")
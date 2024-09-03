using Test, BenchmarkTools, DelimitedFiles

# Default absolute tolerance for self-comparison
SELFatol = 1e-4

include("staticStructuralTests.jl")
include("trimStructuralTests.jl")
include("modalTests.jl")
include("dynamicStructuralTests.jl")
include("aerodynamicTests.jl")
include("steadyAeroelasticTests.jl")
include("trimAeroelasticTests.jl")
include("eigenAeroelasticTests.jl")
include("dynamicAeroelasticTests.jl")

println("Finished tests")
using Test, BenchmarkTools, DelimitedFiles

# Default absolute tolerance for self-comparison
SELFatol = 1e-4

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
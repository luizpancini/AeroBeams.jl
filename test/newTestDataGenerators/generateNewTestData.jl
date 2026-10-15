using DelimitedFiles

include("../testData.jl")

include("newAerodynamicTestsData.jl")
include("newDynamicAeroelasticTestsData.jl")
include("newDynamicStructuralTestsData.jl")
include("newEigenAeroelasticTestsData.jl")
include("newEigenStructuralTestsData.jl")
include("newStaticStructuralTestsData.jl")
include("newSteadyAeroelasticTestsData.jl")
include("newTrimAeroelasticTestsData.jl")
include("newTrimStructuralTestsData.jl")

println("Finished generating new test data")
"""
    Module for the AeroBeams package
"""
module AeroBeams

include("Utilities.jl")
include("UnitsSystem.jl")
include("PointInertia.jl")
include("Atmosphere.jl")
include("Airfoil.jl")
include("AeroSolver.jl")
include("AeroSurface.jl")
include("AeroProperties.jl")
include("Gust.jl")
include("Beam.jl")
include("Element.jl")
include("BC.jl")
include("SpecialNode.jl")
include("TrimLink.jl")
include("Model.jl")
include("Problem.jl")
include("SystemSolver.jl")
include("Aerodynamics.jl")
include("Core.jl")
include("SampleBeams.jl")

end

"""
    Module for the AeroBeams package
"""
module AeroBeams

include("Utilities.jl")
include("PointInertia.jl")
include("Beam.jl")
include("Element.jl")
include("BC.jl")
include("SpecialNode.jl")
include("Model.jl")
include("Problem.jl")
include("SystemSolver.jl")
include("Core.jl")

end

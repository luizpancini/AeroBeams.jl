module AeroBeams

using DSP, FFTW, FiniteDifferences, ForwardDiff, LinearAlgebra, LinearInterpolations, Interpolations, Parameters, QuadGK, Random, Revise, SparseArrays, Statistics, Hungarian

include("UnitsSystem.jl")
include("PointInertia.jl")
include("Spring.jl")
include("Atmosphere.jl")
include("Airfoil.jl")
include("AeroSolver.jl")
include("AeroSurface.jl")
include("AeroProperties.jl")
include("Gust.jl")
include("Utilities.jl")
include("Beam.jl")
include("HingeAxisConstraint.jl")
include("Element.jl")
include("BC.jl")
include("SpecialNode.jl")
include("Links.jl")
include("Model.jl")
include("Problem.jl")
include("SystemSolver.jl")
include("Aerodynamics.jl")
include("Core.jl")
include("SampleModels.jl")
include("Visualization.jl")

end

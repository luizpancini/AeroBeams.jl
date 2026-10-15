module AeroBeams

using ColorSchemes
using DelimitedFiles
using ForwardDiff
using LinearAlgebra
using Parameters
using Plots
using PrettyTables
using Random
using SparseArrays
using StaticArrays

import ADTypes
import ArnoldiMethod
import Arpack
import DSP
import FFTW
import FiniteDifferences
import FunctionWrappers
import Hungarian
import KrylovKit
import LinearMaps
import LineSearches
import LinearInterpolations
import NLsolve
import Plots: annotate!
import QuadGK
import Roots

include("UnitsSystem.jl")
include("PointInertia.jl")
include("Spring.jl")
include("Atmosphere.jl")
include("AeroSolver.jl")
include("Airfoil.jl")
include("AeroSurface.jl")
include("AeroProperties.jl")
include("Gust.jl")
include("Element.jl")
include("Beam.jl")
include("HingeAxisConstraint.jl")
include("BC.jl")
include("SpecialNode.jl")
include("Links.jl")
include("Model.jl")
include("SystemSolver.jl")
include("Problem.jl")
include("Aerodynamics.jl")
include("Core.jl")
include("SampleModels.jl")
include("Visualization.jl")
include("Utilities.jl")

# ColorSchemes exports
export cgrad

# DelimitedFiles exports
export readdlm, writedlm

# LinearAlgebra exports
export cross, Diagonal, diagm, dot, I, norm, tr

# LinearInterpolations exports
export interpolate, Interpolate

# Plots exports
export @animate, @layout, annotate!, font, gif, gr, heatmap, plot, plot!, plot_color, savefig, scatter, scatter!, Shape, surface, surface!, text, xlabel!, ylabel!

# PrettyTables exports
export pretty_table

# Random exports
export rand, randn

# Initialization utility
function __init__()
    Plots.default(dpi=300)
end

end

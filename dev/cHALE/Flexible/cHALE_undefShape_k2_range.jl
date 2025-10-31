using AeroBeams

# Stiffness factor
λ = 1

# Options for stabilizers
stabilizersAero = true
includeVS = true

# Discretization
nElemWing = 20
nElemTailBoom = 5
nElemHorzStabilizer = 4
nElemVertStabilizer = 2

# Gravity
g = 0

# Bending pre-curvature range
k2Range = range(-0.015, 0.045, 5)

# System solver
maxIter = 100
σ0 = 1.0
NR = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=true)

# Initialize outputs
steadyProblem = Array{SteadyProblem}(undef,length(k2Range))

# Sweep bending pre-curvature
for (i,k2) in enumerate(k2Range)
    println("Solving for k2 = $k2")
    # Model for steady problem
    model,_ = create_conventional_HALE(stiffnessFactor=λ,k2=k2,g=g,airspeed=1e-4,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,rootClamped=true)
    # Create and solve steady problem
    steadyProblem[i] = create_SteadyProblem(model=model,systemSolver=NR)
    solve!(steadyProblem[i])
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALE_undefShape_k2_range"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot undeformed shape across k2 range
plt_undefShapes = plot_steady_deformations(steadyProblem,backendSymbol=:gr,plotBCs=false,plotDistLoads=false,plotAxes=false,plotGrid=false,legendEntries=["\$k_2 = $(k2) \$" for k2 in k2Range],legendPos=:top,view=(30,15),plotLimits=([-L*2/3,L*2/3],[-L*2/3,L*2/3],[-L*2/3,L*2/3]),save=true,savePath=string(relPath,"/cHALE_undefShape_k2_range_lambda",λ,".pdf"))
display(plt_undefShapes)

println("Finished cHALE_undefShape_k2_range.jl")
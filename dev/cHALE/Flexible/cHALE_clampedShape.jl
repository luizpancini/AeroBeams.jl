using AeroBeams

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 5

# Altitude
h = 20e3

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Discretization
nElemWing = 40
nElemTailBoom = 5
nElemHorzStabilizer = 4
nElemVertStabilizer = 2

# Airspeed
U = 20

# Bending pre-curvature
k2 = 0.0

# System solver
maxIter = 100
σ0 = 1.0
NR = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)

# Model for clamped trim problem
model,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,k2=k2,hasInducedDrag=hasInducedDrag,rootClamped=true)

# Create and solve steady problem
problem = create_SteadyProblem(model=model,systemSolver=NR)
solve!(problem)

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALE_clampedShape"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot clamped and undeformed shapes
plt_trimmedShape = plot_steady_deformation(problem,backendSymbol=:gr,plotBCs=false,plotDistLoads=false,plotAxes=false,plotGrid=false,legendPos=(0.47,0.7),view=(30,30),plotLimits=([-L*2/3,L*2/3],[-L*2/3,L*2/3],[-L*2/3,L*2/3]),save=true,savePath=string(relPath,"/cHALE_clampedShape_lambda",λ,"_k2",k2,"_U",U,".pdf"))
display(plt_trimmedShape)

println("Finished cHALE_clampedShape.jl")
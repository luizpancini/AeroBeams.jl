using AeroBeams

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 1

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
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,pseudoInverseMethod=:dampedLeastSquares,displayStatus=false)

# Model for clamped trim problem
model,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)

# Create and solve trim problem
trimProblem = create_TrimProblem(model=model,systemSolver=NR)
solve!(trimProblem)

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALE_trimShape"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot trimmed shape and undeformed shape
plt_trimmedShape = plot_steady_deformation(trimProblem,backendSymbol=:gr,plotBCs=false,plotDistLoads=false,plotAxes=false,plotGrid=false,legendPos=(0.47,0.7),view=(30,30),plotLimits=([-L*2/3,L*2/3],[-L*2/3,L*2/3],[-L*2/3,L*2/3]),save=true,savePath=string(relPath,"/cHALE_trimShape_lambda",λ,"_k2",k2,"_U",U,".pdf"))
display(plt_trimmedShape)

println("Finished cHALE_trimShape.jl")
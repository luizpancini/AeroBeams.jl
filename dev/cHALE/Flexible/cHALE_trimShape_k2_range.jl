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
nElemWing = 80
nElemTailBoom = 10
nElemHorzStabilizer = 10
nElemVertStabilizer = 5

# Airspeed
U = 20

# Bending pre-curvature range
k2Range = range(-0.015, 0.045, 5)

# System solver
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,pseudoInverseMethod=:dampedLeastSquares,displayStatus=false)

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(k2Range))

# Sweep bending pre-curvature
for (i,k2) in enumerate(k2Range)
    println("Solving for k2 = $k2, U = $U m/s")
    # Model for trim problem
    model,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)
    # Create and solve trim problem
    trimProblem[i] = create_TrimProblem(model=model,systemSolver=NR)
    solve!(trimProblem[i])
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALE_trimShape_k2_range"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot trimmed shape across k2 range
plt_trimmedShapes = plot_steady_deformations(trimProblem,backendSymbol=:gr,plotBCs=false,plotDistLoads=false,plotAxes=false,plotGrid=false,legendEntries=["\$k_2 = $(k2) \$" for k2 in k2Range],legendPos=:top,view=(30,30),plotLimits=([-L*2/3,L*2/3],[-L*2/3,L*2/3],[-L*2/3,L*2/3]),save=true,savePath=string(relPath,"/cHALE_trimShape_k2_range_lambda",λ,"_U",U,".pdf"))
display(plt_trimmedShapes)

println("Finished cHALE_trimShape_k2_range.jl")
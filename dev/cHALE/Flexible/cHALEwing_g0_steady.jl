using AeroBeams, LinearInterpolations, JLD2

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 1

# Altitude
h = 20e3

# No gravity
g = 0

# Root pitch angle [rad]
θ = 0*π/180

# Option to include induced drag
hasInducedDrag = true

# Parasite drag
cd0 = 1e-2

# Discretization
nElem = 40

# System solver
maxIter = 50
σ0 = 1
NR = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)

# Bending pre-curvature
k2 = 0.045

# Airspeed
U = 1e-4

# Tip in-plane force and root moment (for illustration)
dummyBeam = create_Beam(length=1,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e12)])
tipIPforce = create_BC(name="tipIPforce",beam=dummyBeam,node=nElem+1,types=["F2b"],values=[-1])
rootMoment = create_BC(name="rootMoment",beam=dummyBeam,node=2,types=["M1b"],values=[-1])

# Model
model,L = create_SMW(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElem=nElem,altitude=h,g=g,θ=θ,k2=k2,cd0=cd0,hasInducedDrag=hasInducedDrag,additionalBCs=[tipIPforce,rootMoment])

# Create and solve steady problem
problem = create_SteadyProblem(model=model,systemSolver=NR)
solve!(problem)

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALEwing_g0_steady"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Steady plot
plt_steady = plot_steady_deformation(problem,backendSymbol=:gr,plotUndeformed=false,plotBCs=true,plotDistLoads=false,plotLegend=false,plotAxes=false,plotGrid=false,view=(45,30),loadsSizeScaler=2,plotLimits=([0,L],[-L/2,L/2],[-L/2,L/2]),save=true,savePath=string(relPath,"/cHALEwing_steady.pdf"))
display(plt_steady)

println("Finished cHALEwing_g0_steady.jl")
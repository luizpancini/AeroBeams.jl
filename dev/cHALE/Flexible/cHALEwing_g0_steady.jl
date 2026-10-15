using AeroBeams

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 1

# Altitude
h = 20e3

# No gravity
g = 0

# Root pitch angle [rad]
θ = deg2rad(0)

# Option to include induced drag
hasInducedDrag = true

# Parasite drag
cd0 = 1e-2

# Discretization
nElem = 40

# System solver
systemSolver = create_NLsolver()

# Bending pre-curvature
k2 = 0.045

# Airspeed
U = 1e-4

# Tip in-plane force and root moment (for illustration)
dummyBeam = create_Beam(length=1,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e12)])
tipIPforce = create_BC(name="tipIPforce",beam=dummyBeam,node=nElem+1,conditions=Dict(:F2b=>-1))
rootMoment = create_BC(name="rootMoment",beam=dummyBeam,node=2,conditions=Dict(:M1b=>-1))

# Model
model,L = create_SMW(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElem=nElem,altitude=h,g=g,pitchAngle=θ,k2=k2,cd0=cd0,hasInducedDrag=hasInducedDrag,additionalBCs=[tipIPforce,rootMoment])

# Create and solve steady problem
problem = create_SteadyProblem(model=model,systemSolver=systemSolver)
solve!(problem)

## Plots
# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALEwing_g0_steady"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Steady plot
plt_steady = plot_steady_deformation(problem,plotUndeformed=false,plotBCs=true,plotDistLoads=false,plotLegend=false,plotAxes=false,plotGrid=false,view=(45,30),loadsSizeScaler=2,plotLimits=([0,L],[-L/2,L/2],[-L/2,L/2]),save=true,savePath=string(relPath,"/cHALEwing_steady.pdf"))
display(plt_steady)

println("Finished cHALEwing_g0_steady.jl")
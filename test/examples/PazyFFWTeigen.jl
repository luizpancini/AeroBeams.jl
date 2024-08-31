using AeroBeams

# Hinge node, hinge angle [rad] and flare angle [deg]
hingeNode = 14
hingeAngle = -1*π/2
flareAngle = 10

# Spring stiffness
kSpring = 1e-1

# Root pitch angle
θ = 7*pi/180

# Airspeed
U = 50

# Gravity
g = 9.80665

# Pazy wing with flared folding tip
pazyFFWT,_ = create_PazyFFWT(hingeNode=hingeNode,hingeAngle=hingeAngle,flareAngle=flareAngle,kSpring=kSpring,airspeed=U,p0=[0;0;θ],g=g)

# System solver
σ0 = 1.0
maxIter = 50
NR = create_NewtonRaphson(displayStatus=true,initialLoadFactor=σ0,maximumIterations=maxIter)

# Number of modes
nModes = 6

# Create and solve problem
problem = create_EigenProblem(model=pazyFFWT,systemSolver=NR,nModes=nModes)
solve!(problem)

# Get outputs
freqs = problem.frequenciesOscillatory
damps = round_off!(problem.dampingsOscillatory,1e-8)

println("Finished PazyFFWTeigen.jl")
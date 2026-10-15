using AeroBeams

# Beam
R,θ = 100, π/4
A = 1
Iy = Iz = 1/12
J = Iy+Iz
L = R*θ
E,ν = 1e7, 0
G = E/(2*(1+ν))
nElem = 16
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=E*A,GAy=G*A,GAz=G*A,GJ=G*J,EIy=E*Iy,EIz=E*Iz)],k=[0;0;1/R])

# BCs
F = 600
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,conditions=Dict(:F3A=>F))

# Model
curvedCantilever = create_Model(name="curvedCantilever",beams=[beam],BCs=[clamp,tipForce])

# Set system solver options
σ0 = 0.05
σstep = 0.05
NR = create_NewtonRaphson(initLoadFac=σ0,maxLoadFacStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=curvedCantilever,systemSolver=NR)
solve!(problem)

# Get tip displacements
tipNodeFinalDisp = problem.loadstepHistory.states[end][nElem].nodal.u2
tip_u1 = tipNodeFinalDisp[1] 
tip_u2 = tipNodeFinalDisp[2]
tip_u3 = tipNodeFinalDisp[3]

# Reference solution by Bathe & Bolourchi - Large displacement analysis of three-dimensional beam structures (1979)
u1Ref = -23.5
u2Ref = -13.4
u3Ref = 53.4

# Compute and print relative differences
ϵu1 = tip_u1/u1Ref - 1
ϵu2 = tip_u2/u2Ref - 1
ϵu3 = tip_u3/u3Ref - 1
println("Relative differences: u1 = $(round(100*ϵu1,digits=4)) %, u2 = $(round(100*ϵu2,digits=4)) %, u3 = $(round(100*ϵu3,digits=4)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/curvedCantileverDeadLoad"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/curvedCantileverDeadLoad_deformation.pdf"))
display(deformationPlot)

println("Finished curvedCantileverDeadLoad.jl")
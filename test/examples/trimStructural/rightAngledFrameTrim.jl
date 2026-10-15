using AeroBeams

# Beam
L = 1
∞ = 1e6
nElem1 = 2
nElem2 = 1
beam1 = create_Beam(name="beam1",length=L,nElements=nElem1,S=[isotropic_stiffness_matrix(∞=∞)])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem2,S=[isotropic_stiffness_matrix(∞=∞)],p0=transform_to_WM([0;π/2;0],from=:E321))

# BCs - balanceLoads are the equivalent reactions at the pin
F = 1
roller = create_BC(name="roller",beam=beam1,node=1,conditions=Dict(:u3A=>0))
pin = create_BC(name="pin",beam=beam1,node=nElem1+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
verticalForce = create_BC(name="verticalForce",beam=beam1,node=div(nElem1,2)+1,conditions=Dict(:F3A=>F))
horizontalForce = create_BC(name="horizontalForce",beam=beam2,node=nElem2+1,conditions=Dict(:F1A=>F))
balanceLoads = create_BC(name="balanceLoads",beam=beam2,node=1,conditions=Dict(:F1A=>0,:F3A=>0),toBeTrimmed=trues(2))

# Model
rightAngledFrameTrim = create_Model(name="rightAngledFrameTrim",beams=[beam1,beam2],BCs=[roller,pin,verticalForce,horizontalForce,balanceLoads])

# Set NR system solver
solver = create_NewtonRaphson(displayStatus=true)

# Create and solve the problem
problem = create_TrimProblem(model=rightAngledFrameTrim,systemSolver=solver,getInertiaMatrix=false)
solve!(problem)

# Get solution 
balanceHorizontalForce = problem.trimSolution.trimLoads[1]
balanceVerticalForce = problem.trimSolution.trimLoads[2]

# Compare to analytical solution 
balanceHorizontalForceAnalytical = -F
balanceVerticalForceAnalytical = -(F*L/2+F*L)/L

ϵ_rel_F1 = balanceHorizontalForce/balanceHorizontalForceAnalytical - 1
ϵ_rel_F3 = balanceVerticalForce/balanceVerticalForceAnalytical - 1

println("Relative errors:\nF1: $(round(100ϵ_rel_F1, digits=4)) % \nF3: $(round(100ϵ_rel_F3, digits=4)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/rightAngledFrameTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,legendPos=:bottomleft,save=true,savePath=string(relPath,"/rightAngledFrameTrim_deformation.pdf"))
display(deformationPlot)

println("Finished rightAngledFrameTrim.jl")
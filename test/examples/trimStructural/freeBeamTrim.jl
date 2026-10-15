using AeroBeams

# Beam
L = 1
EI = 1e3/48
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EI)])

# BCs
F = 1
trimFGuess = 0
midForce = create_BC(name="midForce",beam=beam,node=div(nElem,2)+1,conditions=Dict(:F3A=>-F))
trimTipLoad1 = create_BC(name="trimTipLoad1",beam=beam,node=1,conditions=Dict(:F3A=>trimFGuess),toBeTrimmed=trues(1))
trimTipLoad2 = create_BC(name="trimTipLoad2",beam=beam,node=nElem+1,conditions=Dict(:F3A=>trimFGuess),toBeTrimmed=trues(1))

# Trim link
trimLink = create_TrimLoadsLink(masterBC=trimTipLoad1,slaveBCs=[trimTipLoad2])

# Model
freeBeamTrim = create_Model(name="freeBeamTrim",beams=[beam],BCs=[midForce,trimTipLoad1,trimTipLoad2],trimLoadsLinks=[trimLink])

# Set up a nonlinear system solver
solver = create_NLsolver(method=:newton,linesearch=false,displayStatus=true)

# Create and solve the problem
problem = create_TrimProblem(model=freeBeamTrim,systemSolver=solver,getInertiaMatrix=false)
solve!(problem)

# Get solution
hist = problem.loadstepHistory
x1 = reduce(vcat,[beam.elements[e].geometry.x1_n1,
                  beam.elements[e].geometry.x1_n2] for e in 1:nElem)
F3 = reduce(vcat,[hist.states[end][e].nodal.F1[3],
                  hist.states[end][e].nodal.F2[3]] for e in 1:nElem)
M2 = reduce(vcat,[hist.states[end][e].nodal.M1[2],
                  hist.states[end][e].nodal.M2[2]] for e in 1:nElem)
trimForce = problem.trimSolution.trimLoads[1]

# Compare to analytical solution 
trimForceAnalytical = F/2
ϵ_rel = trimForce/trimForceAnalytical - 1
println("Trim force relative error: $(round(100ϵ_rel, digits=4)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/freeBeamTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,scale=1e2,showScale=true,scalePos=[0.2,0.22],save=true,savePath=string(relPath,"/freeBeamTrim_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# F3
plt1 = plot(x1/L, F3/F, lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3/F\$")
display(plt1)
savefig(string(absPath,"/freeBeamTrim_F.pdf"))

# M2
plt2 = plot(x1/L, M2/(F*L/4), lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2/(FL/4)\$")
display(plt2)
savefig(string(absPath,"/freeBeamTrim_M.pdf"))

println("Finished freeBeamTrim.jl")
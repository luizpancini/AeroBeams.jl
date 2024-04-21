using AeroBeams, LinearAlgebra, Plots

# Beam
L = 1
EI = 1e3/48
∞ = 1e14
stiffnessMatrix = diagm([∞,∞,∞,∞,EI,∞])
nElem = 30
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix])

# BCs
F = 1
trimFGuess = 0.
midForce = create_BC(name="midForce",beam=beam,node=div(nElem,2)+1,types=["F3A"],values=[-F])
rollerTip1 = create_BC(name="rollerTip1",beam=beam,node=1,types=["u3A"],values=[0])
pinTip2 = create_BC(name="pinTip2",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
trimTipLoad1 = create_BC(name="trimTipLoad1",beam=beam,node=1,types=["F3A"],values=[trimFGuess],toBeTrimmed=trues(1))
trimTipLoad2 = create_BC(name="trimTipLoad2",beam=beam,node=nElem+1,types=["F3A"],values=[trimFGuess],toBeTrimmed=trues(1))

# Trim link
trimLink = create_TrimLoadsLink(masterBC=trimTipLoad1,slaveBCs=[trimTipLoad2])

# Model
freeBeamTrim = create_Model(name="freeBeamTrim",beams=[beam],BCs=[midForce,trimTipLoad1,trimTipLoad2,rollerTip1,pinTip2],trimLoadsLinks=[trimLink])

# Set NR system solver with increased number of maximum iterations
NR = create_NewtonRaphson(maximumIterations=50,displayStatus=true)

# Create and solve the problem
problem = create_TrimProblem(model=freeBeamTrim,systemSolver=NR)
solve!(problem)

# Get solution 
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)
trimForce = problem.x[end]*problem.model.forceScaling

# Compare to analytical solution 
trimForceAnalytical = F/2
ϵ_rel = trimForce/trimForceAnalytical - 1

println("Trim force relative error: $ϵ_rel")

# Plots
plt2 = plot(x1/L, F3/F, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3/F\$ [N]")
display(plt2)
plt3 = plot(x1/L, M2/(F*L/4), lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2/(FL/4)\$")
display(plt3)

println("Finished freeBeamTrim.jl")
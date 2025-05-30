using AeroBeams, LinearAlgebra

# Beam
L = 1
EI = 1e3/48
∞ = 1e14
stiffnessMatrix = diagm([∞,∞,∞,∞,EI,∞])
nElem = 30
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix])

# BCs
F = 1
u3_mid = -F*L/(48*EI)
midDisp = create_BC(name="midDisp",beam=beam,node=div(nElem,2)+1,types=["u3A"],values=[u3_mid])
pin = create_BC(name="pin",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
roller = create_BC(name="roller",beam=beam,node=nElem+1,types=["u2A","u3A","p1A","p3A"],values=[0,0,0,0])
rollerReaction = create_BC(name="rollerReaction",beam=beam,node=nElem+1,types=["F3A"],values=[0],toBeTrimmed=[true])

# Model
midLoadedBeamTrim = create_Model(name="midLoadedBeamTrim",beams=[beam],BCs=[pin,roller,midDisp,rollerReaction])

# Set NR system solver with increased number of maximum iterations
NR = create_NewtonRaphson(maximumIterations=50,displayStatus=true)

# Create and solve the problem
problem = create_TrimProblem(model=midLoadedBeamTrim,systemSolver=NR)
solve!(problem)

# Get solution 
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)
roller_F3 = problem.nodalStatesOverσ[end][nElem].F_n2[3] 

# Compare to analytical solution (valid for small displacements)
roller_F3_analytical = F/2
ϵ_rel = roller_F3/roller_F3_analytical - 1
println("Relative errors:\nF3: $ϵ_rel")

println("Finished midLoadedBeamTrim.jl")
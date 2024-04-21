using AeroBeams, LinearAlgebra, Plots

# Beam
L = 1
EI = 333.333
∞ = 1e14
stiffnessMatrix = diagm([∞,∞,∞,∞,EI,∞])
nElem = 30
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix])

# BCs
F = 1
trimF3Guess,trimM2Guess = 0,0
tipForce = create_BC(name="tipForce",beam=beam,node=1,types=["F3A"],values=[-F])
clamp = create_BC(name="clamp",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clampReactions = create_BC(name="clampReactions",beam=beam,node=nElem+1,types=["F3A","M2A"],values=[trimF3Guess,trimM2Guess],toBeTrimmed=[true,true])

# Model
tipLoadedCantileverTrim = create_Model(name="tipLoadedCantileverTrim",beams=[beam],BCs=[clamp,tipForce,clampReactions])

# Set NR system solver with increased number of maximum iterations
NR = create_NewtonRaphson(maximumIterations=50,displayStatus=true)

# Create and solve the problem
problem = create_TrimProblem(model=tipLoadedCantileverTrim,systemSolver=NR)
solve!(problem)

# Get solution 
x1 = vcat([vcat(problem.model.beams[1].elements[e].x1_n1,problem.model.beams[1].elements[e].x1_n2) for e in 1:nElem]...)
u3 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
F3 = vcat([vcat(problem.nodalStatesOverσ[end][e].F_n1[3],problem.nodalStatesOverσ[end][e].F_n2[3]) for e in 1:nElem]...)
M2 = vcat([vcat(problem.nodalStatesOverσ[end][e].M_n1[2],problem.nodalStatesOverσ[end][e].M_n2[2]) for e in 1:nElem]...)
tip_u3 = problem.nodalStatesOverσ[end][1].u_n1[3] 
root_F3 = problem.nodalStatesOverσ[end][nElem].F_n2[3] 
root_M2 = problem.nodalStatesOverσ[end][nElem].M_n2[2]

# Compare to analytical solution (valid for small displacements)
tip_u3_analytical = -F*L/(3*EI)
root_F3_analytical = F
root_M2_analytical = F*L

ϵ_rel_u3 = tip_u3/tip_u3_analytical - 1
ϵ_rel_F3 = root_F3/root_F3_analytical - 1
ϵ_rel_M2 = root_M2/root_M2_analytical - 1

println("Relative errors:\nu3: $ϵ_rel_u3 \nF3: $ϵ_rel_F3 \nM2: $ϵ_rel_M2")

# Plots
plt1 = plot(x1/L, u3/(F*L/(3*EI)), lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3 / (FL/3EI)\$")
display(plt1)
plt2 = plot(x1/L, F3/F, lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3/F\$ [N]")
display(plt2)
plt3 = plot(x1/L, M2/(F*L), lw=2, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2/(FL)\$")
display(plt3)


println("Finished tipLoadedCantileverTrim.jl")
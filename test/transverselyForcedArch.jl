using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam 
R,θ = 0.5,π
L = R*θ
A,Iy = 1e-4,0.5e-8
E = 72e9
∞ = 1e12
EA,GAy,GAz,GJ,EIy,EIz = E*A,∞,∞,∞,E*Iy,∞
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
nElem = 20
beam = Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],rotationParametrization="E321",p0=[0;-π/2;0],k=[0;1/R;0])

# BCs
F = -6e3
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipFollowerForce = create_BC(name="tipFollowerForce",beam=beam,node=nElem+1,types=["Ff1A"],values=[F])

# Model
transverselyForcedArch = Model(name="transverselyForcedArch",beams=[beam],BCs=[clamp,tipFollowerForce])

# Set system solver options
σ0 = 0
σstep = 0.02
NR = NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = SteadyProblem(model=transverselyForcedArch,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u1 = [problem.nodalStatesOverσ[i][end].u_n2[1] for i in 1:length(σVector)]
tip_u3 = [problem.nodalStatesOverσ[i][end].u_n2[3] for i in 1:length(σVector)]
tip_angle = [problem.nodalStatesOverσ[i][end].θ_n2 for i in 1:length(σVector)]

# Plot normalized displacements over load steps
x = [-tip_u1/R, -tip_u3/R, -tip_angle/(π/2)]
labels = ["\$-u_1/R\$" "\$-u_3/R\$" "\$-\\theta/(\\pi/2)\$"]
XLabel = "\$-u_1/R, -u_3/R, -\\theta/(\\pi/2)\$"
colors = [:blue,:orange,:green]
plt1 = plot(xlabel=XLabel, ylabel="\$F\$ [kN]", title="Tip generalized displacements")
plot!(x, σVector*abs(F)/(1e3), linewidth=2, label=labels)
display(plt1)

println("Finished transverselyForcedArch.jl")
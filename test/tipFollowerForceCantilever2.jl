using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Beam
L = 100
E = 420e6
A,Iy = 1,1/12
∞ = 1e14
stiffnessMatrix = diagm([E*A,∞,∞,∞,E*Iy,∞])
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix])

# BCs - separate load in two parts
F = 130e3
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipFollowerForce1 = create_BC(name="tipFollowerForce1",beam=beam,node=nElem+1,types=["Ff3A"],values=[F/2])
tipFollowerForce2 = create_BC(name="tipFollowerForce2",beam=beam,node=nElem+1,types=["Ff3A"],values=[F/2])

# Model
tipFollowerForceCantilever = create_Model(name="tipFollowerForceCantilever",beams=[beam],BCs=[clamp,tipFollowerForce1,tipFollowerForce2])

# Set system solver options
σ0 = 0.0
σstep = 0.01
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=tipFollowerForceCantilever,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u1 = [problem.nodalStatesOverσ[i][nElem].u_n2[1] for i in 1:length(σVector)]
tip_u3 = [problem.nodalStatesOverσ[i][nElem].u_n2[3] for i in 1:length(σVector)]

# Load reference solution
u1_ref = readdlm(string(pwd(),"/test/referenceData/tipFollowerForceCantilever/u1.txt"))
u3_ref = readdlm(string(pwd(),"/test/referenceData/tipFollowerForceCantilever/u3.txt"))

# Plot configurations
colors = [:blue,:green]
labels = ["\$-u_1/L\$" "\$u_3/L\$"]
lw = 2
ms = 3
msw = 0
plt1 = plot(xlabel="\$F\$ [kN]", ylabel="\$-u_1/L, u_3/L\$", title="Tip generalized displacements", xticks=collect(0:10:F), yticks=collect(-0.6:0.2:1.2))
plot!([NaN], [NaN], lc=:black,  lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Simo & Vu-Quoc (1986)")
for i=1:2
    plot!([NaN], [NaN], lc=colors[i], m=colors[i],  lw=lw, ms=ms, msw=msw, label=labels[i])
end

# Plot normalized tip displacements over load steps
scatter!([u1_ref[1,:],u3_ref[1,:]], [u1_ref[2,:]/L,u3_ref[2,:]/L], palette=colors, ms=ms, msw=msw, label=false)
plot!(σVector*F/(1e3), [-tip_u1/L, tip_u3/L], palette=colors,  lw=lw, label=false)
display(plt1)

println("Finished tipFollowerForceCantilever.jl")
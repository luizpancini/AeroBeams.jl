using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Beam
R,θ = 100,π/4
L = R*θ
A,Iy,Iz,J = 1,1/12,1/12,1/6
E,G = 1e7,5e6
stiffnessMatrix = diagm([E*A,G*A,G*A,G*J,E*Iy,E*Iz])
nElem = 40
beam = Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],k=[0;0;1/R])

# BCs
F = 3000
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,types=["Ff3A"],values=[F])

# Model
curvedCantileverStaticFollower = Model(name="curvedCantileverStaticFollower",beams=[beam],BCs=[clamp,tipForce])

# Set system solver options
σ0 = 0.0
σstep = 1e-2
NR = NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = SteadyProblem(model=curvedCantileverStaticFollower,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u1 = [problem.nodalStatesOverσ[i][nElem].u_n2[1] for i in 1:length(σVector)]
tip_u2 = [problem.nodalStatesOverσ[i][nElem].u_n2[2] for i in 1:length(σVector)]
tip_u3 = [problem.nodalStatesOverσ[i][nElem].u_n2[3] for i in 1:length(σVector)]

# Reference solution digitalized from Simo and Vu-Quoc (1986)
u1_ref = readdlm(string(pwd(),"/test/referenceData/curvedCantileverStaticFollower/u1.txt"))
u2_ref = readdlm(string(pwd(),"/test/referenceData/curvedCantileverStaticFollower/u2.txt"))
u3_ref = readdlm(string(pwd(),"/test/referenceData/curvedCantileverStaticFollower/u3.txt"))

# Plot configurations
colors = [:blue,:green,:orange]
labels = ["\$-u_1\$" "\$u_2\$" "\$u_3\$"]
lw = 2
ms = 3
msw = 0
plt1 = plot(xlabel="\$F\$ [lb]", ylabel="\$-u_1, u_2, u_3\$ [in]", title="Tip displacements", xticks=collect(0:500:F), yticks=collect(-60:20:80))
plot!([NaN], [NaN], lc=:black,  lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Simo & Vu-Quoc (1986)")
for i=1:3
    plot!([NaN], [NaN], lc=colors[i], m=colors[i],  lw=lw, ms=ms, msw=msw, label=labels[i])
end

# Plot normalized tip displacements over load steps
scatter!([u1_ref[1,:],u2_ref[1,:],u3_ref[1,:]], [u1_ref[2,:],u2_ref[2,:],u3_ref[2,:]], palette=colors, ms=ms, msw=msw, label=false)
plot!(σVector*F, [-tip_u1, tip_u2, tip_u3], palette=colors, lw=lw, label=false)
display(plt1)

println("Finished curvedCantileverStaticFollower.jl")

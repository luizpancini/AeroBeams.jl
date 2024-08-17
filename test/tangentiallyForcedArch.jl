using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Beam 
R,θ = 0.5,π
L = R*θ
A,Iy = 1e-4,0.5e-8
E = 72e9
∞ = 1e12
EA,GAy,GAz,GJ,EIy,EIz = E*A,∞,∞,∞,E*Iy,∞
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],rotationParametrization="E321",p0=[0;-π/2;0],k=[0;1/R;0])

# BCs
F = 2.5e3
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipFollowerForce = create_BC(name="tipFollowerForce",beam=beam,node=nElem+1,types=["Ff3A"],values=[F])

# Model
tangentiallyForcedArch = create_Model(name="tangentiallyForcedArch",beams=[beam],BCs=[clamp,tipFollowerForce])

# Set system solver options
σ0 = 0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=tangentiallyForcedArch,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u1 = [problem.nodalStatesOverσ[i][end].u_n2[1] for i in 1:length(σVector)]
tip_u3 = [problem.nodalStatesOverσ[i][end].u_n2[3] for i in 1:length(σVector)]
tip_angle = [problem.nodalStatesOverσ[i][end].θ_n2 for i in 1:length(σVector)]

# Load reference solution
u1Ref = readdlm(string(pwd(),"/test/referenceData/tangentiallyForcedArch/u1.txt"))
u3Ref = readdlm(string(pwd(),"/test/referenceData/tangentiallyForcedArch/u3.txt"))
θRef = readdlm(string(pwd(),"/test/referenceData/tangentiallyForcedArch/theta.txt"))

# Plot deformed shape
relPath = "/test/outputs/figures/tangentiallyForcedArch"
absPath = string(pwd(),relPath)
mkpath(absPath)
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/tangentiallyForcedArch_deformation.pdf"))
display(deformationPlot)

# Plot normalized displacements over load steps
gr()
lw = 2
ms = 4
msw = 0
x = [-tip_u1/R, tip_u3/R, -tip_angle/(π/2)]
labels = ["\$-u_1/R\$" "\$u_3/R\$" "\$-\\theta/(\\pi/2)\$"]
XLabel = "\$-u_1/R, u_3/R, -\\theta/(\\pi/2)\$"
colors = [:blue,:orange,:green]
plt1 = plot(xlabel=XLabel, ylabel="\$F\$ [kN]", title="Tip generalized displacements",legend=:bottomright)
plot!([NaN], [NaN], lc=:black,  lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Argyris & Symeonidis (1981)")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(x, σVector*abs(F)/(1e3), lw=lw,palette=colors,label=false)
scatter!([u1Ref[1,:],u3Ref[1,:],θRef[1,:]], [u1Ref[2,:]/1e3,u3Ref[2,:]/1e3,θRef[2,:]/1e3], palette=colors,ms=ms,msw=msw,label=false)
display(plt1)
savefig(string(absPath,"/tangentiallyForcedArch_disp.pdf"))

println("Finished tangentiallyForcedArch.jl")
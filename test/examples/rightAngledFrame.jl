using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam frame
L = 0.24
A,Iy = 0.18e-4,0.135e-8
E = 7.124e10
ν = 0.3
G = E/(2*(1+ν))
∞ = 1e12
EA,GAy,GAz,GJ,EIy,EIz = E*A,∞,∞,∞,E*Iy,∞
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
nElem = 20
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,C=[stiffnessMatrix])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,C=[stiffnessMatrix],rotationParametrization="E321",p0=[0;π/2;0])

# BCs
F = 5e3
clamp = create_BC(name="clamp",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
tipFollowerForce = create_BC(name="tipFollowerForce",beam=beam2,node=nElem+1,types=["Ff3b"],values=[-F])

# Model
rightAngledFrame = create_Model(name="rightAngledFrame",beams=[beam1,beam2],BCs=[clamp,tipFollowerForce])

# Set system solver options
σ0 = 0.0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=rightAngledFrame,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
tip_u1 = [problem.nodalStatesOverσ[i][end].u_n2[1] for i in 1:length(σVector)]
tip_u3 = [problem.nodalStatesOverσ[i][end].u_n2[3] for i in 1:length(σVector)]
tip_angle = [problem.nodalStatesOverσ[i][end].θ_n2 for i in 1:length(σVector)]

# Load reference solution
u1Ref = readdlm(string(pwd(),"/test/referenceData/rightAngledFrame/u1.txt"))
u3Ref = readdlm(string(pwd(),"/test/referenceData/rightAngledFrame/u3.txt"))
θRef = readdlm(string(pwd(),"/test/referenceData/rightAngledFrame/theta.txt"))

# Plot deformed shape
relPath = "/test/outputs/figures/rightAngledFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)
deformationPlot = plot_steady_deformation(problem,save=true,legendPos=:bottomleft,savePath=string(relPath,"/rightAngledFrame_deformation.pdf"))
display(deformationPlot)

# Plot normalized displacements over load steps
lw = 2
ms = 4
gr()
x = [-tip_u1/L, tip_u3/L, tip_angle/π]
labels = ["\$-u_1/L\$" "\$u_3/L\$" "\$-\\theta/\\pi\$"]
colors = [:blue,:orange,:green]
plt1 = plot(xlabel="\$-u_1/L, u_3/L, -\\theta/\\pi\$", ylabel="\$F\$ [N]", title="Tip generalized displacements")
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Argyris & Symeonidis (1981)")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(x, σVector*F, lw=lw,palette=colors,label=false)
scatter!([u1Ref[1,:],u3Ref[1,:],θRef[1,:]], [u1Ref[2,:],u3Ref[2,:],θRef[2,:]], palette=colors,ms=ms,msw=msw,label=false)
display(plt1)
savefig(string(absPath,"/rightAngledFrame_disp.pdf"))

println("Finished rightAngledFrame.jl")
using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Beam frame
L = 120
A,Iy = 6,2
E = 7.2e6
ν = 0.0
G = E/(2*(1+ν))
∞ = 1e10
EA,GAy,GAz,GJ,EIy,EIz = E*A,∞,∞,∞,E*Iy,∞
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
nElem = 20
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,C=[stiffnessMatrix],rotationParametrization="E321",p0=[0;-π/2;0])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,C=[stiffnessMatrix])

# BCs
F = 35e3
elemForce = div(nElem,5)
support1 = create_BC(name="support1",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
support2 = create_BC(name="support2",beam=beam2,node=nElem+1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
force = create_BC(name="force",beam=beam2,node=elemForce+1,types=["Ff3A"],values=[-F])

# Model
LeeFrameFollowerLoad = create_Model(name="LeeFrameFollowerLoad",beams=[beam1,beam2],BCs=[support1,support2,force],units=create_UnitsSystem(length="in",force="lbf"))

# Set system solver options
σ0 = 0.0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=LeeFrameFollowerLoad,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
u1_atForce = [problem.nodalStatesOverσ[i][nElem+elemForce].u_n2[1] for i in 1:length(σVector)]
u3_atForce = [problem.nodalStatesOverσ[i][nElem+elemForce].u_n2[3] for i in 1:length(σVector)]

# Load reference solution
u1Ref = readdlm(string(pwd(),"/test/referenceData/LeeFrameFollowerLoad/u1.txt"))
u3Ref = readdlm(string(pwd(),"/test/referenceData/LeeFrameFollowerLoad/u3.txt"))

# Plot deformed shape
relPath = "/test/outputs/figures/LeeFrameFollowerLoad"
absPath = string(pwd(),relPath)
mkpath(absPath)
deformationPlot = plot_steady_deformation(problem,legendPos=:bottomright,save=true,savePath=string(relPath,"/LeeFrameFollowerLoad_deformation.pdf"))
display(deformationPlot)

# Plot normalized displacements over load steps
gr()
lw = 2
ms = 4
x = [u1_atForce/L, -u3_atForce/L]
labels = ["\$u_1/L\$" "\$-u_3/L\$"]
colors = [:blue,:orange]
plt1 = plot(xlabel="\$u_1/L, -u_3/L,\$", ylabel="\$F\$ [kip]", title="Displacements at point of force application")
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Simo and Vu-Quoc (1986)")
for i=1:2
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(x, σVector*F/(1e3), lw=lw,palette=colors,label=false)
scatter!([u1Ref[1,:],u3Ref[1,:]]/L, [u1Ref[2,:],u3Ref[2,:]], palette=colors,ms=ms,msw=msw,label=false)
display(plt1)
savefig(string(absPath,"/LeeFrameFollowerLoad_disp.pdf"))

println("Finished LeeFrameFollowerLoad.jl")
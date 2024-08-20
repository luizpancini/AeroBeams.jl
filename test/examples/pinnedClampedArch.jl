using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam 
R,θ= 100,215*π/180
L = R*θ
A,Iy = 2,2/3
E = 1.5e6
∞ = 1e10
EA,GAy,GAz,GJ,EIy,EIz = E*A,∞,∞,∞,E*Iy,∞
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
nElem = 80
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],rotationParametrization="E321",p0=[0;-θ/2;0],k=[0;1/R;0])

# BCs
λ = 8.9
elemForce = div(nElem,2)
hinge = create_BC(name="hinge",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
clamp = create_BC(name="clamp",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
force = create_BC(name="force",beam=beam,node=elemForce+1,types=["F3A"],values=[-λ*E*Iy/R^2])

# Model
pinnedClampedArch = create_Model(name="pinnedClampedArch",beams=[beam],BCs=[hinge,clamp,force],units=create_UnitsSystem(length="in",force="lbf"))

# Set system solver options
σ0 = 0.0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=pinnedClampedArch,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
u1_atForce = [problem.nodalStatesOverσ[i][elemForce].u_n2[1] for i in 1:length(σVector)]
u3_atForce = [problem.nodalStatesOverσ[i][elemForce].u_n2[3] for i in 1:length(σVector)]

# Load reference solution
u1Ref = readdlm(string(pwd(),"/test/referenceData/pinnedClampedArch/u1.txt"))
u3Ref = readdlm(string(pwd(),"/test/referenceData/pinnedClampedArch/u3.txt"))

# Plot deformed shape
relPath = "/test/outputs/figures/pinnedClampedArch"
absPath = string(pwd(),relPath)
mkpath(absPath)
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/pinnedClampedArch_deformation.pdf"))
display(deformationPlot)

# Plot normalized displacements over load steps
gr()
lw = 2
ms = 4
x = [-u1_atForce/R, -u3_atForce/R]
labels = ["\$-u_1/R\$" "\$-u_3/R\$"]
colors = [:blue,:orange]
plt1 = plot(xlabel="\$-u_1/L, -u_3/L,\$", ylabel="\$\\lambda\$", title="Displacements at point of force application",legend=:bottomright)
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="DaDeppo & Schmidt (1975)")
for i=1:2
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(x, σVector*λ, lw=lw,palette=colors,label=false)
scatter!([u1Ref[1,:],u3Ref[1,:]], [u1Ref[2,:],u3Ref[2,:]], palette=colors,ms=ms,msw=msw,label=false)
display(plt1)
savefig(string(absPath,"/pinnedClampedArch_disp.pdf"))

println("Finished pinnedClampedArch.jl")
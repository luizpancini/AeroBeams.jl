using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Choose force type (dead or follower)
forceType = "follower"

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
elemForce = div(nElem,5)
support1 = create_BC(name="support1",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
support2 = create_BC(name="support2",beam=beam2,node=nElem+1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
if forceType == "dead"
    F = 18e3
    force = create_BC(name="force",beam=beam2,node=elemForce+1,types=["F3A"],values=[-F])
elseif forceType == "follower"
    F = 35e3
    force = create_BC(name="force",beam=beam2,node=elemForce+1,types=["Ff3A"],values=[-F])
else
    error("Wrong force type")
end

# Model
LeeFrame = create_Model(name="LeeFrame",beams=[beam1,beam2],BCs=[support1,support2,force],units=create_UnitsSystem(length="in",force="lbf"))

# Set system solver options
σ0 = 0.0
σstep = 0.02
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=LeeFrame,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
σVector = problem.savedσ
u1_atForce = [problem.nodalStatesOverσ[i][nElem+elemForce].u_n2[1] for i in 1:length(σVector)]
u3_atForce = [problem.nodalStatesOverσ[i][nElem+elemForce].u_n2[3] for i in 1:length(σVector)]

# Plot deformed shape
relPath = "/test/outputs/figures/hingedSpringedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)
deformationPlot = plot_steady_deformation(problem,legendPos=:bottomright,save=true,savePath=string(relPath,"/LeeFrame_deformation.pdf"))
display(deformationPlot)

# Plot normalized displacements over load steps
gr()
x = [u1_atForce/L, -u3_atForce/L]
labels = ["\$u_1/L\$" "\$-u_3/L\$"]
colors = [:blue,:orange]
plt1 = plot()
plot!(x, σVector*F/(1e3), linewidth=2, label=false, xlabel="\$u_1/L, -u_3/L,\$", ylabel="\$F\$ [kip]", title="Displacements at point of force application")
halfNσ = round(Int,length(σVector)/2)
annotate!(x[1][halfNσ], σVector[halfNσ]*F/(1e3), text(labels[1], :top, :left, colors[1]))
annotate!(x[2][halfNσ], σVector[halfNσ]*F/(1e3), text(labels[2], :top, :left, colors[2]))
display(plt1)
savefig(string(absPath,"/LeeFrame_disp.pdf"))

println("Finished LeeFrame.jl")
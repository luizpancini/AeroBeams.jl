using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Beam frame
L = 10
EA,GAy,GAz,GJ,EIy,EIz = 1e6,1e6,1e6,1e3,1e3,1e3
ρA,ρI = 1,10
stiffnessMatrix = diagm([EA,GAy,GAz,GJ,EIy,EIz])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
nElem = 20
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[π/2;0;0])

# BCs
F₀ = 50
F = t -> ifelse.(t.<=1, F₀*t/1, ifelse.(t.<=2, F₀*(2-t), 0.0))
clamp = create_BC(name="clamp",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
elbowForce = create_BC(name="elbowForce",beam=beam2,node=1,types=["F3A"],values=[t->F(t)])

# Model
elbowFrame = create_Model(name="elbowFrame",beams=[beam1,beam2],BCs=[clamp,elbowForce])

# Time variables
tf = 30
Δt = 5e-2

# Create and solve the problem
problem = create_DynamicProblem(model=elbowFrame,finalTime=tf,Δt=Δt)
solve!(problem)

# Get solution over time
t = problem.timeVector
u3_elbow = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]
u3_tip = [problem.nodalStatesOverTime[i][end].u_n2[3] for i in 1:length(t)]

# Load reference solution
u3ElbowRef = readdlm(string(pwd(),"/test/referenceData/elbowFrame/u3_elbow.txt"))
u3TipRef = readdlm(string(pwd(),"/test/referenceData/elbowFrame/u3_tip.txt"))

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 4
relPath = "/test/outputs/figures/elbowFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Animation
plot_dynamic_deformation(problem,plotFrequency=5,view=(30,30),plotLimits=[(0,L),(0,L),(-L/2,L/2)],save=true,savePath=string(relPath,"/elbowFrame_deformation.gif"),displayProgress=true)
# Plot displacements over time
gr()
y = [u3_elbow, u3_tip]
labels = ["Elbow" "Tip"]
colors = [:blue,:orange]
plt1 = plot( xlabel="\$t\$ [s]", ylabel="\$u_3\$ [in]", title="OOP displacements",legend=:bottomleft)
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Simo & Vu-Quoc (1987)")
for i=1:2
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(t, y, lw=lw, palette=colors, label=false)
scatter!([u3ElbowRef[1,:],u3TipRef[1,:]], [u3ElbowRef[2,:],u3TipRef[2,:]], palette=colors,ms=ms,msw=msw,label=false)
display(plt1)
savefig(string(absPath,"/elbowFrame_disp.pdf"))

println("Finished elbowFrame.jl")
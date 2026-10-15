using AeroBeams

# Beam frame
L = 10
EA = GA = 1e6
GJ = EI = 1e3
ρA,ρI = 1, 10
S = isotropic_stiffness_matrix(EA=EA,GAy=GA,GAz=GA,GJ=GJ,EIy=EI,EIz=EI)
I = inertia_matrix(ρA=ρA,ρIy=ρI,ρIz=ρI)
nElem = 20
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,S=[S],I=[I])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,S=[S],I=[I],p0=transform_to_WM([π/2;0;0],from=:E321))

# BCs
F₀ = 50
F = t -> ifelse.(t.<=1, F₀*t/1, ifelse.(t.<=2, F₀*(2-t), 0.0))
clamp = create_BC(name="clamp",beam=beam1,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
elbowForce = create_BC(name="elbowForce",beam=beam2,node=1,conditions=Dict(:F3A=>t->F(t)))

# Model
elbowFrame = create_Model(name="elbowFrame",beams=[beam1,beam2],BCs=[clamp,elbowForce])

# Time variables
tf = 30
Δt = 5e-2

# Create and solve the problem
problem = create_DynamicProblem(model=elbowFrame,finalTime=tf,Δt=Δt)
solve!(problem)

# Get solution over time
sol = problem.solutionHistory
t = sol.t
states = sol.states
u3_elbow = [states[i][nElem].nodal.u2[3] for i in 1:length(t)]
u3_tip = [states[i][end].nodal.u2[3] for i in 1:length(t)]

# Load reference solution
u3ElbowRef = readdlm(pkgdir(AeroBeams,"test","referenceData","elbowFrame", "u3_elbow.txt"))
u3TipRef = readdlm(pkgdir(AeroBeams,"test","referenceData","elbowFrame", "u3_tip.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/elbowFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,plotFrequency=5,view=(30,30),plotLimits=([0,10/9*L],[-L/5,L],[-L/2-L/10,L/2+L/10]),save=true,savePath=string(relPath,"/elbowFrame_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
y = [u3_elbow, u3_tip]
labels = ["Elbow" "Tip"]
colors = [:blue,:orange]
lw = 2
ms = 4
msw = 0
gr()

# Plot displacements over time
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
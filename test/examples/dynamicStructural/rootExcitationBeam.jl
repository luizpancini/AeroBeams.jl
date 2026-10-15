using AeroBeams

# Frequency to excite
frequencyMode = 1

# Beam
L,b,H = 479e-3, 50.8e-3, 0.45e-3
A,Iy,Iz = b*H, b*H^3/12, H*b^3/12
J = Iy+Iz
Ksy = Ksz = 5/6
E,ν,ρ = 127e9, 0.36, 4.43e3
G = E/(2*(1+ν))
nElem = 50
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=E*A,GAy=G*A*Ksy,GAz=G*A*Ksz,GJ=G*J,EIy=E*Iy,EIz=E*Iz)],I=[inertia_matrix(ρA=ρ*A,ρIy=ρ*Iy,ρIz=ρ*Iz)],p0=transform_to_WM([0,-π/2,0],from=:E321))

# BCs
if frequencyMode == 1
    ω = Hz2rads(9)
    V = 0.1399
elseif frequencyMode == 2
    ω = Hz2rads(32)
    V = 0.3414
else
    error("Set frequencyMode as 1 or 2")
end
A = V/ω
T = 2π/ω
u₃b = t -> A*sin.(ω*t)
shaker = create_BC(name="shaker",beam=beam,node=1,conditions=Dict(:u1b=>0,:u2b=>0,:u3b=>t->u₃b(t),:p1b=>0,:p2b=>0,:p3b=>0))

# Model
rootExcitationBeam = create_Model(name="rootExcitationBeam",beams=[beam],BCs=[shaker],gravityVector=[0,0,-standard_atmosphere(0).g])

# Time variables
cycles = 5
tf = cycles * T
Δt = T/100

# Create and solve the problem
problem = create_DynamicProblem(model=rootExcitationBeam,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
x1 = [e.geometry.x1 for e in rootExcitationBeam.assemblyData.elements]
u3b_root = [states[i][1].nodal.u1_b[3] for i in 1:length(t)]
u3b_tip = [states[i][nElem].nodal.u2_b[3] for i in 1:length(t)]
V3_root = [states[i][1].midpoint.V[3] for i in 1:length(t)]
V3_tip = [states[i][nElem].midpoint.V[3] for i in 1:length(t)]
V3 = [[states[i][e].midpoint.V[3] for e = 1:nElem] for i in 1:length(t)]

## Plots
# Set paths
plotName = string("rootExcitationBeam",frequencyMode)
relPath = "/test/outputs/figures/$(plotName)"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,plotFrequency=5,scalePos=[1;-0.05;0],timeStampPos=[1;-0.15;0],plotLimits=([-0.1,0.1],[0,L],[0,L]),plotDistLoads=false,save=true,savePath=string(relPath,"/$(plotName)_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
ms = 3
gr()

# Normalized V3/V over the beam, over time
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$V_3/V_b\$")
for i in 1:length(t)
    plot!(x1/L,V3[i]/V, c=:black, lw=lw, label=false)
end
display(plt1)
savefig(string(absPath,"/$(plotName)_V3oV.pdf"))

# Normalized u3 (basis b) at the root and at the tip
plt2 = plot()
plot!(t/T,u3b_root/A, c=:auto, lw=lw, label="Root", xlabel="\$t/T\$")
plot!(t/T,u3b_tip/A, c=:auto, lw=lw, label="Tip", xlabel="\$t/T\$", ylabel="\$u_3^{+}/A\$")
display(plt2)
savefig(string(absPath,"/$(plotName)_u3.pdf"))

# Normalized V3 at the root and at the tip
plt3 = plot()
plot!(t/T,V3_root/V, c=:auto, lw=lw, label="Root", xlabel="\$t/T\$")
plot!(t/T,V3_tip/V, c=:auto, lw=lw, label="Tip", xlabel="\$t/T\$", ylabel="\$V_3^{*}/A\$")
display(plt3)
savefig(string(absPath,"/$(plotName)_V3.pdf"))

println("Finished rootExcitationBeam.jl")
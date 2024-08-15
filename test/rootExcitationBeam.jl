using AeroBeams, LinearAlgebra, Plots

# Choose excitation mode (1 or 2)
mode = 1

# Beam
L,b,H = 479e-3,50.8e-3,0.45e-3
A,Iy,Iz = b*H,b*H^3/12,H*b^3/12
J = Is = Iy+Iz
Ksy = Ksz = 5/6
E,ν,ρ = 127e9,0.36,4.43e3
G = E/(2*(1+ν))
nElem = 60
stiffnessMatrix = diagm([E*A,G*A*Ksy,G*A*Ksz,G*J,E*Iy,E*Iz])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,ρ*Is,ρ*Iy,ρ*Iz])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[0,-π/2,0])

# BCs
if mode == 1
    ω = 9*(2*π)
    V = 0.1399
elseif mode == 2
    ω = 32*(2*π)
    V = 0.3414
end
A = V/ω
T = 2*π/ω
u₃b = t -> A*sin.(ω*t)
shaker = create_BC(name="shaker",beam=beam,node=1,types=["u1b","u2b","u3b","p1b","p2b","p3b"],values=[0,0,t->u₃b(t),0,0,0])

# Model
rootExcitationBeam = create_Model(name="rootExcitationBeam",beams=[beam],BCs=[shaker],gravityVector=[0,0,-9.81])

# Time variables
cycles = 5
tf = cycles*T
Δt = T/100

# Create and solve the problem
problem = create_DynamicProblem(model=rootExcitationBeam,finalTime=tf,Δt=Δt)
solve!(problem)
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
x1 = [v.x1 for v in rootExcitationBeam.elements]
u3b_root = [problem.nodalStatesOverTime[i][1].u_n1_b[3] for i in 1:length(t)]
u3b_tip = [problem.nodalStatesOverTime[i][nElem].u_n2_b[3] for i in 1:length(t)]
V3_root = [problem.elementalStatesOverTime[i][1].V[3] for i in 1:length(t)]
V3_tip = [problem.elementalStatesOverTime[i][nElem].V[3] for i in 1:length(t)]
V3 = Vector{Vector{Float64}}()
for i in 1:length(t)
    push!(V3,[problem.elementalStatesOverTime[i][e].V[3] for e = 1:nElem])
end

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 3
relPath = "/test/outputs/figures/rootExcitationBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Animation
plot_dynamic_deformation(problem,plotFrequency=5,scalePos=[1;-0.05;0],timeStampPos=[1;-0.15;0],plotLimits=[(-0.1,0.1),(0,L),(0,L)],plotDistLoads=false,save=true,savePath=string(relPath,"/rootExcitationBeam_deformation.gif"),displayProgress=true)
# Normalized V3/V over the beam, over time
gr()
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$V_3/V_b\$")
for i in 1:length(t)
    plot!(x1/L,V3[i]/V, c=:black, lw=lw, label=false)
end
display(plt1)
savefig(string(absPath,"/rootExcitationBeam_V3oV.pdf"))
# Normalized u3 (basis b) at the root and at the tip
plt2 = plot()
plot!(t/T,u3b_root/A, c=:auto, lw=lw, label="Root", xlabel="\$t/T\$")
plot!(t/T,u3b_tip/A, c=:auto, lw=lw, label="Tip", xlabel="\$t/T\$", ylabel="\$u_3^{+}/A\$")
display(plt2)
savefig(string(absPath,"/rootExcitationBeam_u3.pdf"))
# Normalized V3 at the root and at the tip
plt3 = plot()
plot!(t/T,V3_root/V, c=:auto, lw=lw, label="Root", xlabel="\$t/T\$")
plot!(t/T,V3_tip/V, c=:auto, lw=lw, label="Tip", xlabel="\$t/T\$", ylabel="\$V_3^{*}/A\$")
display(plt3)
savefig(string(absPath,"/rootExcitationBeam_V3.pdf"))

println("Finished rootExcitationBeam.jl")
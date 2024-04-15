using AeroBeams, LinearAlgebra, Plots

# Choose excitation mode (1 or 2)
mode = 1

# Beam
L,b,H = 479e-3,50.8e-3,0.45e-3
A,Iy,Iz = b*H,b*H^3/12,H*b^3/12
J = Is = Iy+Iz
Ksy = Ksz = 5/6
E = 127e9
ν = 0.36
G = E/(2*(1+ν))
ρ = 4.43e3
nElem = 60
stiffnessMatrix = diagm([E*A,G*A*Ksy,G*A*Ksz,G*J,E*Iy,E*Iz])
inertiaMatrix = diagm([ρ*A,ρ*A,ρ*A,ρ*Is,ρ*Iy,ρ*Iz])
beam = Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[0,-π/2,0])

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
rootExcitationBeam = Model(name="rootExcitationBeam",beams=[beam],BCs=[shaker],gravityVector=[0,0,-9.81])

# Time variables
cycles = 5
tf = cycles*T
Δt = T/100

# Create and solve the problem
problem = DynamicProblem(model=rootExcitationBeam,finalTime=tf,Δt=Δt)
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
# Normalized V3/V over the beam, over time
plt1 = plot()
for i in 1:length(t)
    plot!(x1/L,V3[i]/V, c=:black, linewidth=2, label=false)
end
xlabel!("\$x_1/L\$")
ylabel!("\$V_3/V_b\$")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/rootExcitationBeam_1.pdf"))
# Normalized u3 (basis b) at the root and at the tip
plt2 = plot()
plot!(t/T,u3b_root/A, c=:auto, linewidth=2, label="Root", xlabel="\$t/T\$")
plot!(t/T,u3b_tip/A, c=:auto, linewidth=2, label="Tip", xlabel="\$t/T\$", ylabel="\$u_3^{+}/A\$")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/rootExcitationBeam_2.pdf"))
# Normalized V3 at the root and at the tip
plt3 = plot()
plot!(t/T,V3_root/V, c=:auto, linewidth=2, label="Root", xlabel="\$t/T\$")
plot!(t/T,V3_tip/V, c=:auto, linewidth=2, label="Tip", xlabel="\$t/T\$", ylabel="\$V_3^{*}/A\$")
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/rootExcitationBeam_3.pdf"))

println("Finished rootExcitationBeam.jl")
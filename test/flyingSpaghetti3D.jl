using AeroBeams, LinearAlgebra, Plots

# Beam
L = 10
EA,GA,GJ,EI = 1e4,1e4,5e2,5e2
ρA,ρI = 1,10
θ₀ = atan(4/3)
nElem = 10
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[0;θ₀;0])

# BCs - stable up to 18 s with M₀ = 100
M₀ = 200
τ = 2.5
M2 = t -> ifelse.(t.<=τ, M₀*t/τ, ifelse.(t.<=2*τ, 2*M₀*(1-t/(2*τ)), 0))
M3 = t -> M2(t)/2
F1 = t -> M2(t)/10
forces = create_BC(name="forces",beam=beam,node=nElem+1,types=["F1A","M2A","M3A"],values=[t->F1(t),t->M2(t),t->M3(t)])

# Model
flyingSpaghetti3D = create_Model(name="flyingSpaghetti3D",beams=[beam],BCs=[forces])

# Time variables
tf = 3.75
Δt = 1e-3

# Create and solve the problem
problem = create_DynamicProblem(model=flyingSpaghetti3D,finalTime=tf,Δt=Δt)
solve!(problem)
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u₁_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
u₂_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[2] for i in 1:length(t)]
u₃_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]
θ_tip = [problem.nodalStatesOverTime[i][nElem].θ_n2 for i in 1:length(t)]

# Plots
# ------------------------------------------------------------------------------
# Nomalized tip displacements
labels = ["\$u_1/L\$" "\$u_2/L\$" "\$u_3/L\$"]
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,[u₁_tip/L, u₂_tip/L, u₃_tip/L], lw=2, label=labels)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/flyingSpaghetti3D_1.pdf"))

# Nomalized tip angle
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$\\theta/(2\\pi)\$")
plot!(t,θ_tip/(2*π), lw=2, label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/flyingSpaghetti3D_2.pdf"))

println("Finished flyingSpaghetti3D.jl")
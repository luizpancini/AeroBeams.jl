using AeroBeams, LinearAlgebra, Plots

# Beam
L = 10
EA,GA,GJ,EI = 1e6,1e6,1e6,1e4
ρA,ρI1,ρI2 = 1,1,10
θ₀ = atan(4/3)
nElem = 20
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix1 = diagm([ρA,ρA,ρA,2*ρI1,ρI1,ρI1])
inertiaMatrix2 = diagm([ρA,ρA,ρA,2*ρI2,ρI2,ρI2])
inertiaMatrices = vcat([inertiaMatrix2 for _ in 1:div(nElem,2)],[inertiaMatrix1 for _ in 1:div(nElem,2)])
beam = Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=inertiaMatrices,rotationParametrization="E321",p0=[0;θ₀;0],hingedNodes=[div(nElem,2)+1],hingedNodesDoF=[[false,true,false]])

# BCs 
M₀ = 160
τ = 0.5
M2 = t -> ifelse.(t.<=τ, M₀, 0)
F1 = t -> M2(t)/4
forces = create_BC(name="forces",beam=beam,node=nElem+1,types=["F1A","M2A"],values=[t->F1(t),t->M2(t)])

# Model
flyingScissors = Model(name="flyingScissors",beams=[beam],BCs=[forces])

# Time variables
tf = 5
Δt = 5e-2

# Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2, Δt=Δt/10)

# Create and solve the problem
problem = DynamicProblem(model=flyingScissors,finalTime=tf,Δt=Δt,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
solve!(problem)
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u₁_tipA = [problem.nodalStatesOverTime[i][1].u_n2[1] for i in 1:length(t)]
u₃_tipA = [problem.nodalStatesOverTime[i][1].u_n2[3] for i in 1:length(t)]
u₁_tipB = [problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
u₃_tipB = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]
u₁_hinge = [problem.nodalStatesOverTime[i][div(nElem,2)].u_n2[1] for i in 1:length(t)]
u₃_hinge = [problem.nodalStatesOverTime[i][div(nElem,2)].u_n2[3] for i in 1:length(t)]

# Plots
# ------------------------------------------------------------------------------
# Nomalized tip displacements
labels = ["Tip A" "Hinge" "Tip B"]
plt1 = plot(xlabel="\$t\$ [s]", ylabel="\$u_1/L\$")
plot!(t,[u₁_tipA/L, u₁_hinge/L, u₁_tipB/L], lw=2, label=labels)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/flyingScissors_1.pdf"))
plt2 = plot(xlabel="\$t\$ [s]", ylabel="\$u_3/L\$")
plot!(t,[u₃_tipA/L, u₃_hinge/L, u₃_tipB/L], lw=2, label=labels)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/flyingScissors_2.pdf"))

println("Finished flyingScissors.jl")
using AeroBeams, LinearAlgebra, Plots, ColorSchemes

# Beam frame
L = 4*3.586301959539938
α = 1.032057439482664
β = 0.22606447252882428
γ = 1.3318e-01
stiffnessMatrix = diagm(1.0 ./ [2.93944738387698e-10, 8.42991725049126e-10, 3.38313996669689e-08, 4.69246721094557e-08, 6.79584100559513e-08, 1.37068861370898e-09])
inertiaMatrix = diagm([4.86e-2, 4.86e-2, 4.86e-2, 1.0632465e-2, 2.10195e-4, 1.042227e-2])
nElem = 8
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[α;-β;γ+π])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[-α;-β;-γ])
beam3 = create_Beam(name="beam3",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[-π+α;β;-γ])
beam4 = create_Beam(name="beam4",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[π-α;β;γ+π],connectedBeams=[beam1],connectedNodesThis=[nElem+1],connectedNodesOther=[1])

# BCs
Fₗ = 1e6
Fₛ = 1e4
τ = 0.04
F₁ = F₂ = t -> ifelse.(t.<=τ/4, Fₗ*t, ifelse.(t.<=τ/2, Fₗ*(τ/2-t), 0.0))
F₃ = t -> ifelse.(t.<=τ/2, Fₛ/2*(1-cos.(2*π*t/τ)), Fₛ)
clamp1 = create_BC(name="clamp1",beam=beam1,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=beam3,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
force1 = create_BC(name="force1",beam=beam2,node=1,types=["F1A","F2A","F3A"],values=[t->F₁(t),t->F₂(t),t->F₃(t)])
force2 = create_BC(name="force2",beam=beam4,node=1,types=["F1A","F2A","F3A"],values=[t->F₁(t),t->F₂(t),t->F₃(t)])

# Model
joinedBeams = create_Model(name="joinedBeams",beams=[beam1,beam2,beam3,beam4],BCs=[clamp1,clamp2,force1,force2])

# Time variables
tf = τ
Δt = τ/1e3

# Create and solve the problem
problem = create_DynamicProblem(model=joinedBeams,finalTime=tf,Δt=Δt)
solve!(problem)

# Get solution over time
t = problem.timeVector
u1 = [problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
u2 = [problem.nodalStatesOverTime[i][nElem].u_n2[2] for i in 1:length(t)]
u3 = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]

# Plot u1 of left tip over time
plt1 = plot(t, u1, linewidth=2, label=false, xlabel="\$t\$ [s]", ylabel="\$u_1\$ [m]")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/joinedBeams_1.pdf"))
# Plot u1 of left tip over time
plt2 = plot(t, u2, linewidth=2, label=false, xlabel="\$t\$ [s]", ylabel="\$u_2\$ [m]")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/joinedBeams_2.pdf"))
# Plot u1 of left tip over time
plt3 = plot(t, u3, linewidth=2, label=false, xlabel="\$t\$ [s]", ylabel="\$u_3\$ [m]")
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/joinedBeams_3.pdf"))

println("Finished joinedBeams.jl")
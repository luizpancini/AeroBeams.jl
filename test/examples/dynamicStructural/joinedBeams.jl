using AeroBeams

# Beam frame
L = 4 * 3.586301959539938
α = 1.032057439482664
β = 0.22606447252882428
γ = 1.3318e-01
EA = 3.4020000000171847e9
GAy = 1.186251264734211e9
GAz = 2.955833958523284e7
GJ = 2.1310750934336137e7
EIy = 1.4714882222475234e7
EIz = 7.295603027547412e8
ρA = 4.86e-2
ρIy = 2.10195e-4
ρIz = 1.042227e-2
S = isotropic_stiffness_matrix(EA=EA,GAy=GAy,GAz=GAz,GJ=GJ,EIy=EIy,EIz=EIz)
I = inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz)
nElem = 8
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,S=[S],I=[I],p0=transform_to_WM([α,-β,γ+π],from=:E321))
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,S=[S],I=[I],p0=transform_to_WM([-α,-β,-γ],from=:E321))
beam3 = create_Beam(name="beam3",length=L,nElements=nElem,S=[S],I=[I],p0=transform_to_WM([-π+α,β,-γ],from=:E321))
beam4 = create_Beam(name="beam4",length=L,nElements=nElem,S=[S],I=[I],p0=transform_to_WM([π-α,β,γ+π],from=:E321))

# BCs
Fₗ = 1e6
Fₛ = 1e4
τ = 0.04
F₁ = F₂ = t -> ifelse.(t.<=τ/4, Fₗ*t, ifelse.(t.<=τ/2, Fₗ*(τ/2-t), 0.0))
F₃ = t -> ifelse.(t.<=τ/2, Fₛ/2*(1-cos.(2*π*t/τ)), Fₛ)
clamp1 = create_BC(name="clamp1",beam=beam1,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
clamp2 = create_BC(name="clamp2",beam=beam3,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
force1 = create_BC(name="force1",beam=beam2,node=1,conditions=Dict(:F1A=>t->F₁(t),:F2A=>t->F₂(t),:F3A=>t->F₃(t)))
force2 = create_BC(name="force2",beam=beam4,node=1,conditions=Dict(:F1A=>t->F₁(t),:F2A=>t->F₂(t),:F3A=>t->F₃(t)))

# Model
joinedBeams = create_Model(name="joinedBeams",beams=[beam1,beam2,beam3,beam4],BCs=[clamp1,clamp2,force1,force2],connections=[(beam4=>beam1,nElem+1=>1)])

# Time variables
tf = τ
Δt = τ/1e3

# Create and solve the problem
problem = create_DynamicProblem(model=joinedBeams,finalTime=tf,Δt=Δt)
solve!(problem)

# Get solution over time
sol = problem.solutionHistory
t = sol.t
states = sol.states
u1 = [states[i][nElem].nodal.u2[1] for i in 1:length(t)]
u2 = [states[i][nElem].nodal.u2[2] for i in 1:length(t)]
u3 = [states[i][nElem].nodal.u2[3] for i in 1:length(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/joinedBeams"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,scale=100,plotFrequency=10,view=(90,0),plotLimits=([-L,L],[-L,L],[-L,L]),save=true,savePath=string(relPath,"/joinedBeams_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
gr()

# Plot u1 of left tip over time
plt1 = plot(t, u1, lw=lw, label=false, xlabel="\$t\$ [s]", ylabel="\$u_1\$ [m]")
display(plt1)
savefig(string(absPath,"/joinedBeams_u1.pdf"))

# Plot u2 of left tip over time
plt2 = plot(t, u2, lw=lw, label=false, xlabel="\$t\$ [s]", ylabel="\$u_2\$ [m]")
display(plt2)
savefig(string(absPath,"/joinedBeams_u2.pdf"))

# Plot u3 of left tip over time
plt3 = plot(t, u3, lw=lw, label=false, xlabel="\$t\$ [s]", ylabel="\$u_3\$ [m]")
display(plt3)
savefig(string(absPath,"/joinedBeams_u3.pdf"))

println("Finished joinedBeams.jl")
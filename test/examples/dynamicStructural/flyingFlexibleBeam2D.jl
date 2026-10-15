using AeroBeams

# Beam
L = 10
EA = GA = 1e4 
GJ = EI = 5e2
ρA,ρI = 1, 10
θ₀ = atan(4/3)
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA,GAy=GA,GAz=GA,GJ=GJ,EIy=EI,EIz=EI)],I=[inertia_matrix(ρA=ρA,ρIy=ρI,ρIz=ρI)],p0=transform_to_WM([0,θ₀,0],from=:E321))

# BCs 
M₀ = 80
τ = 2.5
M2 = t -> ifelse.(t.<=τ, M₀, 0)
F1 = t -> M2(t)/10
loads = create_BC(name="loads",beam=beam,node=nElem+1,conditions=Dict(:F1A=>t->F1(t),:M2A=>t->M2(t)))

# Model
flyingFlexibleBeam2D = create_Model(name="flyingFlexibleBeam2D",beams=[beam],BCs=[loads])

# Time variables
tf = 13
Δt = 1e-2

# Create and solve the problem
problem = create_DynamicProblem(model=flyingFlexibleBeam2D,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
u1_tip = [states[i][nElem].nodal.u2[1] for i in 1:length(t)]
u3_tip = [states[i][nElem].nodal.u2[3] for i in 1:length(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/flyingFlexibleBeam2D"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,plotFrequency=10,plotLimits=([0,3*L],[-L,L],[-L,L/2]),save=true,savePath=string(relPath,"/flyingFlexibleBeam2D_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
labels = ["\$u_1/L\$" "\$u_3/L\$"]
gr()

# Nomalized tip displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,[u1_tip/L, u3_tip/L], lw=lw, label=labels)
display(plt1)
savefig(string(absPath,"/flyingFlexibleBeam2D_disp.pdf"))

println("Finished flyingFlexibleBeam2D.jl")
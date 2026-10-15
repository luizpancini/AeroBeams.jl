using AeroBeams

# Beam
L = 10
EA = GA = 1e4
GJ,EI = 5e2, 1e2
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
flyingSpaghetti2D = create_Model(name="flyingSpaghetti2D",beams=[beam],BCs=[loads])

# Time variables
tf = 7.5
Δt = 1e-2

# Create and solve the problem
problem = create_DynamicProblem(model=flyingSpaghetti2D,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
u1_tip = [states[i][nElem].nodal.u2[1] for i in 1:length(t)]
u3_tip = [states[i][nElem].nodal.u2[3] for i in 1:length(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/flyingSpaghetti2D"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,plotFrequency=10,plotLimits=([0,2*L],[-L/2,L/2],[-L,0]),save=true,savePath=string(relPath,"/flyingSpaghetti2D_deformation.gif"),displayProgress=true,showScale=false,showTimeStamp=false,plotGrid=false,plotAxes=false)
display(anim)

# Snapshots
plt_snap = plot_snapshots(problem,plotBCs=false,plotAxes=false,plotGrid=false,snapshots=vcat(0:0.5:7.5),plotLimits=([0,2*L],[-L/2,L/2],[-L,0]),save=true,savePath=string(relPath,"/flyingSpaghetti2D_snapshots.pdf"))
display(plt_snap)

# Plot configurations
labels = ["\$u_1/L\$" "\$u_3/L\$"]
lw = 2
gr()

# Nomalized tip displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,[u1_tip/L, u3_tip/L], lw=lw, label=labels)
display(plt1)
savefig(string(absPath,"/flyingSpaghetti2D_disp.pdf"))

println("Finished flyingSpaghetti2D.jl")
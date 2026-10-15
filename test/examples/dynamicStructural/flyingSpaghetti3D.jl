using AeroBeams

# Beam
L = 10
EA = GA = 1e4
GJ = EI = 5e2
ρA,ρI = 1, 10
θ₀ = atan(4/3)
nElem = 10
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA,GAy=GA,GAz=GA,GJ=GJ,EIy=EI,EIz=EI)],I=[inertia_matrix(ρA=ρA,ρIy=ρI,ρIz=ρI)],p0=transform_to_WM([0,θ₀,0],from=:E321))

# BCs - stable for much longer with M₀ = 100 than with M₀ = 200
M₀ = 200
τ = 2.5
M2 = t -> ifelse.(t.<=τ, M₀*t/τ, ifelse.(t.<=2*τ, 2*M₀*(1-t/(2*τ)), 0))
M3 = t -> M2(t)/2
F1 = t -> M2(t)/10
loads = create_BC(name="loads",beam=beam,node=nElem+1,conditions=Dict(:F1A=>t->F1(t),:M2A=>t->M2(t),:M3A=>t->M3(t)))

# Model
flyingSpaghetti3D = create_Model(name="flyingSpaghetti3D",beams=[beam],BCs=[loads])

# System solver
systemSolver = create_NLsolver(method=:newton)

# Time variables
tf = 10
Δt = 5e-3

# Create and solve the problem
problem = create_DynamicProblem(model=flyingSpaghetti3D,finalTime=tf,Δt=Δt,systemSolver=systemSolver)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
u1_tip = [states[i][nElem].nodal.u2[1] for i in eachindex(t)]
u2_tip = [states[i][nElem].nodal.u2[2] for i in eachindex(t)]
u3_tip = [states[i][nElem].nodal.u2[3] for i in eachindex(t)]
p1_tip = [states[i][nElem].nodal.p2[1] for i in eachindex(t)]
p2_tip = [states[i][nElem].nodal.p2[2] for i in eachindex(t)]
p3_tip = [states[i][nElem].nodal.p2[3] for i in eachindex(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/flyingSpaghetti3D"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,view=(30,30),plotFrequency=10,plotLimits=([0,2*L],[-L/2,L/2],[-L,L]),save=true,savePath=string(relPath,"/flyingSpaghetti3D_deformation.gif"),displayProgress=true,showScale=false)
display(anim)

# Plot configurations
ulabels = ["\$u_1/L\$" "\$u_2/L\$" "\$u_3/L\$"]
plabels = ["\$p_1\$" "\$p_2\$" "\$p_3\$"]
lw = 2
gr()

# Nomalized tip displacements
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,[u1_tip/L, u2_tip/L, u3_tip/L], lw=lw, label=ulabels)
display(plt1)
savefig(string(absPath,"/flyingSpaghetti3D_disp.pdf"))

# Tip rotation parameters
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Tip rotations")
plot!(t,[p1_tip, p2_tip, p3_tip], lw=lw, label=plabels)
display(plt2)
savefig(string(absPath,"/flyingSpaghetti3D_rot.pdf"))

println("Finished flyingSpaghetti3D.jl")
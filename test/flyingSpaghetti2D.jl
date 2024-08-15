using AeroBeams, LinearAlgebra, Plots

# Beam
L = 10
EA,GA,GJ,EI = 1e4,1e4,5e2,1e2
ρA,ρI = 1,10
θ₀ = atan(4/3)
nElem = 20
stiffnessMatrix = diagm([EA,GA,GA,GJ,EI,EI])
inertiaMatrix = diagm([ρA,ρA,ρA,2*ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix],rotationParametrization="E321",p0=[0;θ₀;0])

# BCs 
M₀ = 80
τ = 2.5
M2 = t -> ifelse.(t.<=τ, M₀, 0)
F1 = t -> M2(t)/10
forces = create_BC(name="forces",beam=beam,node=nElem+1,types=["F1A","M2A"],values=[t->F1(t),t->M2(t)])

# Model
flyingSpaghetti2D = create_Model(name="flyingSpaghetti2D",beams=[beam],BCs=[forces])

# Time variables
tf = 7.5
Δt = 1e-2

# Create and solve the problem
problem = create_DynamicProblem(model=flyingSpaghetti2D,finalTime=tf,Δt=Δt)
solve!(problem)
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
u₁_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
u₃_tip = [problem.nodalStatesOverTime[i][nElem].u_n2[3] for i in 1:length(t)]

# Plots
# ------------------------------------------------------------------------------
lw = 2
relPath = "/test/outputs/figures/flyingSpaghetti2D"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Animation
plot_dynamic_deformation(problem,refBasis="I",plotFrequency=10,plotLimits=[(0,2*L),(-L,0),(-L/2,L/2)],save=true,savePath=string(relPath,"/flyingSpaghetti2D_deformation.gif"),displayProgress=true)
# Nomalized tip displacements
gr()
labels = ["\$u_1/L\$" "\$u_3/L\$"]
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip normalized displacements")
plot!(t,[u₁_tip/L, u₃_tip/L], lw=lw, label=labels)
display(plt1)
savefig(string(absPath,"/flyingSpaghetti2D_disp.pdf"))

println("Finished flyingSpaghetti2D.jl")
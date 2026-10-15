using AeroBeams
 
# Beam
L = 0.5
E = 2e10
ρ = 8e9
A = 1e-6
∞ = 1e12
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=E*A,∞=∞)],I=[inertia_matrix(ρA=ρ*A)])

# Time variables
tf = 1.0
Δt = 2e-3

# BCs
F₀ = 1.0
F = t -> ifelse(t>=Δt, F₀, 0.0)
axialForce = create_BC(name="axialForce",beam=beam,node=nElem+1,conditions=Dict(:F1A=>t->F(t)))
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
axialTractionCantilever = create_Model(name="axialTractionCantilever",beams=[beam],BCs=[clamp,axialForce])

# System solver
systemSolver = create_NewtonRaphson()
systemSolver = create_NLsolver()

# Create and solve the problem
problem = create_DynamicProblem(model=axialTractionCantilever,finalTime=tf,Δt=Δt,displayProgress=false,systemSolver=systemSolver)
# solve!(problem)
# @profview_allocs solve!(problem) sample_rate=0.01
@time solve!(problem)
# using Profile
# Profile.clear()
# Profile.init(n=10^7)
# @profview solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
x1 = reduce(vcat,[beam.elements[e].geometry.x1_n1,
                  beam.elements[e].geometry.x1_n2] for e in 1:nElem)
ind_t_10 = length(t)
ind_t_08 = round(Int,0.8*ind_t_10)
u1_08 = vcat([vcat(states[ind_t_08][i].nodal.u1[1], states[ind_t_08][i].nodal.u2[1]) for i in 1:nElem]...)
u1_10 = vcat([vcat(states[ind_t_10][i].nodal.u1[1], states[ind_t_10][i].nodal.u2[1]) for i in 1:nElem]...)

# Load reference solution (by Reddy)
x1_ref = collect(LinRange(0,L,21))
u1_08_ref = vec(readdlm(pkgdir(AeroBeams,"test","referenceData","axialTractionCantilever", "t08.txt")))
u1_10_ref = vec(readdlm(pkgdir(AeroBeams,"test","referenceData","axialTractionCantilever", "t10.txt")))

## Plots
# Set paths
relPath = "/test/outputs/figures/axialTractionCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# # Animation
# anim = plot_dynamic_deformation(problem,scale=1e4,plotFrequency=10,plotLimits=([0,2*L],[-L,L],[-L,L]),view=(30,30),lw=5,save=true,savePath=string(relPath,"/axialTractionCantilever_deformation.gif"),displayProgress=true)
# display(anim)

# Plot configurations
lw = 2
ms = 5
msw = 0
colors = [:blue,:orange]
labels = ["\$t\$ = 0.8 s" "\$t\$ = 1.0 s"]
gr()

# Axial displacement at time instances 0.8 s and 1.0 s
plt1 = plot(xlabel="\$x_1\$ [m]", ylabel="\$u_1\$ [mm]")
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, msw=msw, label="Reddy (2005)")
for i=1:2
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(x1,u1_08*1e3, c=colors[1], lw=lw, label=false)
plot!(x1,u1_10*1e3, c=colors[2], lw=lw, label=false)
scatter!(x1_ref,u1_08_ref, c=colors[1], ms=ms, msw=msw, label=false)
scatter!(x1_ref,u1_10_ref, c=colors[2], ms=ms, msw=msw, label=false)
display(plt1)
savefig(string(absPath,"/axialTractionCantilever_u1.pdf"))

println("Finished axialTractionCantilever.jl")
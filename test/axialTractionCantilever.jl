using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles
 
# Beam
L = 0.5
EA,GA,GJ,EIy,EIz = 2e10*1e-6,1e12,1e12,1e12,1e12
ρA,ρI = 8e9*1e-6,0
nElem = 20
stiffnessMatrix = diagm([EA,GA,GA,GJ,EIy,EIz])
inertiaMatrix = diagm([ρA,ρA,ρA,ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])

# Time variables
tf = 1.0
Δt = 2e-3

# BCs
F₀ = 1.0
F = t -> ifelse.(t.>=Δt, F₀, 0.0)
axialForce = create_BC(name="axialForce",beam=beam,node=nElem+1,types=["F1A"],values=[t->F(t)])
clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
axialTractionCantilever = create_Model(name="axialTractionCantilever",beams=[beam],BCs=[clamp,axialForce])

# Create and solve the problem
problem = create_DynamicProblem(model=axialTractionCantilever,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
x1 = [axialTractionCantilever.r_n[i][1] for i in 1:nElem+1]
ind_t_10 = length(problem.timeVector)
ind_t_08 = round(Int,0.8*ind_t_10)
u1_08 = [problem.nodalStatesOverTime[ind_t_08][i].u_n2[1] for i in 1:nElem]
u1_10 = [problem.nodalStatesOverTime[ind_t_10][i].u_n2[1] for i in 1:nElem]
insert!(u1_08, 1, 0)
insert!(u1_10, 1, 0)

# Reference solution by Reddy
x1_ref = collect(LinRange(0,L,21))
u1_08_ref = vec(readdlm(string(pwd(),"/test/referenceData/axialTractionCantilever/t0.8.txt")))
u1_10_ref = vec(readdlm(string(pwd(),"/test/referenceData/axialTractionCantilever/t1.0.txt")))

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 5
relPath = "/test/outputs/figures/axialTractionCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Animation
plot_dynamic_deformation(problem,scale=1e4,plotLimits=[(0,2*L),(-0.1,0.1),(0,1)],save=true,savePath=string(relPath,"/axialTractionCantilever_deformation.gif"))
# Axial displacement at time instances 0.8 s and 1.0 s
gr()
plt1 = plot(xlabel="\$x_1\$ [m]", ylabel="\$u_1\$ [mm]")
plot!(x1,u1_08*1e3, c=:black, lw=lw, label="AeroBeams")
scatter!(x1_ref,u1_08_ref, c=:black, ms=ms, msw=0, label="Reddy (2005)")
plot!([NaN], [NaN], lc=:black, m=:black, lw=2, ms=ms, msw=0, label="\$t\$ = 0.8 s")
plot!(x1,u1_10*1e3, c=:blue, lw=lw, label=false)
scatter!(x1_ref,u1_10_ref, c=:blue, ms=ms, msw=0, label=false)
plot!([NaN], [NaN], lc=:blue, m=:blue, lw=lw, ms=ms, msw=0, label="\$t\$ = 1.0 s")
display(plt1)
savefig(string(absPath,"/axialTractionCantilever_u1.pdf"))

println("Finished axialTractionCantilever.jl")
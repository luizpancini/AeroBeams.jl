using AeroBeams

# Beam 
L = 1
EA = 1
nElem = 2
midElem = div(nElem,2)
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA)])

# Spring (attached at middle and tip nodes)
k = 2
spring = create_Spring(ku=k)
add_spring_to_beams!(beams=[beam,beam],spring=spring,nodesLocalIDs=[midElem+1,nElem+1])

# BCs
u₀ = 1/3
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipDisp = create_BC(name="tipDisp",beam=beam,node=nElem+1,conditions=Dict(:u1A=>u₀))

# Model
axialDoublyAttachedSpringCantilever2 = create_Model(name="axialDoublyAttachedSpringCantilever2",beams=[beam],BCs=[clamp,tipDisp])

# Create and solve problem
problem = create_SteadyProblem(model=axialDoublyAttachedSpringCantilever2)
solve!(problem)

# Get outputs
hist = problem.loadstepHistory
x1 = reduce(vcat,[beam.elements[e].geometry.x1_n1,
                  beam.elements[e].geometry.x1_n2] for e in 1:nElem)
u1 = reduce(vcat,[hist.states[end][e].nodal.u1[1],
                  hist.states[end][e].nodal.u2[1]] for e in 1:nElem)
F1 = reduce(vcat,[hist.states[end][e].nodal.F1[1],
                  hist.states[end][e].nodal.F2[1]] for e in 1:nElem)

# Analytical solution
FrootAnalytical = u₀/L / (1/(2EA) + 1/(2EA+k*L))
FsAnalytical = k * (FrootAnalytical/(1+k*L/(2EA)) * L/(2EA))

# Relative errors
ϵFr = F1[1]/FrootAnalytical - 1
ϵFs = norm(spring.state.Fs)/FsAnalytical - 1

# Display results
println("Relative error in root reaction = $(round(100*ϵFr, digits=2))%")
println("Relative error in spring force = $(round(100*ϵFs, digits=2))%")

## Plots
# Set paths
relPath = "/test/outputs/figures/axialDoublyAttachedSpringCantilever2"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,plotUndeformed=true,plotDistLoads=true,save=true,savePath=string(relPath,"/axialDoublyAttachedSpringCantilever2_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Axial displacement
plt_u1 = plot(x1/L, u1/u₀, lw=lw, label=false, xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$u_1/u_0\$")
display(plt_u1)
savefig(string(absPath,"/axialDoublyAttachedSpringCantilever2_u1.pdf"))

# Axial force
plt_F1 = plot(x1/L, F1/F1[1], lw=lw, label=false, xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$F_1^{\\star}(x_1)/F_1^{\\star}(0)\$")
display(plt_F1)
savefig(string(absPath,"/axialDoublyAttachedSpringCantilever2_F1.pdf"))

println("Finished axialDoublyAttachedSpringCantilever2.jl")
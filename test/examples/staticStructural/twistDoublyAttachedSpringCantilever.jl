using AeroBeams

# Beam 
L = 1
GJ = 1
nElem = 2
midElem = div(nElem,2)
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(GJ=GJ)])

# Spring (attached at middle and tip nodes)
k = 2
spring = create_Spring(kp=k)
add_spring_to_beams!(beams=[beam,beam],spring=spring,nodesLocalIDs=[midElem+1,nElem+1])

# BCs
M₀ = 1e-1
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
moment = create_BC(name="moment",beam=beam,node=nElem+1,conditions=Dict(:M1A=>M₀))

# Model
twistDoublyAttachedSpringCantilever = create_Model(name="twistDoublyAttachedSpringCantilever",beams=[beam],BCs=[clamp,moment])

# Create and solve problem
problem = create_SteadyProblem(model=twistDoublyAttachedSpringCantilever)
@time solve!(problem)
# @profview_allocs solve!(problem) sample_rate=1

# Get outputs
hist = problem.loadstepHistory
x1 = reduce(vcat,[beam.elements[e].geometry.x1_n1,
                  beam.elements[e].geometry.x1_n2] for e in 1:nElem)
p1 = reduce(vcat,[hist.states[end][e].nodal.p1[1],
                  hist.states[end][e].nodal.p2[1]] for e in 1:nElem)
M1 = reduce(vcat,[hist.states[end][e].nodal.M1[1],
                  hist.states[end][e].nodal.M2[1]] for e in 1:nElem)

# Analytical solution
pTipAnalytical = M₀*L * (1/(2GJ) + 1/(2GJ+k*L))
MsAnalytical = M₀*(1-1/(1+k*L/(2GJ)))

# Relative differences
ϵp = p1[end]/pTipAnalytical - 1
ϵMs = norm(spring.state.Ms)/MsAnalytical - 1

# Display results
println("Relative difference in tip twist rotation = $(round(100ϵp,digits=4)) %")
println("Relative difference in spring moment = $(round(100ϵMs,digits=4)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/twistDoublyAttachedSpringCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,plotUndeformed=true,plotDistLoads=true,save=true,savePath=string(relPath,"/twistDoublyAttachedSpringCantilever_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Twist rotation
plt_p1 = plot(x1/L, p1, lw=lw, label=false, xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$p_1\$")
display(plt_p1)
savefig(string(absPath,"/twistDoublyAttachedSpringCantilever_p1.pdf"))

# Twisting moment
plt_M1 = plot(x1/L, M1/M₀, lw=lw, label=false, xlims=[0,1], ylims=[0,1], xlabel="\$x_1/L\$", ylabel="\$M_1^{\\star}/M_0\$")
display(plt_M1)
savefig(string(absPath,"/twistDoublyAttachedSpringCantilever_M1.pdf"))

println("Finished twistDoublyAttachedSpringCantilever.jl")
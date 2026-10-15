using AeroBeams

# Beam
L = 100
E = 420e6
A,Iy = 1, 1/12
∞ = 1e14
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,EA=E*A,EIy=E*Iy)])

# BCs
F = 130e3
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipFollowerForce = create_BC(name="tipFollowerForce",beam=beam,node=nElem+1,conditions=Dict(:Ff3A=>F))

# Model
tipFollowerForceCantilever = create_Model(name="tipFollowerForceCantilever",beams=[beam],BCs=[clamp,tipFollowerForce],units=create_UnitsSystem(length="in",force="lbf"))

# Set system solver options
σ0 = 0.0
σstep = 0.02
NR = create_NewtonRaphson(initLoadFac=σ0,maxLoadFacStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=tipFollowerForceCantilever,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
hist = problem.loadstepHistory
σVec = hist.σ
states = hist.states
tip_u1 = [states[i][nElem].nodal.u2[1] for i in 1:length(σVec)]
tip_u3 = [states[i][nElem].nodal.u2[3] for i in 1:length(σVec)]

# Load reference solution
u1Ref = readdlm(pkgdir(AeroBeams,"test","referenceData","tipFollowerForceCantilever", "u1.txt"))
u3Ref = readdlm(pkgdir(AeroBeams,"test","referenceData","tipFollowerForceCantilever", "u3.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/tipFollowerForceCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,legendPos=:bottomright,save=true,savePath=string(relPath,"/tipFollowerForceCantilever_deformation.pdf"))
display(deformationPlot)

# Plot configurations
colors = [:blue,:green]
labels = ["\$-u_1/L\$" "\$u_3/L\$"]
lw = 2
ms = 3
msw = 0
gr()

# Plot normalized tip displacements over load steps
plt1 = plot(xlabel="\$F\$ [kip]", ylabel="\$-u_1/L, u_3/L\$", title="Tip generalized displacements", xticks=collect(0:10:F), yticks=collect(-0.6:0.2:1.2))
plot!([NaN], [NaN], lc=:black,  lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Simo & Vu-Quoc (1986)")
for i=1:2
    plot!([NaN], [NaN], lc=colors[i], m=colors[i],  lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(σVec*F/(1e3), [-tip_u1/L, tip_u3/L], palette=colors,  lw=lw, label=false)
scatter!([u1Ref[1,:],u3Ref[1,:]], [u1Ref[2,:]/L,u3Ref[2,:]/L], palette=colors, ms=ms, msw=msw, label=false)
display(plt1)
savefig(string(absPath,"/tipFollowerForceCantilever_disp.pdf"))

println("Finished tipFollowerForceCantilever.jl")
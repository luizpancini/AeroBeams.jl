using AeroBeams

# Beam 
R,θ = 0.5, π
L = R*θ
A,Iy = 1e-4, 0.5e-8
E = 72e9
∞ = 1e12
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=E*A,EIy=E*Iy,∞=∞)],p0=transform_to_WM([0;-π/2;0],from=:E321),k=[0;1/R;0])

# BCs
F = 2.5e3
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipFollowerForce = create_BC(name="tipFollowerForce",beam=beam,node=nElem+1,conditions=Dict(:Ff3A=>F))

# Model
tangentiallyForcedArch = create_Model(name="tangentiallyForcedArch",beams=[beam],BCs=[clamp,tipFollowerForce])

# Set system solver options
σ0 = 0
σstep = 0.02
NR = create_NewtonRaphson(initLoadFac=σ0,maxLoadFacStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=tangentiallyForcedArch,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
hist = problem.loadstepHistory
σVec = hist.σ
states = hist.states
tip_u1 = [states[i][nElem].nodal.u2[1] for i in 1:length(σVec)]
tip_u3 = [states[i][nElem].nodal.u2[3] for i in 1:length(σVec)]
tip_angle = [states[i][nElem].nodal.θ[2] for i in 1:length(σVec)]

# Load reference solution
u1Ref = readdlm(pkgdir(AeroBeams,"test","referenceData","tangentiallyForcedArch", "u1.txt"))
u3Ref = readdlm(pkgdir(AeroBeams,"test","referenceData","tangentiallyForcedArch", "u3.txt"))
θRef = readdlm(pkgdir(AeroBeams,"test","referenceData","tangentiallyForcedArch", "theta.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/tangentiallyForcedArch"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/tangentiallyForcedArch_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
ms = 4
msw = 0
labels = ["\$-u_1/R\$" "\$u_3/R\$" "\$-\\theta/(\\pi/2)\$"]
XLabel = "\$-u_1/R, u_3/R, -\\theta/(\\pi/2)\$"
colors = [:blue,:orange,:green]
gr()

# Plot normalized displacements over load steps
plt1 = plot(xlabel=XLabel, ylabel="\$F\$ [kN]", title="Tip generalized displacements",legend=:bottomright)
plot!([NaN], [NaN], lc=:black,  lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Argyris & Symeonidis (1981)")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-tip_u1/R, tip_u3/R, -tip_angle/(π/2)], colors)
    plot!(x, σVec*abs(F)/(1e3), c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1Ref[1,:], u3Ref[1,:], θRef[1,:]], [u1Ref[2,:]/1e3, u3Ref[2,:]/1e3, θRef[2,:]/1e3], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/tangentiallyForcedArch_disp.pdf"))

println("Finished tangentiallyForcedArch.jl")
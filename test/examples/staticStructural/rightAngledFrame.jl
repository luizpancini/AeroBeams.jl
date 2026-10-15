using AeroBeams

# Beam frame
L = 0.24
E, A,Iy = 7.124e10, 0.18e-4, 0.135e-8
∞ = 1e12
S = isotropic_stiffness_matrix(EA=E*A,EIy=E*Iy,∞=∞)
nElem = 20
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,S=[S])
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,S=[S],p0=transform_to_WM([0;π/2;0],from=:E321))

# BCs
F = 5e3
clamp = create_BC(name="clamp",beam=beam1,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipFollowerForce = create_BC(name="tipFollowerForce",beam=beam2,node=nElem+1,conditions=Dict(:Ff3b=>-F))

# Model
rightAngledFrame = create_Model(name="rightAngledFrame",beams=[beam1,beam2],BCs=[clamp,tipFollowerForce])

# Set system solver options
σ0 = 0.0
σstep = 0.02
NR = create_NewtonRaphson(initLoadFac=σ0,maxLoadFacStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=rightAngledFrame,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
hist = problem.loadstepHistory
σVec = hist.σ
states = hist.states
tip_u1 = [states[i][end].nodal.u2[1] for i in 1:length(σVec)]
tip_u3 = [states[i][end].nodal.u2[3] for i in 1:length(σVec)]
tip_angle = [states[i][end].nodal.θ[2] for i in 1:length(σVec)]

# Load reference solution
u1Ref = readdlm(pkgdir(AeroBeams,"test","referenceData","rightAngledFrame", "u1.txt"))
u3Ref = readdlm(pkgdir(AeroBeams,"test","referenceData","rightAngledFrame", "u3.txt"))
θRef = readdlm(pkgdir(AeroBeams,"test","referenceData","rightAngledFrame", "theta.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/rightAngledFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,legendPos=:bottomleft,savePath=string(relPath,"/rightAngledFrame_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
ms = 4
msw = 0
labels = ["\$-u_1/L\$" "\$u_3/L\$" "\$-\\theta/\\pi\$"]
colors = [:blue,:orange,:green]
gr()

# Plot normalized displacements over load steps
plt1 = plot(xlabel="\$-u_1/L, u_3/L, -\\theta/\\pi\$", ylabel="\$F\$ [kN]", title="Tip generalized displacements")
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Argyris & Symeonidis (1981)")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-tip_u1/L, tip_u3/L, tip_angle/π], colors)
    plot!(x, σVec*F/1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1Ref[1,:], u3Ref[1,:], θRef[1,:]], [u1Ref[2,:], u3Ref[2,:], θRef[2,:]], colors)
    scatter!(x, y/1e3, c=c, ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/rightAngledFrame_disp.pdf"))

println("Finished rightAngledFrame.jl")
using AeroBeams

# Beam frame
L = 120
A = 6
Iy = 2
E = 7.2e6
∞ = 1e10
S = isotropic_stiffness_matrix(EA=E*A,EIy=E*Iy,∞=∞)
nElem = 20
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,S=[S],p0=transform_to_WM([0;-π/2;0],from=:E321))
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,S=[S])

# BCs
F = 40e3
elemForce = div(nElem,5)
support1 = create_BC(name="support1",beam=beam1,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
support2 = create_BC(name="support2",beam=beam2,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
force = create_BC(name="force",beam=beam2,node=elemForce+1,conditions=Dict(:F3A=>-F))

# Model
LeeFrameDeadLoad = create_Model(name="LeeFrameDeadLoad",beams=[beam1,beam2],BCs=[support1,support2,force],units=create_UnitsSystem(length="in",force="lbf"))

# Set system solver options (arclength continuation is required to trace the snap-through path)
systemSolver = create_ArclengthSolver(maxArcLength=5)

# Create and solve the problem
problem = create_SteadyProblem(model=LeeFrameDeadLoad,systemSolver=systemSolver)
solve!(problem)

# Get solution at load steps
hist = problem.loadstepHistory
σVec = hist.σ
states = hist.states
u1_atForce = [states[i][nElem+elemForce].nodal.u2[1] for i in 1:length(σVec)]
u3_atForce = [states[i][nElem+elemForce].nodal.u2[3] for i in 1:length(σVec)]

# Load reference solution
u1Ref = readdlm(pkgdir(AeroBeams,"test","referenceData","LeeFrameDeadLoad", "u1.txt"))
u3Ref = readdlm(pkgdir(AeroBeams,"test","referenceData","LeeFrameDeadLoad", "u3.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/LeeFrameDeadLoad"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,legendPos=:bottomright,save=true,savePath=string(relPath,"/LeeFrameDeadLoad_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
ms = 3
msw = 0
x = [u1_atForce/L, -u3_atForce/L]
labels = ["\$u_1/L\$" "\$-u_3/L\$"]
colors = [:blue,:orange]
gr()

# Plot normalized displacements over load steps
plt1 = plot(xlabel="\$u_1/L, -u_3/L,\$", ylabel="\$F\$ [kip]", ylims=[-10,40], title="Displacements at point of force application")
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Simo and Vu-Quoc (1986)")
for i=1:2
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(x, σVec*F/(1e3), lw=lw,palette=colors,label=false)
scatter!([u1Ref[1,:],u3Ref[1,:]]/L, [u1Ref[2,:],u3Ref[2,:]], palette=colors,ms=ms,msw=msw,label=false)
display(plt1)
savefig(string(absPath,"/LeeFrameDeadLoad_disp.pdf"))

println("Finished LeeFrameDeadLoad.jl")
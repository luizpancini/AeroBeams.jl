using AeroBeams

# Beam 
R,θ = 100, deg2rad(215)
L = R*θ
E,A,Iy = 1.5e6, 2, 2/3
∞ = 1e10
nElem = 80
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=E*A,EIy=E*Iy,∞=∞)],p0=transform_to_WM([0;-θ/2;0],from=:E321),k=[0;1/R;0])

# BCs
λ = 9.0
elemForce = div(nElem,2)
pin = create_BC(name="pin",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
clamp = create_BC(name="clamp",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
force = create_BC(name="force",beam=beam,node=elemForce+1,conditions=Dict(:F3A=>-λ*E*Iy/R^2))

# Model
pinnedClampedArch = create_Model(name="pinnedClampedArch",beams=[beam],BCs=[pin,clamp,force],units=create_UnitsSystem(length="in",force="lbf"))

# Set system solver options (arclength continuation is required to trace the snap-through path)
systemSolver = create_ArclengthSolver(maxArcLength=5)

# Create and solve the problem
problem = create_SteadyProblem(model=pinnedClampedArch,systemSolver=systemSolver)
solve!(problem)

# Get solution at partial load steps
hist = problem.loadstepHistory
σVec = hist.σ
states = hist.states
u1_atForce = [states[i][elemForce].nodal.u2[1] for i in 1:length(σVec)]
u3_atForce = [states[i][elemForce].nodal.u2[3] for i in 1:length(σVec)]

# Load reference solution
u1Ref = readdlm(pkgdir(AeroBeams,"test","referenceData","pinnedClampedArch", "u1.txt"))
u3Ref = readdlm(pkgdir(AeroBeams,"test","referenceData","pinnedClampedArch", "u3.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/pinnedClampedArch"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/pinnedClampedArch_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
ms = 4
msw = 0
x = [-u1_atForce/R, -u3_atForce/R]
labels = ["\$-u_1/R\$" "\$-u_3/R\$"]
colors = [:blue,:orange]
gr()

# Plot normalized displacements over load steps
plt1 = plot(xlabel="\$-u_1/L, -u_3/L,\$", ylabel="\$\\lambda\$",title="Displacements at point of force application",legend=:bottomright)
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="DaDeppo & Schmidt (1975)")
for i=1:2
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
plot!(x, σVec*λ, lw=lw,palette=colors,label=false)
scatter!([u1Ref[1,:],u3Ref[1,:]], [u1Ref[2,:],u3Ref[2,:]], palette=colors,ms=ms, msw=msw,label=false)
display(plt1)
savefig(string(absPath,"/pinnedClampedArch_disp.pdf"))

println("Finished pinnedClampedArch.jl")
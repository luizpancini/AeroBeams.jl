using AeroBeams

# Beam
L = 1
E = 210e6
A,Iy = 20e-4, 5/3*1e-8
∞ = 1e12
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=E*A,EIy=E*Iy,∞=∞)])

# BCs
q = 100
add_distributed_loads_to_beam!(beam=beam,loads=Dict(:ff_b_of_x1t => (x1,t)->[0; 0; q]))
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
distributedLoadCantilever = create_Model(name="distributedLoadCantilever",beams=[beam],BCs=[clamp],units=create_UnitsSystem(length="in",force="lbf"))

# Set system solver options
σ0 = 0
σstep = 0.02
NR = create_NewtonRaphson(initLoadFac=σ0,maxLoadFacStep=σstep)

# Create and solve the problem
problem = create_SteadyProblem(model=distributedLoadCantilever,systemSolver=NR)
solve!(problem)

# Get solution at partial load steps
hist = problem.loadstepHistory
σVec = hist.σ
states = hist.states
tip_u1 = [states[i][nElem].nodal.u2[1] for i in 1:length(σVec)]
tip_u3 = [states[i][nElem].nodal.u2[3] for i in 1:length(σVec)]
tip_angle = [states[i][nElem].nodal.θ[2] for i in 1:length(σVec)]

# Load reference solution
u1Ref = readdlm(pkgdir(AeroBeams,"test","referenceData","distributedLoadCantilever", "u1.txt"))
u3Ref = readdlm(pkgdir(AeroBeams,"test","referenceData","distributedLoadCantilever", "u3.txt"))
θRef = readdlm(pkgdir(AeroBeams,"test","referenceData","distributedLoadCantilever", "theta.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/distributedLoadCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,plotLimits=([-0.5,1],[-1,1],[-0.2,0.6]),save=true,savePath=string(relPath,"/distributedLoadCantilever_deformation.pdf"))
display(deformationPlot)

# Plot configurations
ts = 10
fs = 16
lw = 2
ms = 4
msw = 0
labels = ["\$-u_1/L\$" "\$u_3/L\$" "\$-\\theta/\\pi\$"]
colors = [:blue,:orange,:green]
gr()

# Plot normalized displacements over load steps
plt1 = plot(xlabel="Tip generalized displacements", ylabel="\$q\$ [kip]",legend=:bottomright, tickfont=font(ts), guidefont=font(fs), legendfontsize=8)
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Argyris & Symeonidis (1981)")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-tip_u1/L, tip_u3/L, -tip_angle/π], colors)
    plot!(x, σVec*q, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1Ref[1,:], u3Ref[1,:], θRef[1,:]], [u1Ref[2,:], u3Ref[2,:], θRef[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/distributedLoadCantilever_disp.pdf"))

println("Finished distributedLoadCantilever.jl")
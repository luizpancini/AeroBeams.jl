using AeroBeams

# Beam
L = 1
EIy = 1
nElem = 40
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EIy)])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipMoment = create_BC(name="tipMoment",beam=beam,node=nElem+1,conditions=Dict(:M2A=>2*π*EIy/L))

# Model
tipMomentCantilever = create_Model(name="tipMomentCantilever",beams=[beam],BCs=[clamp,tipMoment],units=create_UnitsSystem(length="in",force="lbf"))

# Set system NLsol options
σ0 = 0.0
σstep = 0.02
NLsol = create_NLsolver(method=:newton,initLoadFac=σ0,maxLoadFacStep=σstep,displayStatus=false)

# # Create and solve the problem
# problem = create_SteadyProblem(model=tipMomentCantilever,systemSolver=NLsol)
# solve!(problem)

# println("Solving with allocation profiler:")
# problem = create_SteadyProblem(model=tipMomentCantilever,systemSolver=NLsol)
# @profview_allocs solve!(problem) sample_rate=0.1

# println("Solving with time profiler:")
# problem = create_SteadyProblem(model=tipMomentCantilever,systemSolver=NLsol)
# using Profile
# Profile.clear()
# Profile.init(n=10^8, delay=1e-4)
# @profview solve!(problem)

println("Solving problem with time measurement:")
problem = create_SteadyProblem(model=tipMomentCantilever,systemSolver=NLsol)
@time solve!(problem)

# Get solution at partial load steps
hist = problem.loadstepHistory
σVec = hist.σ
states = hist.states
tip_u1 = [states[i][nElem].nodal.u2[1] for i in 1:length(σVec)]
tip_u3 = [states[i][nElem].nodal.u2[3] for i in 1:length(σVec)]
tip_angle = [states[i][nElem].nodal.θ[2] for i in 1:length(σVec)]

## Plots
# Set paths
relPath = "/test/outputs/figures/tipMomentCantilever"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/tipMomentCantilever_deformation.pdf"))
display(deformationPlot)

# Plot internal bending moment and strain
plot_steady_outputs(problem,outputs=["M2","κ2"],save=true,saveFolder=string(relPath,"/"))

# Plot configurations
lw = 2
y = [1.0 .+ tip_u1/L, tip_u3/L, tip_angle/π]
labels = ["\$1+u_1/L\$" "\$u_3/L\$" "\$\\theta/\\pi\$"]
colors = [:blue,:orange,:green]
gr()

# Plot normalized displacements over load steps
plt1 = plot(xlabel="\$ML/(2\\pi EI)\$", ylabel="\$1+u_1/L, u_3/L, \\theta/L\$", title="Tip generalized displacements")
plot!(σVec, y, palette=colors, lw=lw, label=false)
halfNσ = round(Int,length(σVec)/2)
for i=1:3
    if i==1
        isfinite(y[i][halfNσ]) && annotate!(σVec[halfNσ], y[i][halfNσ], text(labels[i], :bottom, :left, colors[i]))
    else
        isfinite(y[i][halfNσ]) && annotate!(σVec[halfNσ], y[i][halfNσ], text(labels[i], :bottom, :right, colors[i]))
    end
end
display(plt1)
savefig(string(absPath,"/tipMomentCantilever_disp.pdf"))

println("Finished tipMomentCantilever.jl")
using AeroBeams

# Beam
R,θ = 2.54, deg2rad(120)
k2 = 1/R
L = R*θ
A,Iy = 4.05e-4, 13.1e-8
E = 70.4e9
EA,EIy = E*A,E*Iy
nElem = 80
beam = create_Beam(name="arch",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e12,EA=EA,EIy=EIy)],p0=transform_to_WM([0;-θ/2;0]; from=:E321),k=[0;k2;0]);

# BCs
λ = 8.7
q = -λ*EIy/R^2
add_distributed_loads_to_beam!(beam=beam,loads=Dict(:f_A_of_x1t => (x1,t)->[0; 0; q]))
clamp1 = create_BC(name="clamp1",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
clamp2 = create_BC(name="clamp2",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
archUnderDeadPressure = create_Model(name="archUnderDeadPressure",beams=[beam],BCs=[clamp1,clamp2])

# System solver
σ0 = 0
σstep = 0.02
NR = create_NewtonRaphson(initLoadFac=σ0,maxLoadFacStep=σstep)

# Create and solve the steady problem
problem = create_SteadyProblem(model=archUnderDeadPressure,systemSolver=NR)
solve!(problem)

# Outputs
hist = problem.loadstepHistory
σVec = hist.σ
states = hist.states
mid_u3 = [states[i][div(nElem,2)].nodal.u2[3] for i in eachindex(σVec)]

## Plots
# Set paths
relPath = "/test/outputs/figures/archUnderDeadPressure"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/archUnderDeadPressure_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Plot normalized displacements over load steps
plt1 = plot(-mid_u3/R, σVec*λ, color=:black, lw=lw, xlabel="Midpoint \$-u_3/R\$", ylabel="\$\\lambda\$", label=false)
display(plt1)
savefig(string(absPath,"/archUnderDeadPressure_disp.pdf"))

println("Finished archUnderDeadPressure.jl") #src
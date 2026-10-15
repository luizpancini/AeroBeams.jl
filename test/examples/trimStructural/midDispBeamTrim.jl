using AeroBeams

# Beam
L = 1
EI = 1e3/48
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EI)])

# BCs
F = 1
u3_mid = -F*L/(48*EI)
Rguess = 0
midDisp = create_BC(name="midDisp",beam=beam,node=div(nElem,2)+1,conditions=Dict(:u3A=>u3_mid))
pin = create_BC(name="pin",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
roller = create_BC(name="roller",beam=beam,node=nElem+1,conditions=Dict(:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
rollerReaction = create_BC(name="rollerReaction",beam=beam,node=nElem+1,conditions=Dict(:F3A=>Rguess),toBeTrimmed=trues(1))

# Model
midDispBeamTrim = create_Model(name="midDispBeamTrim",beams=[beam],BCs=[pin,roller,midDisp,rollerReaction])

# Set up a nonlinear system solver
solver = create_NewtonRaphson(displayStatus=false,trimSolverMethod=:LM)

# Create and solve the trim problem
problem = create_TrimProblem(model=midDispBeamTrim,systemSolver=solver,getInertiaMatrix=false)
solve!(problem)

# # Create and solve the problem
# println("Solving with time and memory measurement...")
# problem = create_TrimProblem(model=midDispBeamTrim,systemSolver=solver,getInertiaMatrix=false)
# @time solve!(problem)

# println("Solving with memory profiling...")
# problem = create_TrimProblem(model=midDispBeamTrim,systemSolver=solver,getInertiaMatrix=false)
# @profview_allocs solve!(problem) sample_rate=1

# println("Solving with time profiling...")
# problem = create_TrimProblem(model=midDispBeamTrim,systemSolver=solver,getInertiaMatrix=false)
# using Profile
# Profile.clear()
# Profile.init(n=10^7,delay=1e-6)
# @profview begin
#     for i in 1:1000
#         solve!(problem)
#     end
# end

# Get solution
hist = problem.loadstepHistory
x1 = reduce(vcat,[beam.elements[e].geometry.x1_n1,
                  beam.elements[e].geometry.x1_n2] for e in 1:nElem)
u3 = reduce(vcat,[hist.states[end][e].nodal.u1[3],
                  hist.states[end][e].nodal.u2[3]] for e in 1:nElem)                  
F3 = reduce(vcat,[hist.states[end][e].nodal.F1[3],
                  hist.states[end][e].nodal.F2[3]] for e in 1:nElem)
M2 = reduce(vcat,[hist.states[end][e].nodal.M1[2],
                  hist.states[end][e].nodal.M2[2]] for e in 1:nElem)
roller_F3 = problem.trimSolution.trimLoads[1]

# Compare to analytical solution (valid for small displacements)
roller_F3_analytical = F/2
ϵ_rel = roller_F3/roller_F3_analytical - 1
println("Relative errors:\nF3: $(round(100ϵ_rel, digits=4)) %")

## Plots
# Set paths
relPath = "/test/outputs/figures/midLoadedBeamTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed state
deformationPlot = plot_steady_deformation(problem,scale=1e2,scalePos=[0.2;0.2],showScale=true,save=true,savePath=string(relPath,"/midLoadedBeamTrim_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# u3
plt1 = plot(x1/L, u3/(F*L/(48*EI)), lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$u_3 / (FL/48EI)\$")
display(plt1)
savefig(string(absPath,"/midLoadedBeamTrim_u3.pdf"))

# F3
plt2 = plot(x1/L, F3/(F/2), lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$F_3/(F/2)\$ [N]")
display(plt2)
savefig(string(absPath,"/midLoadedBeamTrim_F3.pdf"))

# M2
plt3 = plot(x1/L, M2/(F*L/4), lw=lw, label=false, xlabel="\$x_1/L\$", ylabel="\$M_2/(FL/4)\$")
display(plt3)
savefig(string(absPath,"/midLoadedBeamTrim_M2.pdf"))

println("Finished midDispBeamTrim.jl")
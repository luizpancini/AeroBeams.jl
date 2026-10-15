using AeroBeams

# Beam
L = 1
EI = 1000/3
nElem = 20
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EIy=EI)])

# BCs
F = 1
trimF3Guess = 0
trimM2Guess = 0
tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,conditions=Dict(:F3A=>-F))
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
clampReactions = create_BC(name="clampReactions",beam=beam,node=1,conditions=Dict(:F3A=>trimF3Guess,:M2A=>trimM2Guess),toBeTrimmed=trues(2))

# Model
tipLoadedCantileverTrim = create_Model(name="tipLoadedCantileverTrim",beams=[beam],BCs=[clamp,tipForce,clampReactions])

# Set NR system solver
solver = create_NewtonRaphson(displayStatus=true)

# Create and solve the problem
problem = create_TrimProblem(model=tipLoadedCantileverTrim,systemSolver=solver,getInertiaMatrix=false)
solve!(problem)

# Get numerical solution
hist = problem.loadstepHistory
x1 = reduce(vcat,[beam.elements[e].geometry.x1_n1,
                  beam.elements[e].geometry.x1_n2] for e in 1:nElem)
u3 = reduce(vcat,[hist.states[end][e].nodal.u1[3],
                  hist.states[end][e].nodal.u2[3]] for e in 1:nElem)                  
F3 = reduce(vcat,[hist.states[end][e].nodal.F1[3],
                  hist.states[end][e].nodal.F2[3]] for e in 1:nElem)
M2 = reduce(vcat,[hist.states[end][e].nodal.M1[2],
                  hist.states[end][e].nodal.M2[2]] for e in 1:nElem)

# Analytical solution (valid for small displacements)
u3_analytical(x1) = -F/(6*EI)*x1.^2 .* (3*L .- x1)
F3_analytical(x1) = -F
M2_analytical(x1) = F*L * (1 .- x1)

## Plots
# Set paths
relPath = "/test/outputs/figures/tipLoadedCantileverTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 4
msw = 0
gr()

# u3
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$u_3 / (FL/3EI)\$")
plot!(x1/L, u3/(F*L/(3*EI)), c=:black, lw=lw, label="Numerical")
scatter!(x1/L, u3_analytical(x1)/(F*L/(3*EI)), c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt1)
savefig(string(absPath,"/tipLoadedCantileverTrim_u3.pdf"))

# F3
plt2 = plot(xlabel="\$x_1/L\$", ylabel="\$F_3/F\$ [N]", ylims=[-2,0])
plot!(x1/L, F3/F, c=:black, lw=lw, label="Numerical")
scatter!(x1/L, F3_analytical.(x1)/F, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt2)
savefig(string(absPath,"/tipLoadedCantileverTrim_F3.pdf"))

# M2
plt3 = plot(xlabel="\$x_1/L\$", ylabel="\$M_2/(FL)\$")
plot!(x1/L, M2/(F*L), c=:black, lw=2, label="Numerical")
scatter!(x1/L, M2_analytical(x1)/(F*L), c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt3)
savefig(string(absPath,"/tipLoadedCantileverTrim_M2.pdf"))

println("Finished tipLoadedCantileverTrim.jl")
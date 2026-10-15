using AeroBeams

# Altitude [m]
h = 20e3

# Airspeed [m/s]
U = 30

# Stiffness factor (for the structure)
λ = 1

# Twist curvature [1/m]
k1 = 0.002

# Wing and stabilizers parasite drag
wingCd0 = stabsCd0 = 0

# Flag to include vertical stabilizer
includeVS = false

# Flag to clamp wing root
rootClamped = true

# Discretization
nElemWing = 32
nElemTailBoom = 5
nElemHorzStabilizer = 4
nElemVertStabilizer = 2

# Aerodynamic solver
aeroSolver = Indicial()

# Model
cHALEclampedSteady,leftWing,rightWing,_ = create_cHALE(aeroSolver=aeroSolver,stiffnessFactor=λ,altitude=h,airspeed=U,k1=k1,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,wingCd0=wingCd0,stabsCd0=stabsCd0,includeVS=includeVS,rootClamped=rootClamped)

# Create and solve steady problem
problem = create_SteadyProblem(model=cHALEclampedSteady)
solve!(problem)

# Undeformed nodal positions
x1_0_l = reduce(vcat, [leftWing.elements[e].assemblyData.r_n1[1], leftWing.elements[e].assemblyData.r_n2[1]] for e in 1:div(nElemWing,2))
x1_0_r = reduce(vcat, [rightWing.elements[e].assemblyData.r_n1[1], rightWing.elements[e].assemblyData.r_n2[1]] for e in 1:div(nElemWing,2))
x3_0_l = reduce(vcat, [leftWing.elements[e].assemblyData.r_n1[3], leftWing.elements[e].assemblyData.r_n2[3]] for e in 1:div(nElemWing,2))
x3_0_r = reduce(vcat, [rightWing.elements[e].assemblyData.r_n1[3], rightWing.elements[e].assemblyData.r_n2[3]] for e in 1:div(nElemWing,2))
x1_0 = vcat(x1_0_l,x1_0_r)
x3_0 = vcat(x3_0_l,x3_0_r)

# Undeformed elemental positions
x1_e_l = reverse(-[leftWing.elements[e].geometry.x1 for e in 1:div(nElemWing,2)])
x1_e_r = [rightWing.elements[e].geometry.x1 for e in 1:div(nElemWing,2)]
x1_e = vcat(x1_e_l,x1_e_r)

# Displacements over span
hist = problem.loadstepHistory
u1_of_x1 = reduce(vcat,[hist.states[end][e].nodal.u1[1],
                        hist.states[end][e].nodal.u2[1]] for e in 1:nElemWing)
u3_of_x1 = reduce(vcat,[hist.states[end][e].nodal.u1[3],
                        hist.states[end][e].nodal.u2[3]] for e in 1:nElemWing)

# Deformed nodal positions
x1_def = x1_0 .+ u1_of_x1
x3_def = x3_0 .+ u3_of_x1

# Angle of attack over span
α_of_x1 = [rad2deg(hist.aeroVariables[end][e].unsteadyKin.αₑ) for e in 1:nElemWing]

## Plots
# Set paths
relPath = "/test/outputs/figures/cHALEclampedSteady"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem,view=(30,30),save=true,savePath=string(relPath,"/cHALEclampedSteady_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Normalized deformed wingspan
L = maximum(x1_0)
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$ [% semispan]", xlims=[-1,1], ylims=[-10,30])
plot!(x1_def/L, x3_def/L*100, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/cHALEclampedSteady_def.pdf"))

# Angle of attack over wingspan
plt2 = plot(xlabel="\$x_1/L\$", ylabel="\$\\alpha\$ [deg]", xlims=[-1,1])
plot!(x1_e/L, α_of_x1, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/cHALEclampedSteady_aoa.pdf"))

println("Finished cHALEclampedSteady.jl")